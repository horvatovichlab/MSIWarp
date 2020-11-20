#include "ransac.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>  // RANSAC

#include "util/parallel.hpp"
#include "warp_util.hpp"

using std::pair;
using std::vector;

namespace warp::ransac {

template <class Func, class Params>
recalibration_function align_ransac_template(
    const std::vector<peak_pair>& pairs, const params& ransac_params,
    Func node_func, const Params& node_params_inliers) {
  // find
  util::params_uniform node_params_ransac{node_params_inliers.inst,
                                          ransac_params.n_steps,
                                          ransac_params.n_samples,
                                          ransac_params.n_segments + 1,
                                          ransac_params.mz_begin,
                                          ransac_params.mz_end,
                                          ransac_params.slack};

  const auto pairs_inliers =
      ransac_pairs(pairs, ransac_params, node_params_ransac);
  const auto nodes_inliers = node_func(pairs_inliers, node_params_inliers);

  // find the optimal shifts using nodes placed based on inliers
  auto optimal_shifts =
      find_optimal_warping_pairs(pairs_inliers, nodes_inliers);

  return nodes_to_recal(nodes_inliers, optimal_shifts);
}

recalibration_function align_ransac_uniform(
    const std::vector<peak_pair>& pairs, const params& ransac_params,
    const util::params_uniform& node_params) {
  return align_ransac_template(pairs, ransac_params,
                               util::get_warping_nodes_uniform, node_params);
}

recalibration_function align_ransac_density(
    const std::vector<peak_pair>& pairs, const params& ransac_params,
    const util::params_density& node_params) {
  return align_ransac_template(pairs, ransac_params,
                               util::get_warping_nodes_density, node_params);
}

std::vector<peak_pair> ransac_pairs(const std::vector<peak_pair>& pairs,
                                    const params& ransac_params,
                                    const util::params_uniform& node_params) {
  if (pairs.size() < ransac_params.min_matched_peaks) {
    return pairs;  // early return if too few peak matches
  }

  const auto nodes_ransac = util::get_warping_nodes_uniform(pairs, node_params);

  vector<vector<peak_pair>> pairs_ransac;
  for (size_t i = 0; nodes_ransac.size() - 1; ++i) {
    const auto& n_left = nodes_ransac[i];
    const auto& n_right = nodes_ransac[i + 1];

    pairs_ransac.emplace_back(peak_pairs_between(pairs, n_left.mz, n_right.mz));
  }

  const auto rr =
      ransac(pairs_ransac, nodes_ransac, ransac_params.n_iterations,
             ransac_params.n_samples, ransac_params.distance_threshold);

  // use the model with the most inliers as the best
  const auto it_max_inliers = std::max_element(
      rr.inliers.begin(), rr.inliers.end(), [](const auto& v1, const auto& v2) {
        return std::count(v1.begin(), v1.end(), true) <
               std::count(v2.begin(), v2.end(), true);
      });

  // copy and return the pairs that fit the best model
  vector<peak_pair> pairs_inliers;
  pairs_inliers.reserve(pairs.size());

  for (size_t i = 0; i < pairs.size(); ++i) {
    if ((*it_max_inliers)[i]) {
      pairs_inliers.push_back(pairs[i]);
    }
  }

  return pairs_inliers;
}

ransac_result ransac(const vector<vector<peak_pair>>& peak_pairs,
                     const node_vec& warping_nodes, size_t n_iterations,
                     size_t n_samples, double distance_threshold) {
  size_t n_steps = warping_nodes.front().mz_shifts.size();
  size_t n_segments = peak_pairs.size();
  size_t n_peaks =
      std::accumulate(peak_pairs.begin(), peak_pairs.end(), size_t{0},
                      [](size_t sum, const auto& p) { return sum + p.size(); });

  // setup random sampling
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> udist(0, 100000);

  //
  vector<double> errors(n_iterations, 0.0);
  vector<vector<bool>> inliers(n_iterations, vector<bool>(n_peaks, false));
  vector<vector<size_t>> optimal_paths;
  vector<vector<vector<size_t>>> maybe_inliers;

  // pre-allocate warping surfaces
  vector<vector<double>> warping_surfs(n_segments,
                                       vector<double>(n_steps * n_steps, 0.0));

  for (size_t i = 0; i < n_iterations; ++i) {
    vector<vector<size_t>> indices_pairs;

    for (size_t j = 0; j < n_segments; ++j) {
      vector<size_t> indices_j;

      auto& warping_surf_j = warping_surfs[j];
      std::fill(warping_surf_j.begin(), warping_surf_j.end(), 0.0);

      const auto& peak_pairs_segment = peak_pairs[j];

      // sample m peaks per segment
      vector<peak_pair> maybe_inliers_j;
      for (size_t k = 0; k < n_samples; ++k) {
        size_t x = udist(gen) % peak_pairs_segment.size();
        indices_j.push_back(x);
        maybe_inliers_j.push_back(peak_pairs_segment[x]);
      }

      // store indices of maybe_inliers
      indices_pairs.push_back(indices_j);

      detail::compute_warping_surf_impl(warping_surf_j, maybe_inliers_j,
                                        warping_nodes[j], warping_nodes[j + 1]);
    }

    auto optimal_path = detail::optimal_warping_impl(warping_surfs, n_steps);

    // residual error and set of inliers for this iteration
    auto& inliers_i = inliers[i];
    double e_i = 0.0;
    size_t n_inliers_i = 0;

    size_t peak_index = 0;

    // test how many peaks fit warping
    for (size_t j = 0; j < n_segments; ++j) {
      const auto& peak_pairs_segment = peak_pairs[j];
      const node& left = warping_nodes[j];
      const node& right = warping_nodes[j + 1];

      double mz_left_warped = left.mz + left.mz_shifts[optimal_path[j]];
      double mz_right_warped = right.mz + right.mz_shifts[optimal_path[j + 1]];

      for (const auto& p : peak_pairs_segment) {
        double mz_warped = lerp(p.second.mz, left.mz, right.mz, mz_left_warped,
                                mz_right_warped);
        double d = std::pow((p.first.mz - mz_warped) / p.second.sigma, 2);

        if (d < distance_threshold) {
          n_inliers_i++;
          inliers_i[peak_index] = true;
          e_i += d;
        }

        peak_index++;
      }
    }

    errors[i] = std::sqrt(e_i / n_inliers_i);
    optimal_paths.push_back(optimal_path);
    maybe_inliers.push_back(indices_pairs);
  }

  return {errors, inliers, optimal_paths, maybe_inliers};
}

}  // namespace warp::ransac