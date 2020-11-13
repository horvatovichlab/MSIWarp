#include "warp.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>   // RANSAC

#include "util/parallel.hpp"

using std::pair;
using std::vector;

namespace warp {

/* linear interpolation. */
inline double lerp(double x,
                   double x_min,
                   double x_max,
                   double y_min,
                   double y_max) {
  double a = (x - x_min) / (x_max - x_min);
  return y_min + a * (y_max - y_min);
}

/* utility function */
peak_vec interpolate_peaks(const peak_vec& peaks,
                           double mz_begin,
                           double mz_end,
                           double mz_begin_ref,
                           double mz_end_ref) {
  peak_vec warped_peaks;
  warped_peaks.reserve(peaks.size());

  for (const auto& p : peaks) {
    double mz_new = lerp(p.mz, mz_begin, mz_end, mz_begin_ref, mz_end_ref);
    warped_peaks.push_back({p.id, mz_new, p.height, p.sigma_mz});
  }

  return warped_peaks;
}

vector<size_t> find_optimal_warping(const peak_vec& s_r,
                                    const peak_vec& s_s,
                                    const node_vec& nodes,
                                    double epsilon) {
  // Only use overlapping peaks
  const auto s_m = overlapping_peak_pairs(s_r, s_s, epsilon);

  std::vector<std::vector<double>> warping_surfs;
  for (size_t i = 0; i < nodes.size() - 1; ++i) {
    const auto& n_l = nodes[i];
    const auto& n_r = nodes[i + 1];

    const auto s_m_i = peak_pairs_between(s_m, n_l.mz, n_r.mz);
    const auto w_i = compute_warping_surf(s_m_i, n_l, n_r);
    warping_surfs.push_back(w_i);
  }

  size_t n_steps = nodes.front().mz_shifts.size();
  return optimal_warping(warping_surfs, n_steps);
}

vector<vector<size_t>> find_optimal_warpings(const vector<peak_vec>& spectra,
                                     const peak_vec& s_r,
                                     const node_vec& nodes, double epsilon) {
  auto out = par::for_each(
      spectra.size() / 5, spectra.begin(), spectra.end(),
      [&](const auto& s_i) { return find_optimal_warping(s_r, s_i, nodes, epsilon); });

  return out;
}

peak_vec warp_peaks(const peak_vec& peaks,
                    const node_vec& nodes,
                    const std::vector<size_t>& moves) {
    size_t n = nodes.size();
    
    recalibration_function recal_func(n);    
    for (size_t i = 0; i < n; ++i) {
      recal_func[i] = {nodes[i].mz, nodes[i].mz_shifts[moves[i]]};
    }

    return warp_peaks_unique(peaks, recal_func);
}

peak_vec warp_peaks_unique(const peak_vec& peaks,
                           const recalibration_function& recal_func) {
  peak_vec warped_peaks;
  warped_peaks.reserve(peaks.size());

  for (size_t i = 0; i < recal_func.size() - 1; ++i) {
    const auto& left = recal_func[i];
    const auto& right = recal_func[i + 1];

    double mz_left = left.first;
    double mz_right = right.first;
    double mz_delta_left = left.second;
    double mz_delta_right = right.second;

    peak_vec peaks_i = peaks_between(peaks, mz_left, mz_right);
    peak_vec peaks_warped_i =
        interpolate_peaks(peaks_i, mz_left, mz_right, mz_left + mz_delta_left,
                          mz_right + mz_delta_right);

    for (const auto& p_w : peaks_warped_i) {
      warped_peaks.push_back(p_w);
    }
  }

  return warped_peaks;
}

node_vec init_nodes(const vector<double>& mzs,
                    const vector<double>& sigmas,
                    uint32_t n_steps) {
  size_t n_nodes = mzs.size();
  node_vec nodes;
  nodes.reserve(n_nodes);

  // All possible moves (contraints in the warping function can be added here)
  for (size_t i = 0; i < n_nodes; ++i) {
    double mz_i = mzs[i];
    double sigma_i = sigmas[i];

    size_t m = 2 * n_steps + 1;
    vector<double> mz_shifts_i(m, 0.0);

    for (size_t j = 0; j < m; ++j) {
      int k = j - n_steps;
      mz_shifts_i[j] = sigma_i * k;
    }

    nodes.push_back({mz_i, sigma_i, mz_shifts_i});
  }

  return nodes;
}

template <class T, class Comp>
vector<T> get_range(const vector<T>& v, double begin, double end, Comp comp) {
  vector<T> out;
  if (begin >= end) {
    return out;
  }

  auto first = std::lower_bound(v.begin(), v.end(), begin, comp);
  auto last = std::lower_bound(v.begin(), v.end(), end, comp);

  std::copy(first, last, std::back_inserter(out));
  return out;
}

peak_vec peaks_between(const peak_vec& peaks, double mz_begin, double mz_end) {
  auto comp_mz = [](const peak& p_a, double mz) { return p_a.mz < mz; };
  return get_range(peaks, mz_begin, mz_end, comp_mz);
}

vector<peak_pair> peak_pairs_between(const vector<peak_pair>& peak_pairs,
                                     double mz_begin,
                                     double mz_end) {
  auto comp_mz = [](const peak_pair& p, double mz) { return p.second.mz < mz; };
  return get_range(peak_pairs, mz_begin, mz_end, comp_mz);
}

// TODO:
peak_vec peaks_top_n(const peak_vec& peaks, size_t n) {
  size_t max_n = std::min(n, peaks.size());

  // Sort peaks by height in descending order then take the first n peaks.
  peak_vec out = peaks;
  std::sort(out.begin(), out.end(),
            [](const peak& a, const peak& b) { return a.height > b.height; });
  out.erase(out.begin() + max_n, out.end());

  return out;
}

// TODO:
vector<peak_pair> peak_pairs_top_n(const vector<peak_pair>& peak_pairs,
                                   size_t n) {
  size_t max_n = std::min(n, peak_pairs.size());

  auto peak_pairs_sorted = peak_pairs;
  std::sort(peak_pairs_sorted.begin(), peak_pairs_sorted.end(),
            [](const peak_pair& p_a, const peak_pair& p_b) {
              return (p_a.first.height * p_a.second.height) >
                     (p_b.first.height * p_b.second.height);
            });
  peak_pairs_sorted.erase(peak_pairs_sorted.begin() + max_n,
                          peak_pairs_sorted.end());

  return peak_pairs_sorted;
}

// Calculate the gaussian contribution of the overlap between two points in
// one dimension.
double gaussian_contribution(double x_a,
                             double x_b,
                             double sigma_a,
                             double sigma_b) {
  double var_a = std::pow(sigma_a, 2);
  double var_b = std::pow(sigma_b, 2);

  double a = (var_a + var_b) / (var_a * var_b) *
             std::pow((x_a * var_b + x_b * var_a) / (var_a + var_b), 2);
  double b = (x_a * x_a) / var_a + (x_b * x_b) / var_b;

  return sigma_a * sigma_b * std::exp(0.5 * (a - b)) / std::sqrt(var_a + var_b);
}

/* mutating implementation that can be used to reduce the number of constructions and
 * destructions of warping surfaces in RANSAC */
void compute_warping_surf_impl(vector<double>& warping_surf,
                               const vector<peak_pair>& peaks_pairs,
                               const node& node_left,
                               const node& node_right) {
  double mz_begin = node_left.mz;
  double mz_end = node_right.mz;

  size_t n1 = node_left.mz_shifts.size();
  size_t n2 = node_right.mz_shifts.size();

  // Generate all warpings
  vector<pair<double, double>> warpings;
  warpings.reserve(n1 * n2);
  for (const auto& mz_i : node_left.mz_shifts) {
    for (const auto& mz_j : node_right.mz_shifts) {
      warpings.push_back({mz_begin + mz_i, mz_end + mz_j});
    }
  }

  // Compute warping surfs with all pairs
  for (const auto& pp : peaks_pairs) {
    const auto& p_ref = pp.first;
    const auto& p_warp = pp.second;

    size_t i = 0;
    for (const auto& warping : warpings) {
      // Warp peaks by interpolation
      double mz_warped =
          lerp(p_warp.mz, mz_begin, mz_end, warping.first, warping.second);

      // Compute overlap after warping
      double overlap = p_ref.height * p_warp.height *
                       gaussian_contribution(p_ref.mz, mz_warped,
                                             p_ref.sigma_mz, p_warp.sigma_mz);
      warping_surf[i] += overlap;
      ++i;
    }
  }
}

// pure interface to compute_warping_surf_impl
vector<double> compute_warping_surf(const vector<peak_pair>& peak_pairs,
                                    const node& node_left,
                                    const node& node_right) {
  vector<double> warping_surf(
      node_left.mz_shifts.size() * node_left.mz_shifts.size(), 0.0);
  compute_warping_surf_impl(warping_surf, peak_pairs, node_left, node_right);
  return warping_surf;
}

/* mutating implementation that can be used to reduce the number of constructions and
 * destructions of warping surf data in RANSAC */
vector<size_t> optimal_warping_impl(vector<vector<double>>& cum_surfs,
                                    size_t n_steps) {
  // Dynamic programming to find optimal combination of node moves
  size_t n_levels = cum_surfs.size();

  // Backward pass to compute cummulative levels
  for (size_t i = n_levels - 1; i > 0; --i) {
    const auto& surf_current = cum_surfs[i];
    auto& surf_next = cum_surfs[i - 1];

    // Find best right move for each left move;
    vector<double> row_max(n_steps, 0.0);
    auto it_current = surf_current.begin();
    for (size_t j = 0; j < n_steps; ++j) {
      const auto max_it = std::max_element(it_current, it_current + n_steps);
      row_max[j] = *max_it;
      it_current += n_steps;
    }

    /* Take maximum for each move of left node from current surf and add to
     * right moves of next surf */
    auto it_next = surf_next.begin();
    for (size_t j = 0; j < n_steps; ++j) {
      std::transform(it_next, it_next + n_steps, row_max.begin(), it_next,
                     std::plus<double>{});
      it_next += n_steps;
    }
  }

  // Then do a forward pass to find optimal path
  const auto& cum_first = cum_surfs.front();

  // Find optimal left and right node of first level
  size_t i_max = std::distance(
      cum_first.begin(), std::max_element(cum_first.begin(), cum_first.end()));
  size_t i_left = i_max / n_steps;
  size_t i_right = i_max % n_steps;

  vector<size_t> optimal_path;
  optimal_path.push_back(i_left);
  optimal_path.push_back(i_right);

  // Find optimal right move for each level given optimal left
  for (size_t i = 1; i < n_levels; ++i) {
    const auto& current_surf = cum_surfs[i];
    i_left = i_right;

    const auto it = current_surf.begin() + i_left * n_steps;
    i_right = std::distance(it, std::max_element(it, it + n_steps));

    optimal_path.push_back(i_right);
  }

  return optimal_path;
}

// pure interface to optimal_warping_impl
vector<size_t> optimal_warping(const vector<vector<double>>& warping_surfs,
                               size_t n_steps) {
  auto cum_surfs = warping_surfs;
  auto optimal_path = optimal_warping_impl(cum_surfs, n_steps);
  return optimal_path;
}

// random sample consensus (RANSAC) to detect outliers/inliers among peak
// pairs
ransac_result ransac(const vector<vector<peak_pair>>& peak_pairs,
                     const node_vec& warping_nodes,
                     size_t n_iterations,
                     size_t m,
                     double distance_threshold) {
  size_t n_steps = warping_nodes.front().mz_shifts.size();
  size_t n_segments = peak_pairs.size();
  size_t n_peaks =
      std::accumulate(peak_pairs.begin(), peak_pairs.end(), 0,
                      [](size_t sum, const auto& p) { return sum + p.size(); });

  // setup random sampling
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> udist(0, 100000);

  vector<double> errors(n_iterations, 0.0);
  vector<vector<bool>> inliers(n_iterations, vector<bool>(n_peaks, false));
  vector<vector<size_t>> optimal_paths;
  vector<vector<vector<size_t>>> maybe_inliers;

  // pre-allocate warping surfs
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
      for (size_t k = 0; k < m; ++k) {
        size_t x = udist(gen) % peak_pairs_segment.size();
        indices_j.push_back(x);
        maybe_inliers_j.push_back(peak_pairs_segment[x]);
      }

      // store indices of maybe_inliers
      indices_pairs.push_back(indices_j);

      compute_warping_surf_impl(warping_surf_j, maybe_inliers_j,
                                warping_nodes[j], warping_nodes[j + 1]);
    }

    auto optimal_path = optimal_warping_impl(warping_surfs, n_steps);

    // set of inliers and error for this iteration
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
        double d = std::pow((p.first.mz - mz_warped) / p.second.sigma_mz, 2);

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

vector<double> splat_peaks(const peak_vec& peaks,
                           const vector<double>& xi,
                           double splat_dist) {
  vector<double> y(xi.size(), 0.0);

  // minor optimization
  auto start = xi.begin();

  for (const auto& p_i : peaks) {
    const double d_mz = p_i.sigma_mz * splat_dist;
    const double mz_begin = p_i.mz - d_mz;
    const double mz_end = p_i.mz + d_mz;

    // peaks must be sorted by mz
    start = std::lower_bound(start, xi.end(), mz_begin);
    const auto stop = std::upper_bound(start, xi.end(), mz_end);

    for (auto it = start; it != stop; ++it) {
      const double d = (*it - p_i.mz) / p_i.sigma_mz;
      const double w = std::exp(-0.5 * std::pow(d, 2));

      const size_t i = std::distance(xi.begin(), it);
      y[i] += w * p_i.height;
    }
  }

  return y;
}

peak_vec merge_spectra(const peak_vec& s_a, const peak_vec& s_b) {
  peak_vec out;
  out.reserve(s_a.size() + s_b.size());

  auto comp_mz = [](const peak& p_a, const peak& p_b) {
    return p_a.mz < p_b.mz;
  };
  std::merge(s_a.begin(), s_a.end(), s_b.begin(), s_b.end(),
             std::back_inserter(out), comp_mz);

  return out;
}

vector<peak_pair> overlapping_peak_pairs(const peak_vec& s_a,
                                         const peak_vec& s_b,
                                         double max_dist) {
  if (s_a.empty() || s_b.empty()) {
    return {};
  }

  const auto s_m = merge_spectra(s_a, s_b);

  // s_a is reference
  auto ref_id = s_a.front().id;
  vector<peak_pair> pairs;

  for (size_t i = 0; i < s_m.size() - 1; ++i) {
    const auto& p_current = s_m[i];
    const auto& p_next = s_m[i + 1];

    if (p_current.id != p_next.id) {
      double max_current = p_current.mz + max_dist * p_current.sigma_mz;
      double min_next = p_next.mz - max_dist * p_next.sigma_mz;

      if (max_current > min_next) {  // the peaks overlap
        if (p_current.id == ref_id) {
          pairs.push_back({p_current, p_next});  // emplace back?
        } else {
          pairs.push_back({p_next, p_current});
        }
      }
    }
  }

  // TODO: sort before return?
  return pairs;
}

/* ----- legacy functions ----- */
double peak_overlap(const peak& peak_a, const peak& peak_b) {
  // Early return if the peaks do not intersect in the +/-3 * sigma_mz
  {
    double min_mz_a = peak_a.mz - 3 * peak_a.sigma_mz;
    double max_mz_a = peak_a.mz + 3 * peak_a.sigma_mz;
    double min_mz_b = peak_b.mz - 3 * peak_b.sigma_mz;
    double max_mz_b = peak_b.mz + 3 * peak_b.sigma_mz;

    if (max_mz_a < min_mz_b || max_mz_b < min_mz_a) {
      return 0;
    }
  }

  double mz_contrib = gaussian_contribution(peak_a.mz, peak_b.mz,
                                            peak_a.sigma_mz, peak_b.sigma_mz);

  return mz_contrib * peak_a.height * peak_b.height;
}

double compute_overlap(const peak_vec& peaks_a, const peak_vec& peaks_b) {
  double overlap = 0.0;
  for (const auto& p_a : peaks_a) {
    for (const auto& p_b : peaks_b) {
      overlap += peak_overlap(p_a, p_b);
    }
  }
  return overlap;
}

vector<double> compute_warping_surf(const peak_vec& peaks_src,
                                    const peak_vec& peaks_ref,
                                    const node& node_left,
                                    const node& node_right) {
  // All possible moves (TODO: add constraints)
  double mz_begin = node_left.mz;
  double mz_end = node_right.mz;

  size_t n1 = node_left.mz_shifts.size();
  size_t n2 = node_right.mz_shifts.size();

  vector<double> warping_surf;
  warping_surf.reserve(n1 * n2);

  for (const auto& mz_i : node_left.mz_shifts) {
    for (const auto& mz_j : node_right.mz_shifts) {
      auto peaks_warped = interpolate_peaks(peaks_src, mz_begin, mz_end,
                                            mz_begin + mz_i, mz_end + mz_j);
      warping_surf.push_back(compute_overlap(peaks_warped, peaks_ref));
    }
  }

  return warping_surf;
}

}  // namespace warp