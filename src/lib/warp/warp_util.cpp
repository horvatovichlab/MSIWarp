#include "warp_util.hpp"

#include <cmath>
#include <algorithm>

#include "util/parallel.hpp"
     
namespace warp::util {

node_vec get_warping_nodes_uniform(const std::vector<peak_pair>& pairs,
                                   const params_uniform& p) {
  size_t n_pairs = pairs.size();
  std::vector<double> mzs;

  if (n_pairs < p.n_peaks * 2 || p.max_nodes < 3) {
    // use only one segment
    mzs = {p.mz_begin, p.mz_end};
  } else {
    // use...
    size_t n = std::min(p.max_nodes, n_pairs / p.n_peaks);
    size_t stride = pairs.size() / (n - 1);

    mzs.resize(n);
    mzs.front() = p.mz_begin;
    mzs.back() = p.mz_end;

    size_t j = stride;
    for (size_t i = 1; i < n - 1; ++i) {
      mzs[i] = pairs[j].second.mz;
      j += stride;
    }
  }

  std::vector<double> slacks(mzs.size());
  std::transform(mzs.begin(), mzs.end(), slacks.begin(), [&p](double mz) {
    return p.slack * get_mz_scaling(mz, p.inst);
  });

  return warp::init_nodes(mzs, slacks, p.n_steps);
}

node_vec get_warping_nodes_density(const std::vector<peak_pair>& pairs,
                                   const params_density& p) {
  // TODO: implement
  return {};
}

/* Function template for aligning with a custom node placement function, Func,
   that takes two arguments: a list of peak pairs and the parameter struct
   expected by the function. */
template <class Func, class Params>
recalibration_function find_optimal_warping_unique(const peak_vec& s_r,
                                                   const peak_vec& s_s,
                                                   double epsilon,
                                                   Func node_func,
                                                   const Params& params) {
  const auto pairs = overlapping_peak_pairs(s_r, s_s, epsilon);

  // place nodes uniquely for each spectrum
  const auto nodes = node_func(pairs, params);

  //
  const auto optimal_warping = find_optimal_warping_pairs(pairs, nodes);

  // With uniquely placed nodes we must return the nodes' m/z in addition to
  // their optimal shifts.
  recalibration_function recal_func;
  recal_func.reserve(nodes.size());
  for (size_t i = 0; i < nodes.size(); ++i) {
    recal_func.push_back({nodes[i].mz, nodes[i].mz_shifts[optimal_warping[i]]});
  }

  return recal_func;
}

/* Template specialization for uniform node placement */
std::vector<recalibration_function> find_optimal_warpings_uni(
    const std::vector<peak_vec>& spectra,
    const peak_vec& s_r,
    const params_uniform& params,
    double epsilon,
    size_t n_cores) {
  //
  auto out =
      par::for_each(spectra.size() / n_cores, spectra.begin(), spectra.end(),
                    [&](const auto& s_i) {
                      return find_optimal_warping_unique(
                          s_r, s_i, epsilon, get_warping_nodes_uniform, params);
                    });
  return out;
}

/* TODO: add density based node placement function */

}  // namespace warp::util