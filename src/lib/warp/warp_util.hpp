#ifndef WARP_UTIL_HPP
#define WARP_UTIL_HPP

#include "warp.hpp"

namespace warp::util {

/* Parameters for uniform node placement. Nodes are placed so that all segments
 * contain the same number of peaks. The parameter 'n_peaks' is the lower bound
 * for the number of peaks per segment. If the number of peak matches is
 * greater than 'n_peaks' * 'max_nodes', the number of peaks is equal to the
 * number of peak matches divided by 'max_nodes'. Slack increases with m/z.
 */
struct params_uniform {
  instrument inst;
  size_t n_steps;
  size_t n_peaks;
  size_t max_nodes;
  double mz_begin;
  double mz_end;
  double slack;
};

/* Parameters for density-based node placement. */
struct params_density {
  instrument inst;
  double bandwidth;
  double threshold;
  double mz_begin;
  double mz_end;
  double slack;  
  size_t n_steps;
};

/* */
node_vec get_warping_nodes_uniform(const std::vector<peak_pair>& ps,
                                   const params_uniform& params);

/* */
node_vec get_warping_nodes_density(const std::vector<peak_pair>& ps,
                                   const params_density& params);

/* */
std::vector<recalibration_function> find_optimal_warpings_uni(
    const std::vector<peak_vec>& spectra,
    const peak_vec& s_ref,
    const params_uniform& params,
    double epsilon,
    size_t n_cores);

}  // namespace warp::util

#endif  // WARP_UTIL_HPP