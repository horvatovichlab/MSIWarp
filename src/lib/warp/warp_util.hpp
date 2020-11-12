#ifndef WARP_UTIL_HPP
#define WARP_UTIL_HPP

#include "warp.hpp"

namespace warp::util {

using recalibration_function = std::vector<std::pair<double, double>>;

/* parameters for uniquely placed nodes */
struct node_params {
  instrument inst;
  size_t n_steps;
  size_t n_peaks;
  size_t max_nodes;
  double mz_begin;
  double mz_end;
  double slack;
};

/* TODO: move this to warp.hpp / warp.cpp */   
double get_mz_scaling(double mz, instrument inst);

/* */
node_vec get_warping_nodes_uniform(const std::vector<peak_pair>& ps,
                                   const node_params& params);

/* */
node_vec get_warping_nodes_density(const std::vector<peak_pair>& ps,
                                   double sigma_1,
                                   double epsilon,
                                   size_t n_steps,
                                   instrument inst);

/* */
std::vector<recalibration_function> find_optimal_warpings_uni(
    const std::vector<peak_vec>& spectra,
    const peak_vec& s_ref,
    const node_params& params,
    double epsilon,
    size_t n_cores);

} // namespace warp::util

#endif // WARP_UTIL_HPP