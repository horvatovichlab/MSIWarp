#ifndef RANSAC_HPP
#define RANSAC_HPP

#include <vector>

#include "warp.hpp"
#include "warp_util.hpp"

namespace warp::ransac {

/* TODO: document
example:
*/
struct ransac_params {
  size_t n_segments;
  size_t n_iterations;
  size_t n_samples;
  size_t min_matched_peaks;
  double distance_threshold;
  double mz_begin;
  double mz_end;
};

/* TODO: document
example:
*/
struct ransac_result {
  std::vector<double> errors;
  std::vector<std::vector<bool>> inliers;
  std::vector<std::vector<size_t>> models;
  std::vector<std::vector<std::vector<size_t>>> maybe_inliers;
};

/* */
std::vector<peak_pair> ransac_pairs(const std::vector<peak_pair>& pairs,
                                    const ransac_params& params,
                                    const util::params_uniform& node_params);

/* */
ransac_result ransac(
    const std::vector<std::vector<warp::peak_pair>>& peak_pairs,
    const warp::node_vec& warping_nodes,
    size_t n_iterations,
    size_t m,
    double distance_threshold);

/* */
std::vector<std::pair<ransac_result, warp::node_vec>> align_ransac(
    const std::vector<warp::peak_vec>& spectra,
    const warp::peak_vec s_r,
    const ransac_params& params);



}  // namespace warp::ransac

#endif