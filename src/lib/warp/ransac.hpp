#ifndef RANSAC_HPP
#define RANSAC_HPP

#include <vector>

#include "warp.hpp"
#include "warp_util.hpp"

namespace warp::ransac {

/* TODO: document
 * example: */
struct params {
  size_t n_segments;
  size_t n_iterations;
  size_t n_samples;
  size_t min_matched_peaks;
  size_t n_steps;
  double slack;
  double mz_begin;
  double mz_end;
  double distance_threshold;
};

/* TODO: document
 * example: */
struct ransac_result {
  std::vector<double> errors;
  std::vector<std::vector<bool>> inliers;
  std::vector<std::vector<size_t>> models;
  std::vector<std::vector<std::vector<size_t>>> maybe_inliers;
};

/* Remove spurious peak matches with RANSAC. */
std::vector<peak_pair> ransac_pairs(const std::vector<peak_pair>& pairs,
                                    const params& params,
                                    const util::params_uniform& node_params);

/* Align a list of spectra with RANSAC and unique node placement. */
std::vector<std::pair<ransac_result, warp::node_vec>> align_ransac(
    const std::vector<warp::peak_vec>& spectra,
    const warp::peak_vec s_r,
    const params& params);

/* RANSAC implementation. Peak pairs must be grouped by warping segment. */
ransac_result ransac(
    const std::vector<std::vector<warp::peak_pair>>& peak_pairs,
    const warp::node_vec& warping_nodes,
    size_t n_iterations,
    size_t n_samples,
    double distance_threshold);

/* */
recalibration_function align_ransac_uniform(
    const std::vector<peak_pair>& pairs, const params& ransac_params,
    const util::params_uniform& node_params);

/* Find optimal warpings of specta relative to the reference spectrum, s_ref,
 * with RANSAC and uniform node placement.  */
std::vector<recalibration_function> find_optimal_warpings_uni(
    const std::vector<peak_vec>& spectra,
    const peak_vec& s_ref,
    const params& ransac_params,
    const util::params_uniform& params,
    double epsilon,
    size_t n_cores);

}  // namespace warp::ransac

#endif