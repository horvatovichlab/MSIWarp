#ifndef WARP_HPP
#define WARP_HPP

#include <cstdint>
#include <vector>

namespace warp {

/* */
enum class instrument { TOF, Orbitrap, FT_ICR, Quadrupole };

/* piecewise linear recalibration function composed of the m/z and m/z shift for
 * each warping node */
using recalibration_function = std::vector<std::pair<double, double>>;


/* */
struct peak {
  uint32_t id;
  double mz;
  double height;
  double sigma_mz;
};

/* */
struct node {
  double mz;
  double sigma_mz;
  std::vector<double> mz_shifts;
};

// TODO: add struct for warp parameters ?

using peak_pair = std::pair<peak, peak>;
using peak_vec = std::vector<peak>; // often a spectrum, but not always
using node_vec = std::vector<node>;

/* Finds the optimal combination node moves that align the source spectrum to
 * the reference spectrum. */
std::vector<size_t> find_optimal_warping(const peak_vec& s_r,
                                         const peak_vec& s_s,
                                         const node_vec& nodes,
                                         double epsilon);

/* Aligns a list of spectra to reference spectrum, s_r, without RANSAC and with
 constant nodes. */
std::vector<std::vector<size_t>> find_optimal_warpings(
    const std::vector<peak_vec>& spectra,
    const peak_vec& s_r,
    const node_vec& nodes,
    double epsilon);

/* Warps peak masses with piecewise linear interpolation between pairs of shifted
   warping nodes. */
peak_vec warp_peaks(const peak_vec& peaks,
                    const node_vec& nodes,
                    const std::vector<size_t>& moves);

/* Warps peak masses with piecewise linear interpolation between pairs of shifted
   warping nodes. */
peak_vec warp_peaks_unique(const peak_vec& peaks,
                           const recalibration_function& recal_func);

/* the search space of each warping node is its original mz +- its sigma times
 * n_steps  */
node_vec init_nodes(const std::vector<double>& mzs,
                    const std::vector<double>& sigmas,
                    uint32_t n_steps);

/* */
peak_vec peaks_between(const peak_vec& peaks, double mz_begin, double mz_end);

/* Pairs of matched peaks between bounds. Second peak in pair is assumed to be
   warped. */
std::vector<peak_pair> peak_pairs_between(const std::vector<peak_pair>& peaks,
                                          double mz_begin,
                                          double mz_end);

/* returns the peaks with top n highest intensity */
peak_vec peaks_top_n(const peak_vec& peaks, size_t n);

/* returns the peak pairs with top n highest intensity product */
std::vector<peak_pair> peak_pairs_top_n(
    const std::vector<peak_pair>& peak_pairs,
    size_t n);

/* computes the sum overlap between two spectra */
double compute_overlap(const peak_vec& peaks_a, const peak_vec& peaks_b);

/* computes (sum overlap) */
std::vector<double> compute_warping_surf(const peak_vec& peaks_src,
                                         const peak_vec& peaks_ref,
                                         const node& node_left,
                                         const node& node_right);

/* Alternative, more efficient warping */
std::vector<double> compute_warping_surf(
    const std::vector<std::pair<peak, peak>>& peaks_pairs,
    const node& node_left,
    const node& node_right);

/* finds optimal combination of node shifts with dynamic programming */
std::vector<size_t> optimal_warping(
    const std::vector<std::vector<double>>& warping_surfs,
    size_t n_steps);

/* Splats peaks on m/z sampling points xi */
std::vector<double> splat_peaks(const peak_vec& peaks,
                                const std::vector<double>& xi,
                                double splat_dist);

/* Merges two peak lists. Peaks lists must be sorted by m/z. */
peak_vec merge_spectra(const peak_vec&, const peak_vec&);

/* Finds overlapping peaks between spectrum s_a and s_b. Linear complexity. */
std::vector<std::pair<peak, peak>> overlapping_peak_pairs(const peak_vec& s_a,
                                                          const peak_vec& s_b,
                                                          double max_dist);

/* Linear interpolation. */
inline double lerp(double x,
                   double x_min,
                   double x_max,
                   double y_min,
                   double y_max) {
  double a = (x - x_min) / (x_max - x_min);
  return y_min + a * (y_max - y_min);
}

namespace detail {

/* Mutating implementation of compute_warping_surfs that can be used to reduce
 * the number of constructions and destructions of warping surfaces in RANSAC.
 */
void compute_warping_surf_impl(std::vector<double>& warping_surf,
                               const std::vector<peak_pair>& peaks_pairs,
                               const node& node_left,
                               const node& node_right);

/* Mutating implementationof optimal_warping that can be used to reduce the
 * number of constructions and destructions of warping surf data in RANSAC. */
std::vector<size_t> optimal_warping_impl(
    std::vector<std::vector<double>>& cum_surfs,
    size_t n_steps);

}  // namespace detail

}  // namespace warp

#endif