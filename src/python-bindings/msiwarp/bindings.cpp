#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <fstream>
#include <vector>

#include "util/msi_triplet.hpp"
#include "warp/warp.hpp"

namespace python_api {

/* ------ bindings to read/write functions of binary data ------ */

int write_triplets(const std::string& fname,
                   const std::vector<triplet::triplet>& ts) {
  std::ofstream ofs{fname, std::ios::binary};
  if (triplet::write_triplets(ofs, ts)) {
    return 0;
  }
  return -1;
}

int sort_write_triplets(const std::string& fname,
                        std::vector<triplet::triplet> ts) {
  std::sort(ts.begin(), ts.end(),
            [](const auto& t_a, const auto& t_b) { return t_a.mz < t_b.mz; });

  std::ofstream ofs{fname, std::ios::binary};
  if (triplet::write_triplets(ofs, ts)) {
    return 0;
  }
  return -1;
}

std::vector<triplet::triplet> read_triplets(const std::string& fname) {
  std::vector<triplet::triplet> ts;
  std::ifstream ifs{fname, std::ios::binary};
  if (triplet::read_triplets(ifs, &ts)) {
    return ts;
  }
  return {};
}

std::vector<triplet::triplet> get_triplets_range(const std::string& fname,
                                                 double begin,
                                                 double end) {
  std::ifstream ifs{fname, std::ios::binary};
  return triplet::get_range(ifs, begin, end);
}

/* ------ bindings to warp functions and data structures ------ */
std::vector<std::pair<warp::ransac_result, std::vector<warp::node>>>
align_ransac(const std::vector<std::vector<warp::peak>>& spectra,
             const std::vector<warp::peak> s_r,
             const warp::ransac_params& params) {
  return warp::align_ransac(spectra, s_r, params);
}

double compute_overlap(const std::vector<warp::peak>& reference_peaks,
                       const std::vector<warp::peak>& source_peaks) {
  return warp::compute_overlap(source_peaks, reference_peaks);
}

std::vector<warp::node> init_nodes(const std::vector<double>& mzs,
                                   const std::vector<double>& sigmas,
                                   uint32_t n_steps) {
  return warp::init_nodes(mzs, sigmas, n_steps);
}

std::vector<warp::peak> peaks_between(const std::vector<warp::peak>& peaks,
                                      double mz_begin,
                                      double mz_end) {
  return warp::peaks_between(peaks, mz_begin, mz_end);
}

std::vector<warp::peak_pair> peak_pairs_between(
    const std::vector<warp::peak_pair>& pairs,
    double mz_begin,
    double mz_end) {
  return warp::peak_pairs_between(pairs, mz_begin, mz_end);
}

std::vector<warp::peak> peaks_top_n(const std::vector<warp::peak>& peaks,
                                    size_t n) {
  return warp::peaks_top_n(peaks, n);
}

std::vector<warp::peak_pair> peak_pairs_top_n(
    const std::vector<warp::peak_pair>& peak_pairs,
    size_t n) {
  return warp::peak_pairs_top_n(peak_pairs, n);
}

// Naive implementation
std::vector<double> compute_warping_surf(
    const std::vector<warp::peak>& peaks_src,
    const std::vector<warp::peak>& peaks_ref,
    const warp::node& node_left,
    const warp::node& node_right) {
  return warp::compute_warping_surf(peaks_src, peaks_ref, node_left,
                                    node_right);
}

// Alternative and more efficient warping
std::vector<double> compute_warping_surf_pairs(
    const std::vector<std::pair<warp::peak, warp::peak>>& peaks_pairs,
    const warp::node& node_left,
    const warp::node& node_right) {
  return warp::compute_warping_surf(peaks_pairs, node_left, node_right);
}

// Find optimal combination of warpings with dynamic programming
std::vector<size_t> optimal_warping(
    const std::vector<std::vector<double>>& warping_surfs,
    size_t n_steps) {
  return warp::optimal_warping(warping_surfs, n_steps);
}

//
warp::ransac_result ransac(
    const std::vector<std::vector<warp::peak_pair>>& peak_pairs,
    const std::vector<warp::node>& warping_nodes,
    size_t n_iterations,
    size_t m,
    double distance_threshold) {
  if ((peak_pairs.size() + 1) != warping_nodes.size()) {
    return {{}, {}};
  }

  return warp::ransac(peak_pairs, warping_nodes, n_iterations, m,
                      distance_threshold);
}

// Splats peaks on m/z sampling points xi
std::vector<double> splat_peaks(const std::vector<warp::peak>& peaks,
                                const std::vector<double>& xi,
                                double splat_dist) {
  return warp::splat_peaks(peaks, xi, splat_dist);
}

std::vector<warp::peak> merge_spectra(const std::vector<warp::peak>& s_a,
                                      const std::vector<warp::peak>& s_b) {
  return warp::merge_spectra(s_a, s_b);
}

/* Finds overlapping peaks between spectrum s_a and s_b. Linear complexity. */
std::vector<std::pair<warp::peak, warp::peak>> overlapping_peak_pairs(
    const std::vector<warp::peak>& s_a,
    const std::vector<warp::peak>& s_b,
    double max_dist) {
  return warp::overlapping_peak_pairs(s_a, s_b, max_dist);
}

}  // namespace python_api

namespace py = pybind11;

PYBIND11_MODULE(msiwarp, m) {
  m.doc() = R"pbdoc(
            Pybind11 example plugin
            -----------------------
            .. currentmodule:: cmake_example
            .. autosummary::
               :toctree: _generate
               add
               subtract
        )pbdoc";

  /* ---- python bindings to utility functions ---- */
  m.def("write_triplets", &python_api::write_triplets, R"pbdoc(
            write triplets (index, mz, height) to binary file
            Some other explanation about the subtract function.
        )pbdoc");

  m.def("sort_write_triplets", &python_api::sort_write_triplets, R"pbdoc(
            sort triplets by m/z and then write them to binary file
            Some other explanation about the subtract function.
        )pbdoc");

  m.def("read_triplets", &python_api::read_triplets, R"pbdoc(
            read triplets (index, mz, height) from binary file
        )pbdoc");

  m.def("get_triplets_range", &python_api::get_triplets_range, R"pbdoc(
            get all triplets within range (triplets must be sorted by m/z)
        )pbdoc");

  py::class_<triplet::triplet>(m, "msi_triplet")
      .def(py::init<uint32_t, double, float>())
      .def_readonly("index", &triplet::triplet::index)
      .def_readonly("mz", &triplet::triplet::mz)
      .def_readonly("height", &triplet::triplet::height)
      .def("__repr__", [](const triplet::triplet& t) {
        return "triplet <index: " + std::to_string(t.index) +
               ", mz: " + std::to_string(t.mz) +
               ", height: " + std::to_string(t.height) + ">";
      });

  /* ---- python bindings to warp functions ---- */
  m.def("find_optimal_spectrum_warping", &warp::find_optimal_warping,
        R"pbdoc(
        Find optimal warpings from sample spectrum to reference spectrum.
    )pbdoc");

  m.def("find_optimal_spectra_warpings", &warp::find_optimal_warpings, R"pbdoc(
        Find optimal warpings from sample spectra to reference spectrum.
        Warpings nodes are constant for all spectra.
    )pbdoc");

  m.def("warp_peaks", &warp::warp_peaks, R"pbdoc(
        Warp peaks with fixed warping nodes and optimal moves.
    )pbdoc");

  m.def("initialize_nodes", &python_api::init_nodes, R"pbdoc(
        Initializes warping nodes. Warping will be performed
         between each pair of nodes.
    )pbdoc");

  m.def("peaks_between", &python_api::peaks_between, R"pbdoc(
        Finds all peaks between mz_begin and mz_end in peaks.
        Peaks are required to be sorted by mz.
    )pbdoc");

  m.def("peak_pairs_between", &python_api::peak_pairs_between, R"pbdoc(
        Finds all pairs in range [mz_begin, mz_end).
        Pairs are required to be sorted by mz.
    )pbdoc");

  m.def("peaks_top_n", &python_api::peaks_top_n, R"pbdoc(
        Returns the top n peaks.
    )pbdoc");

  m.def("peak_pairs_top_n", &python_api::peak_pairs_top_n, R"pbdoc(
        Returns the top n peaks.
    )pbdoc");

  m.def("compute_overlap", &python_api::compute_overlap, R"pbdoc(
        Computes the total overlap between peaks in a and b.
    )pbdoc");

  m.def("compute_warping_surf", &python_api::compute_warping_surf, R"pbdoc(
        Computes the overlap surface for all possible node moves.
    )pbdoc");

  m.def("compute_warping_surf_pairs", &python_api::compute_warping_surf_pairs,
        R"pbdoc(
        Computes the overlap surface for all possible node moves.
    )pbdoc");

  m.def("optimal_warping", &python_api::optimal_warping, R"pbdoc(
        Find the optimal combination of warpings.
    )pbdoc");

  m.def("ransac", &python_api::ransac, R"pbdoc(
        Random sampling consensus of preliminary peak matches.
    )pbdoc");

  m.def("splat_peaks", &python_api::splat_peaks, R"pbdoc(
        Splats peaks on sampling points, xi. Peaks must be sorted by m/z.
    )pbdoc");

  m.def("merge_spectra", &python_api::merge_spectra, R"pbdoc(
        Merges peak list p_a and p_b. p_a and p_b must be sorted by m/z.
    )pbdoc");

  m.def("overlapping_peak_pairs", &python_api::overlapping_peak_pairs, R"pbdoc(
        Get all overlapping peak pairs between s_a and s_b. Peaks must be sorted by m/z.
    )pbdoc");

  py::class_<warp::peak>(m, "peak")
      .def_readonly("id", &warp::peak::id)
      .def_readonly("mz", &warp::peak::mz)
      .def_readonly("height", &warp::peak::height)
      .def_readonly("sigma_mz", &warp::peak::sigma_mz)
      .def(py::init<uint64_t, double, double, double>())
      .def("__repr__", [](const warp::peak& p) {
        return "Peak <entity_id: " + std::to_string(p.id) +
               ", mz: " + std::to_string(p.mz) +
               ", height: " + std::to_string(p.height) +
               ", sigma_mz: " + std::to_string(p.sigma_mz) + ">";
      });

  py::class_<warp::node>(m, "node")
      .def_readonly("mz", &warp::node::mz)
      .def_readonly("sigma_mz", &warp::node::sigma_mz)
      .def_readonly("mz_shifts", &warp::node::mz_shifts)
      .def("__repr__", [](const warp::node& n) {
        return "Node <mz: " + std::to_string(n.mz) +
               ", sigma_mz: " + std::to_string(n.sigma_mz) + ", mz_shifts: [" +
               std::to_string(n.mz_shifts.front()) + ", " +
               std::to_string(n.mz_shifts.back()) + "]>";
      });

  py::class_<warp::ransac_result>(m, "ransac_result")
      .def_readonly("errors", &warp::ransac_result::errors)
      .def_readonly("inliers", &warp::ransac_result::inliers)
      .def_readonly("models", &warp::ransac_result::models)
      .def_readonly("maybe_inliers", &warp::ransac_result::maybe_inliers)
      .def("__repr__", [](const warp::ransac_result& rr) {
        return "RANSAC results <number of iterations: " +
               std::to_string(rr.errors.size()) + ">";
      });

  py::class_<warp::ransac_params>(m, "ransac_params")
      .def_readonly("n_segments", &warp::ransac_params::n_segments)
      .def_readonly("n_iterations", &warp::ransac_params::n_iterations)
      .def_readonly("n_samples", &warp::ransac_params::n_samples)
      .def_readonly("min_matched_peaks",
                    &warp::ransac_params::min_matched_peaks)
      .def_readonly("distance_threshold",
                    &warp::ransac_params::distance_threshold)
      .def_readonly("mz_begin", &warp::ransac_params::mz_begin)
      .def_readonly("mz_end", &warp::ransac_params::mz_end)
      .def("__repr__", [](const warp::ransac_params& p) {
        return "RANSAC parameters: <n_segments: " +
               std::to_string(p.n_segments) +
               ", n_iterations: " + std::to_string(p.n_iterations) +
               ", n_samples: " + std::to_string(p.n_samples) +
               ", min_matched_peaks: " + std::to_string(p.min_matched_peaks) +
               ", distance_threshold: " + std::to_string(p.distance_threshold) +
               ", mz_begin: " + std::to_string(p.mz_begin) +
               ", mz_end: " + std::to_string(p.mz_end) + ">";
      });

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}