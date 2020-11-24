#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <fstream>
#include <vector>
#include <numeric>

#include "util/msi_triplet.hpp"
#include "warp/warp.hpp"
#include "warp/warp_util.hpp"
#include "warp/ransac.hpp"

namespace python_api {

/* ------ bindings to read/write of binary data ------ */

bool write_triplets(const std::string& fname,
                    const std::vector<triplet::triplet>& ts) {
  std::ofstream ofs{fname, std::ios::binary};
  return triplet::write_triplets(ofs, ts);
}

bool sort_write_triplets(const std::string& fname,
                         std::vector<triplet::triplet>& ts) {
  std::sort(ts.begin(), ts.end(),
            [](const auto& t_a, const auto& t_b) { return t_a.mz < t_b.mz; });

  return write_triplets(fname, ts);
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
                                                 double begin, double end) {
  std::ifstream ifs{fname, std::ios::binary};
  return triplet::get_range(ifs, begin, end);
}

bool spectra_to_triplets(const std::string& fname,
                         const std::vector<warp::peak_vec>& spectra) {
  size_t n_triplets =
      std::accumulate(spectra.begin(), spectra.end(), size_t{0},
                      [](size_t n, const auto& s) { return n + s.size(); });

  std::vector<triplet::triplet> ts;
  ts.reserve(n_triplets);

  for (const auto& s : spectra) {
    for (const auto& p : s) {
      ts.emplace_back(
          triplet::triplet{p.id, p.mz, static_cast<float>(p.height)});
    }
  }

  return sort_write_triplets(fname, ts);
}

/* ------ bindings to warp functions and data structures ------ */

// Compute warping surface for peak lists and warping nodes
std::vector<double> compute_warping_surf(
    const std::vector<warp::peak>& peaks_src,
    const std::vector<warp::peak>& peaks_ref, const warp::node& node_left,
    const warp::node& node_right) {
  return warp::compute_warping_surf(peaks_src, peaks_ref, node_left,
                                    node_right);
}

// Compute warping surface for a list of matched peaks and warping nodes
std::vector<double> compute_warping_surf_pairs(
    const std::vector<std::pair<warp::peak, warp::peak>>& peaks_pairs,
    const warp::node& node_left, const warp::node& node_right) {
  return warp::compute_warping_surf(peaks_pairs, node_left, node_right);
}

// // TODO: unfinished interface to RANSAC
// warp::ransac_result ransac(
//     const std::vector<std::vector<warp::peak_pair>>& peak_pairs,
//     const std::vector<warp::node>& warping_nodes, size_t n_iterations, size_t m,
//     double distance_threshold) {
//   if ((peak_pairs.size() + 1) != warping_nodes.size()) {
//     return {{}, {}};
//   }
//   return {{}, {}};

//   return warp::ransac(peak_pairs, warping_nodes, n_iterations, m,
//                       distance_threshold);
// }

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
        )pbdoc");

  m.def("sort_write_triplets", &python_api::sort_write_triplets, R"pbdoc(
            sort triplets by m/z and then write them to binary file
        )pbdoc");

  m.def("read_triplets", &python_api::read_triplets, R"pbdoc(
            read triplets (index, mz, height) from binary file
        )pbdoc");

  m.def("get_triplets_range", &python_api::get_triplets_range, R"pbdoc(
            get all triplets within range (triplets must be sorted by m/z)
        )pbdoc");

  m.def("spectra_to_triplets", &python_api::spectra_to_triplets, R"pbdoc(
            sort and write all peaks from list of spectra to file
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

  m.def("warp_peaks_unique", &warp::warp_peaks_unique, R"pbdoc(
        Warp peaks with fixed warping nodes and optimal moves.
    )pbdoc");

  m.def("initialize_nodes", &warp::init_nodes, R"pbdoc(
        Initializes warping nodes. Warping will be performed
         between each pair of nodes.
    )pbdoc");

  m.def("peaks_between", &warp::peaks_between, R"pbdoc(
        Finds all peaks between mz_begin and mz_end in peaks.
        Peaks are required to be sorted by mz.
    )pbdoc");

  m.def("peak_pairs_between", &warp::peak_pairs_between, R"pbdoc(
        Finds all pairs in range [mz_begin, mz_end).
        Pairs are required to be sorted by mz.
    )pbdoc");

  m.def("peaks_top_n", &warp::peaks_top_n, R"pbdoc(
        Returns the top n peaks.
    )pbdoc");

  m.def("peak_pairs_top_n", &warp::peak_pairs_top_n, R"pbdoc(
        Returns the top n peaks.
    )pbdoc");

  m.def("compute_overlap", &warp::compute_overlap, R"pbdoc(
        Computes the total overlap between peaks in a and b.
    )pbdoc");

  m.def("compute_warping_surf", &python_api::compute_warping_surf, R"pbdoc(
        Computes the overlap surface for all possible node moves.
    )pbdoc");

  m.def("compute_warping_surf_pairs", &python_api::compute_warping_surf_pairs,
        R"pbdoc(
        Computes the overlap surface for all possible node moves.
    )pbdoc");

  m.def("optimal_warping", &warp::optimal_warping, R"pbdoc(
        Find the optimal combination of warpings.
    )pbdoc");

  // m.def("ransac", &python_api::ransac, R"pbdoc(
  //       Random sampling consensus of preliminary peak matches.
  //   )pbdoc");

  m.def("splat_peaks", &warp::splat_peaks, R"pbdoc(
        Splats peaks on sampling points, xi. Peaks must be sorted by m/z.
    )pbdoc");

  m.def("merge_spectra", &warp::merge_spectra, R"pbdoc(
        Merges peak list p_a and p_b. p_a and p_b must be sorted by m/z.
    )pbdoc");

  m.def("overlapping_peak_pairs", &warp::overlapping_peak_pairs, R"pbdoc(
        Get all overlapping peak pairs between s_a and s_b. Peaks must be sorted by m/z.
    )pbdoc");

  py::enum_<warp::instrument>(m, "Instrument")
      .value("TOF", warp::instrument::TOF)
      .value("Orbitrap", warp::instrument::Orbitrap)
      .value("FT-ICR", warp::instrument::FT_ICR)
      .value("Quadrupole", warp::instrument::Quadrupole)
      .export_values();

  py::class_<warp::peak>(m, "peak")
      .def_readonly("id", &warp::peak::id)
      .def_readonly("mz", &warp::peak::mz)
      .def_readonly("height", &warp::peak::height)
      .def_readonly("sigma_mz", &warp::peak::sigma)
      .def(py::init<uint32_t, double, double, double>())
      .def("__repr__", [](const warp::peak& p) {
        return "Peak <entity_id: " + std::to_string(p.id) +
               ", mz: " + std::to_string(p.mz) +
               ", height: " + std::to_string(p.height) +
               ", sigma_mz: " + std::to_string(p.sigma) + ">";
      });

  py::class_<warp::node>(m, "node")
      .def_readonly("mz", &warp::node::mz)
      .def_readonly("slack", &warp::node::slack)
      .def_readonly("mz_shifts", &warp::node::mz_shifts)
      .def("__repr__", [](const warp::node& n) {
        return "Node <mz: " + std::to_string(n.mz) +
               ", slack: " + std::to_string(n.slack) + ", mz_shifts: [" +
               std::to_string(n.mz_shifts.front()) + ", " +
               std::to_string(n.mz_shifts.back()) + "]>";
      });

  /* ---- python bindings to warp::util functions ---- */
  py::class_<warp::util::params_uniform>(m, "params_uniform")
      .def_readonly("instrument", &warp::util::params_uniform::inst)
      .def_readonly("n_steps", &warp::util::params_uniform::n_steps)
      .def_readonly("n_peaks", &warp::util::params_uniform::n_peaks)
      .def_readonly("max_nodes", &warp::util::params_uniform::max_nodes)
      .def_readonly("mz_begin", &warp::util::params_uniform::mz_begin)
      .def_readonly("mz_end", &warp::util::params_uniform::mz_end)
      .def_readonly("slack", &warp::util::params_uniform::slack)
      .def(py::init<warp::instrument, size_t, size_t, size_t, double, double,
                    double>())
      .def("__repr__", [](const warp::util::params_uniform& p) {
        return "parameters for node placement function";
      });

  m.def("get_warping_nodes_uniform", &warp::util::get_warping_nodes_uniform,
        R"pbdoc(
        Place warping nodes so that all segments have the same number of peak mathces.
    )pbdoc");

  m.def("find_optimal_warpings_uni", &warp::util::find_optimal_warpings_uni,
        R"pbdoc(
        TODO: add documentation.
    )pbdoc");

  /* ---- python bindings to warp::ransac functions ---- */
  py::class_<warp::ransac::ransac_result>(m, "ransac_result")
      .def_readonly("errors", &warp::ransac::ransac_result::errors)
      .def_readonly("inliers", &warp::ransac::ransac_result::inliers)
      .def_readonly("models", &warp::ransac::ransac_result::models)
      .def_readonly("maybe_inliers", &warp::ransac::ransac_result::maybe_inliers)
      .def("__repr__", [](const warp::ransac::ransac_result& rr) {
        return "RANSAC results <number of iterations: " +
               std::to_string(rr.errors.size()) + ">";
      });

  py::class_<warp::ransac::params>(m, "ransac_params")
      .def_readonly("n_segments", &warp::ransac::params::n_segments)
      .def_readonly("n_iterations", &warp::ransac::params::n_iterations)
      .def_readonly("n_samples", &warp::ransac::params::n_samples)
      .def_readonly("min_matched_peaks",
                    &warp::ransac::params::min_matched_peaks)
      .def_readonly("n_steps", &warp::ransac::params::n_steps)
      .def_readonly("slack", &warp::ransac::params::slack)
      .def_readonly("mz_begin", &warp::ransac::params::mz_begin)
      .def_readonly("mz_end", &warp::ransac::params::mz_end)
      .def_readonly("distance_threshold",
                    &warp::ransac::params::distance_threshold)
      .def(py::init<size_t, size_t, size_t, size_t, size_t, double, double,
                    double, double>())
      .def("__repr__", [](const warp::ransac::params& p) {
        return "RANSAC parameters <n_segments:" + std::to_string(p.n_segments) +
               "n_iterations:" + std::to_string(p.n_iterations) +
               "n_samples:" + std::to_string(p.n_samples) +
               "min_peak_matches:" + std::to_string(p.min_matched_peaks) +
               "n_steps:" + std::to_string(p.n_steps) +
               "slack:" + std::to_string(p.slack) +
               "mz_begin:" + std::to_string(p.mz_begin) +
               "mz_end:" + std::to_string(p.mz_end) +
               "distance_threshold:" + std::to_string(p.distance_threshold) +
               ">";
      });

  m.def("ransac_pairs", &warp::ransac::ransac_pairs,
        R"pbdoc(
        Remove spurious peak matches with RANSAC.
    )pbdoc");

  m.def("ransac", &warp::ransac::ransac,
        R"pbdoc(
        Perform RANSAC on pairs and return results for each iteration.
    )pbdoc");

  m.def("align_ransac_uniform", &warp::ransac::align_ransac_uniform,
        R"pbdoc(
        Perform RANSAC on pairs and return results for each iteration.
    )pbdoc");

  m.def("find_optimal_warpings_ransac_uni",
        &warp::ransac::find_optimal_warpings_uni,
        R"pbdoc(
        Perform RANSAC on pairs and return results for each iteration.
    )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}