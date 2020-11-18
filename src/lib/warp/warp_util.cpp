#include "warp_util.hpp"

#include <cmath>
#include <algorithm>

#include "util/parallel.hpp"
     
namespace warp::util {

// TODO: move to warp.cpp?
double get_mz_scaling(double mz, instrument inst) {  
  switch (inst) {
    case instrument::TOF        : return mz;      
    case instrument::Orbitrap   : return std::pow(mz, 1.5);    
    case instrument::FT_ICR     : return std::pow(mz, 2.0);
    case instrument::Quadrupole : return 1.0;
  }
}

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
    size_t stride = pairs.size() / n;

    mzs.resize(n);
    mzs.front() = p.mz_begin;
    mzs.back() = p.mz_end;

    size_t j = stride;
    for (size_t i = 1; i < n - 1; ++i) {
      mzs[i] = pairs[j].first.mz;
      j += stride;
    }
  }

  // slack = ds * n_steps; each node is shifted up and down in m/z
  double ds = p.slack / p.n_steps;

  std::vector<double> deltas(mzs.size());
  std::transform(mzs.begin(), mzs.end(), deltas.begin(), [ds, &p](double mz) {
    return ds * get_mz_scaling(mz, p.inst);
  });

  return warp::init_nodes(mzs, deltas, p.n_steps);
}

// TODO: duplicate of warp.cpp
std::vector<size_t> find_optimal_warping_pairs(const std::vector<peak_pair>& ps,
                                               const node_vec& nodes) {
  std::vector<std::vector<double>> warping_surfs;
  for (size_t i = 0; i < nodes.size() - 1; ++i) {
    const auto& n_l = nodes[i];
    const auto& n_r = nodes[i + 1];

    const auto ps_i = peak_pairs_between(ps, n_l.mz, n_r.mz);
    const auto w_i = compute_warping_surf(ps_i, n_l, n_r);
    warping_surfs.push_back(w_i);
  }

  size_t n_steps = nodes.front().mz_shifts.size();
  return optimal_warping(warping_surfs, n_steps);
}

/* Function template for aligning with a custom node placement function Func
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