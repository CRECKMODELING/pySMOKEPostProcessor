/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Authors: Lorenzo Giardini <lorenzo.giardini@polimi.it>                |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
| License                                                                 |
|                                                                         |
|   Copyright(C) 2016-2012  Alberto Cuoci                                 |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
constexpr double Rgas = 8314.462618;   // universal gas constant, J/(kmol.K)
constexpr double NA = 6.02214076e23;  // #/mol
}  // namespace

Soot::Soot() { data_ = nullptr; }

void Soot::SetDatabase(ProfilesDatabase* data) { data_ = data; }

// local value to print the PSD from
unsigned int Soot::AbscissaIndex(double local_value) const {
  const std::vector<double>& abscissa = data_->additional[0];
  const unsigned int npts = static_cast<unsigned int>(data_->number_of_abscissas_);
  for (unsigned int j = 0; j < npts; j++)
    if (abscissa[j] >= local_value) return j;
  if (npts == 0) return 0;
  return npts - 1;
}

// Shared back half of both distributions: number density per selected bin at
// `point`, then merge near-equal diameters into fixed sections and form
// dN/dlog10(d). `size_m` is the size coordinate (metres) of each entry in `bins`.
std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, double>>
Soot::BuildDistribution(const std::vector<unsigned int>& bins,
                        const std::vector<double>& size_m, const unsigned int point,
                        const double merge_tol) const {
  const double T = data_->additional[data_->index_T][point];
  const double P = data_->additional[data_->index_P][point];
  const double MWmix = data_->additional[data_->index_MW][point];
  const double rho = P * MWmix / (Rgas * T);        // kg/m3, ideal gas
  const double number_factor = rho * 1000.0 * NA;  // * omega_i / MW_i -> #/m3

  const unsigned int n = static_cast<unsigned int>(bins.size());
  std::vector<double> d(n), N(n);
  for (unsigned int k = 0; k < n; k++) {
    const unsigned int sp = data_->soot_bin_index_[bins[k]];  // gas species index
    double Ni = 0.0;
    if (sp < data_->omega.size() && sp < data_->mw_species_.size() &&
        data_->mw_species_[sp] > 0.0)
      Ni = data_->omega[sp][point] * number_factor / data_->mw_species_[sp];
    d[k] = size_m[k];
    N[k] = Ni;
  }

  std::vector<unsigned int> order(n);
  for (unsigned int i = 0; i < n; i++) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](unsigned int a, unsigned int b) { return d[a] < d[b]; });

  // Greedy fixed sections: keep growing the current section while
  // d / section_min <= 1 + merge_tol; representative diameter is the geometric
  // mean, so it does not move with the current concentrations.
  std::vector<double> col_d, col_dmin, col_dmax, col_N, col_nbins;
  double section_min = 0.0, log_sum = 0.0, dmin = 0.0, dmax = 0.0, Nsum = 0.0;
  unsigned int count = 0;
  const auto flush = [&]() {
    if (count == 0) return;
    col_d.push_back(std::pow(10.0, log_sum / count) * 1e9);
    col_dmin.push_back(dmin * 1e9);
    col_dmax.push_back(dmax * 1e9);
    col_N.push_back(Nsum);
    col_nbins.push_back(static_cast<double>(count));
  };
  for (unsigned int oi = 0; oi < n; oi++) {
    const unsigned int i = order[oi];
    if (count > 0 && d[i] / section_min > 1.0 + merge_tol) {
      flush();
      count = 0;
      log_sum = 0.0;
      Nsum = 0.0;
    }
    if (count == 0) {
      section_min = d[i];
      dmin = d[i];
      dmax = d[i];
    }
    dmin = std::min(dmin, d[i]);
    dmax = std::max(dmax, d[i]);
    log_sum += std::log10(d[i]);
    Nsum += N[i];
    count++;
  }
  flush();

  // Centred log spacing on the sorted representative diameters (nm): one-sided at
  // the ends, centred inside.
  const size_t m = col_d.size();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> dlog(m, nan), dndlog(m, nan);
  if (m >= 2) {
    dlog[0] = std::log10(col_d[1]) - std::log10(col_d[0]);
    dlog[m - 1] = std::log10(col_d[m - 1]) - std::log10(col_d[m - 2]);
  }
  for (size_t i = 1; i + 1 < m; i++)
    dlog[i] = 0.5 * (std::log10(col_d[i + 1]) - std::log10(col_d[i - 1]));
  for (size_t i = 0; i < m; i++)
    if (std::fabs(dlog[i]) > 1e-12) dndlog[i] = col_N[i] / dlog[i];

  std::map<std::string, std::vector<double>> cols;
  cols["d_nm"] = col_d;
  cols["d_min_nm"] = col_dmin;
  cols["d_max_nm"] = col_dmax;
  cols["N_per_m3"] = col_N;
  cols["n_bins"] = col_nbins;
  cols["Dlog10_d"] = dlog;
  cols["dN_dlog10_d_per_m3"] = dndlog;

  std::map<std::string, double> meta;
  meta["abscissa"] = data_->additional[0][point];
  meta["T"] = T;
  meta["P"] = P;
  meta["rho"] = rho;
  meta["MW"] = MWmix;
  return std::make_pair(cols, meta);
}

std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, double>>
Soot::ParticleSizeDistribution(const double local_value, const std::string& particle_type,
                               const std::string& diameter_type, const int min_section,
                               const double mobility_exponent,
                               const double merge_tol) const {
  if (data_ == nullptr || !data_->soot_available_)
    throw std::invalid_argument("The mechanism has no <SootProperties> block");
  if (particle_type != "all" && particle_type != "primary")
    throw std::invalid_argument("particle_type must be \"all\" or \"primary\"");
  if (diameter_type != "dmob" && diameter_type != "dpp" && diameter_type != "dcol")
    throw std::invalid_argument("diameter_type must be \"dmob\", \"dpp\" or \"dcol\"");

  const unsigned int point = AbscissaIndex(local_value);
  const double numpp_tol = 1e-6;
  const bool primary_only = (particle_type == "primary");

  std::vector<unsigned int> bins;
  std::vector<double> size_m;
  for (unsigned int b = 0; b < data_->soot_number_of_bins_; b++) {
    const double npp = data_->soot_bin_numpp_[b];
    if (primary_only) {
      if (std::fabs(npp - 1.0) > numpp_tol) continue;
    } else {
      if (data_->soot_bin_section_[b] < min_section) continue;
    }

    bins.push_back(b);
    if (diameter_type == "dpp") {
      size_m.push_back(data_->soot_bin_dpp_[b]);
    } else if (diameter_type == "dcol") {
      size_m.push_back(data_->soot_bin_dcol_[b]);
    } else {  // "dmob": numPP <= 0 (nascent) counts as one spherule
      double npp_eff = npp;
      if (npp_eff <= 0.0) npp_eff = 1.0;
      size_m.push_back(data_->soot_bin_dpp_[b] * std::pow(npp_eff, mobility_exponent));
    }
  }
  if (bins.empty()) {
    if (primary_only) throw std::invalid_argument("No BIN with numPP == 1");
    throw std::invalid_argument("No BIN with Bin_section >= min_section");
  }

  return BuildDistribution(bins, size_m, point, merge_tol);
}
