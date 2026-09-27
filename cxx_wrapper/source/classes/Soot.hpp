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
constexpr double PI = 3.14159265358979323846;
}  // namespace

Soot::Soot() { data_ = nullptr; }

void Soot::SetResults(PostProcessorCore* data) {
  if (!data->HasKinetics()) {
    throw std::invalid_argument(
        "Soot::SetResults: no kinetics mechanism loaded on this PostProcessorCore");
  }
  if (!data->HasOutput()) {
    throw std::invalid_argument(
        "Soot::SetResults: no simulation output loaded on this PostProcessorCore");
  }
  data_ = data;
}

// local value to print the PSD from
unsigned int Soot::AbscissaIndex(double local_value) const {
  const std::vector<double>& abscissa = data_->additional[0];
  const unsigned int npts = static_cast<unsigned int>(data_->number_of_abscissas_);
  for (unsigned int j = 0; j < npts; j++)
    if (abscissa[j] >= local_value) return j;
  if (npts == 0) return 0;
  return npts - 1;
}

std::vector<double> Soot::BinMassFractionToMolarConcentration(
    const unsigned int point) const {
  const unsigned int nbins = data_->gasKinetics()->soot_number_of_bins_;
  const double rho = data_->AdditionalProfile("density")[point];  // kg/m3
  std::vector<double> C(nbins, 0.0);
  for (unsigned int b = 0; b < nbins; b++) {
    const unsigned int sp =
        data_->gasKinetics()->soot_bin_index_[b];  // gas species index
    if (sp < data_->omega.size() && sp < data_->mw_species_.size() &&
        data_->mw_species_[sp] > 0.0)
      C[b] = data_->omega[sp][point] * rho * 1000.0 / data_->mw_species_[sp];
  }
  return C;
}

// Specific surface area of every BIN [m2/kg]: pi * Dpp^2 * numPP / Bin_mass.
// Only BINs with Bin_section >= min_section are soot, the rest stay NaN; among
// those, numPP <= 0 (nascent/liquid-like BINs) counts as numPP = 1.
std::vector<double> Soot::SpecificSurfaceArea(const int min_section) const {
  const unsigned int nbins = data_->gasKinetics()->soot_number_of_bins_;
  const auto& dpp = data_->gasKinetics()->soot_bin_dpp_;
  const auto& npp = data_->gasKinetics()->soot_bin_numpp_;
  const auto& mass = data_->gasKinetics()->soot_bin_mass_;  // kg per particle
  const auto& section = data_->gasKinetics()->soot_bin_section_;
  std::vector<double> ssa(nbins, std::numeric_limits<double>::quiet_NaN());
  for (unsigned int b = 0; b < nbins; b++) {
    if (section[b] < min_section) continue;
    double npp_eff = npp[b];
    if (npp_eff <= 0.0) npp_eff = 1.0;
    if (mass[b] > 0.0) ssa[b] = PI * dpp[b] * dpp[b] * npp_eff / mass[b];
  }
  return ssa;
}

// Shared builder for any distribution: concentration per selected bin at
// `point`, then merge near-equal coordinates into fixed sections.
// `coord` is the coordinate of each entry in `bins`, already in nm.
std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, double>>
Soot::BuildDistribution(const std::vector<unsigned int>& bins,
                        const std::vector<double>& coord, const unsigned int point,
                        const double merge_tol) const {
  const double T = data_->additional[data_->index_T][point];
  const double P = data_->additional[data_->index_P][point];
  const double MWmix = data_->additional[data_->index_MW][point];
  const double rho = data_->AdditionalProfile("density")[point];  // kg/m3
  const auto& bin_index = data_->gasKinetics()->soot_bin_index_;
  const auto& bin_mass = data_->gasKinetics()->soot_bin_mass_;  // kg per particle

  // N_i = rho * omega_i / Bin_mass_i [#/m3]. Uses the BIN particle mass rather than
  // MW_i / N_A, so N matches OpenSMOKE's PolimiSoot N(tot); the two differ by ~6e-5
  // because of the digits carried by the BinProperties.
  const unsigned int n = static_cast<unsigned int>(bins.size());
  std::vector<double> d(n), N(n, 0.0);
  for (unsigned int k = 0; k < n; k++) {
    d[k] = coord[k];
    const unsigned int sp = bin_index[bins[k]];
    if (sp < data_->omega.size() && bin_mass[bins[k]] > 0.0)
      N[k] = rho * data_->omega[sp][point] / bin_mass[bins[k]];
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
    col_d.push_back(std::pow(10.0, log_sum / count));
    col_dmin.push_back(dmin);
    col_dmax.push_back(dmax);
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

  // Centred log spacing on the sorted representative coordinates: one-sided at
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
  cols["x"] = col_d;
  cols["x_min"] = col_dmin;
  cols["x_max"] = col_dmax;
  cols["N_per_m3"] = col_N;
  cols["n_bins"] = col_nbins;
  cols["Dlog10_x"] = dlog;
  cols["dN_dlog10_x_per_m3"] = dndlog;

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
  if (data_ == nullptr || !data_->gasKinetics()->soot_available_)
    throw std::invalid_argument("The mechanism has no <SootProperties> block");
  if (particle_type != "all" && particle_type != "primary")
    throw std::invalid_argument("particle_type must be \"all\" or \"primary\"");
  if (diameter_type != "dmob" && diameter_type != "dpp" && diameter_type != "dcol" &&
      diameter_type != "dva")
    throw std::invalid_argument(
        "diameter_type must be \"dmob\", \"dpp\", \"dcol\" or \"dva\"");

  const unsigned int point = AbscissaIndex(local_value);
  const double numpp_tol = 1e-6;
  const bool primary_only = (particle_type == "primary");

  std::vector<unsigned int> bins;
  std::vector<double> size_m;
  for (unsigned int b = 0; b < data_->gasKinetics()->soot_number_of_bins_; b++) {
    const double npp = data_->gasKinetics()->soot_bin_numpp_[b];
    if (data_->gasKinetics()->soot_bin_section_[b] < min_section) continue;
    if (primary_only && std::fabs(npp - 1.0) > numpp_tol) continue;

    bins.push_back(b);
    if (diameter_type == "dpp") {
      size_m.push_back(data_->gasKinetics()->soot_bin_dpp_[b]);
    } else if (diameter_type == "dcol") {
      size_m.push_back(data_->gasKinetics()->soot_bin_dcol_[b]);
    } else if (diameter_type == "dva") {  // volume-equivalent sphere diameter (Bin_dsph)
      size_m.push_back(data_->gasKinetics()->soot_bin_dsph_[b]);
    } else {  // "dmob": numPP <= 0 (nascent) counts as one spherule
      double npp_eff = npp;
      if (npp_eff <= 0.0) npp_eff = 1.0;
      size_m.push_back(data_->gasKinetics()->soot_bin_dpp_[b] *
                       std::pow(npp_eff, mobility_exponent));
    }
  }
  if (bins.empty()) {
    if (primary_only)
      throw std::invalid_argument(
          "No BIN with numPP == 1 and Bin_section >= min_section");
    throw std::invalid_argument("No BIN with Bin_section >= min_section");
  }

  std::vector<double> size_nm(size_m.size());
  for (size_t k = 0; k < size_m.size(); k++) size_nm[k] = size_m[k] * 1e9;
  return BuildDistribution(bins, size_nm, point, merge_tol);
}

std::map<std::string, std::vector<double>> Soot::AveragedProfile(
    const std::vector<double>& property, const int min_section,
    const std::string& weighting) const {
  if (data_ == nullptr || !data_->gasKinetics()->soot_available_)
    throw std::invalid_argument("The mechanism has no <SootProperties> block");
  if (weighting != "mass" && weighting != "carbon")
    throw std::invalid_argument("weighting must be \"mass\" or \"carbon\"");

  const unsigned int nbins = data_->gasKinetics()->soot_number_of_bins_;
  if (property.size() != nbins)
    throw std::invalid_argument("property must hold one value per BIN");
  const auto& bin_index = data_->gasKinetics()->soot_bin_index_;
  const auto& bin_nc = data_->gasKinetics()->soot_bin_nc_;
  const auto& section = data_->gasKinetics()->soot_bin_section_;
  const std::vector<double>& rho = data_->AdditionalProfile("density");  // kg/m3
  const bool carbon_weighted = (weighting == "carbon");

  std::vector<unsigned int> bins;
  for (unsigned int b = 0; b < nbins; b++)
    if (section[b] >= min_section && std::isfinite(property[b]) &&
        bin_index[b] < data_->omega.size())
      bins.push_back(b);
  if (bins.empty()) throw std::invalid_argument("No BIN with Bin_section >= min_section");

  const unsigned int npts = static_cast<unsigned int>(data_->number_of_abscissas_);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x(npts), mean(npts, nan), mass(npts, 0.0);
  for (unsigned int j = 0; j < npts; j++) {
    x[j] = data_->additional[0][j];
    const std::vector<double> C = BinMassFractionToMolarConcentration(j);
    double weight_sum = 0.0;
    double weighted_sum = 0.0;
    for (const unsigned int b : bins) {
      const double m_b = rho[j] * data_->omega[bin_index[b]][j];  // kg/m3
      mass[j] += m_b;
      double w_b = m_b;
      if (carbon_weighted) w_b = C[b] * bin_nc[b];  // mol of C/m3
      weight_sum += w_b;
      weighted_sum += w_b * property[b];
    }
    if (weight_sum > 0.0) mean[j] = weighted_sum / weight_sum;
  }

  std::map<std::string, std::vector<double>> cols;
  cols["abscissa"] = x;
  cols["mean"] = mean;
  cols["mass_kg_per_m3"] = mass;
  return cols;
}
