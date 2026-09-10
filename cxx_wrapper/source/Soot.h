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

#ifndef SOOT_H
#define SOOT_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ProfilesDatabase.h"

// Soot post-processing of the mechanism's <SootProperties> optional leaf
//
// Currently there are getters for all BinProperties + ParticleSizeDistribution
// function working for both primary and overall PSD

class Soot {
 public:
  Soot();
  void SetDatabase(ProfilesDatabase* data);

  // True when kinetics.xml actually carried a <SootProperties> block. When false,
  // every vector below is empty and numberOfBins() is 0.
  inline bool sootAvailable() const { return data_ != nullptr && data_->soot_available_; }
  
  inline unsigned int numberOfBins() const { return data_->soot_number_of_bins_; }

  inline const std::vector<unsigned int>& index() const { return data_->soot_bin_index_; }
  inline const std::vector<int>& section() const { return data_->soot_bin_section_; }
  inline const std::vector<double>& nc() const { return data_->soot_bin_nc_; }
  inline const std::vector<double>& nh() const { return data_->soot_bin_nh_; }
  inline const std::vector<double>& no() const { return data_->soot_bin_no_; }
  inline const std::vector<double>& htoc() const { return data_->soot_bin_htoc_; }
  inline const std::vector<double>& mw() const { return data_->soot_bin_mw_; }
  inline const std::vector<double>& density() const { return data_->soot_bin_density_; }
  inline const std::vector<double>& volume() const { return data_->soot_bin_volume_; }
  inline const std::vector<double>& mass() const { return data_->soot_bin_mass_; }
  inline const std::vector<double>& numpp() const { return data_->soot_bin_numpp_; }
  inline const std::vector<double>& dsph() const { return data_->soot_bin_dsph_; }
  inline const std::vector<double>& dcol() const { return data_->soot_bin_dcol_; }
  inline const std::vector<double>& dpp() const { return data_->soot_bin_dpp_; }
  inline const std::vector<double>& df() const { return data_->soot_bin_df_; }

  // --- soot particle size distribution at one abscissa location -------------
  // local_value selects the profile point the same way local ROPA does. 
  // (additional[0] - time for a reactor, a spatial coordinate for a flame)
  //
  //   particle_type = "all"     -> every BIN with Bin_section >= min_section
  //                                (numPP <= 0 counts as one spherule)
  //                   "primary" -> free primary particles only, numPP == 1
  //                                (min_section is ignored)
  //   diameter_type = "dmob"    -> mobility, dm = Dpp * numPP^mobility_exponent
  //                   "dpp"     -> primary-particle diameter
  //                   "dcol"    -> collision diameter
  //
  // Returns {columns, meta}:
  //   columns: name -> one per-section vector - d_nm, d_min_nm, d_max_nm,
  //            N_per_m3, n_bins, Dlog10_d, dN_dlog10_d_per_m3
  //   meta:    the scalar gas state actually used - abscissa, T, P, rho, MW
  //
  // Per-bin number density is the mass-fraction route, needing nothing beyond
  // ProfilesDatabase:  N_i = omega_i * rho / MW_i * 1000 * N_A  [#/m3], with
  // rho = P*MW_mix/(R*T). Near-equal diameters are merged into fixed sections
  // (representative = geometric mean) before dN/dlog10(d) is formed.
  std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, double>>
  ParticleSizeDistribution(double local_value, const std::string& particle_type,
                           const std::string& diameter_type, int min_section,
                           double mobility_exponent, double merge_tol) const;

 protected:
  ProfilesDatabase* data_;

  unsigned int AbscissaIndex(double local_value) const;
  std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, double>>
  BuildDistribution(const std::vector<unsigned int>& bins,
                    const std::vector<double>& size_m, unsigned int point,
                    double merge_tol) const;
};

#include "Soot.hpp"
#endif  // SOOT_H
