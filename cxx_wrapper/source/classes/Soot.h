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

#include "core/PostProcessorCore.h"

// Soot post-processing of the mechanism's <SootProperties> optional leaf

class Soot {
 public:
  Soot();
  void SetResults(PostProcessorCore* data);

  // True when kinetics.xml actually carried a <SootProperties> block.
  // Basically true for any CRECK-compiled model 2024+
  inline bool sootAvailable() const {
    return data_ != nullptr && data_->gasKinetics()->soot_available_;
  }

  inline unsigned int numberOfBins() const {
    return data_->gasKinetics()->soot_number_of_bins_;
  }

  // Properties (as wrote into the kinetics.xml from BinProperties.txt at mech
  // compile-time)
  inline const std::vector<unsigned int>& index() const {
    return data_->gasKinetics()->soot_bin_index_;
  }
  inline const std::vector<int>& section() const {
    return data_->gasKinetics()->soot_bin_section_;
  }
  inline const std::vector<double>& nc() const {
    return data_->gasKinetics()->soot_bin_nc_;
  }
  inline const std::vector<double>& nh() const {
    return data_->gasKinetics()->soot_bin_nh_;
  }
  inline const std::vector<double>& no() const {
    return data_->gasKinetics()->soot_bin_no_;
  }
  inline const std::vector<double>& htoc() const {
    return data_->gasKinetics()->soot_bin_htoc_;
  }
  inline const std::vector<double>& mw() const {
    return data_->gasKinetics()->soot_bin_mw_;
  }
  inline const std::vector<double>& density() const {
    return data_->gasKinetics()->soot_bin_density_;
  }
  inline const std::vector<double>& volume() const {
    return data_->gasKinetics()->soot_bin_volume_;
  }
  inline const std::vector<double>& mass() const {
    return data_->gasKinetics()->soot_bin_mass_;
  }
  inline const std::vector<double>& numpp() const {
    return data_->gasKinetics()->soot_bin_numpp_;
  }
  inline const std::vector<double>& dsph() const {
    return data_->gasKinetics()->soot_bin_dsph_;
  }
  inline const std::vector<double>& dcol() const {
    return data_->gasKinetics()->soot_bin_dcol_;
  }
  inline const std::vector<double>& dpp() const {
    return data_->gasKinetics()->soot_bin_dpp_;
  }
  inline const std::vector<double>& df() const {
    return data_->gasKinetics()->soot_bin_df_;
  }

  // Soot Particle Size Distribution at one abscissa location
  // abscissa = additional[0] (time for a reactor, spatial coordinate for a flame)
  //
  //   particle_type = "all"     -> all Soot particles (e.g. liq, pp, aggs)
  //                   "primary" -> pp only (numPP == 1)
  //   both restricted to BINs with Bin_section >= min_section
  //   diameter_type = "dmob"    -> mobility, dm = Dpp * numPP^mobility_exponent
  //                   "dpp"     -> primary-particle diameter
  //                   "dcol"    -> collision diameter
  //                   "dva"     -> volume-equivalent sphere diameter (Bin_dsph)
  //
  // Returns {columns, meta}:
  //   columns: name -> one per-section vector - x, x_min, x_max (diameters in
  //            nm), N_per_m3, n_bins, Dlog10_x, dN_dlog10_x_per_m3
  //   meta:    the scalar gas state actually used - abscissa, T, P, rho, MW
  //
  // Per-bin number density computed from ProfilesDatabase:
  // N_i = omega_i * rho / MW_i * 1000 * N_A  [#/m3].
  // Near-equal diameters are merged into fixed sections
  // (representative = geometric mean) before dN/dlog10(d) is formed.
  std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, double>>
  ParticleSizeDistribution(double local_value, const std::string& particle_type,
                           const std::string& diameter_type, int min_section,
                           double mobility_exponent, double merge_tol) const;

  // --- soot specific surface area (SSA) --------------------------------------
  // Per BIN: SSA_i = pi * Dpp_i^2 * numPP_i / Bin_mass_i  [m2/kg]
  // Only BINs with Bin_section >= min_section (NaN otherwise); among those, liq.
  // particles with numPP = 0 count as 1.
  std::vector<double> SpecificSurfaceArea(int min_section) const;

  // --- Takes the weighted average of any per-BIN property along the abscissa --------
  // Over BINs with Bin_section >= min_section and a finite property value:
  //   weighting = "mass"   -> w_i = N_i m_i   (kg/m3),  N_i = C_i * N_A --- used for fv
  //   and SSA
  //               "carbon" -> w_i = C_i nC_i  (mol of C/m3); for p = H/C --- used for H/C
  //               ratio
  // `property` holds one value per BIN, in BIN order (e.g. htoc(), dpp()).
  // Returns columns abscissa, mean (NaN where there is no soot), mass_kg_per_m3.
  std::map<std::string, std::vector<double>> AveragedProfile(
      const std::vector<double>& property, int min_section,
      const std::string& weighting) const;

 protected:
  PostProcessorCore* data_;

  unsigned int AbscissaIndex(double local_value) const;
  std::vector<double> BinMassFractionToMolarConcentration(unsigned int point) const;
  std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, double>>
  BuildDistribution(const std::vector<unsigned int>& bins,
                    const std::vector<double>& coord, unsigned int point,
                    double merge_tol) const;
};

#include "Soot.hpp"
#endif  // SOOT_H
