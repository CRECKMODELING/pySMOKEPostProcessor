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

// Post-processing of the mechanism's <SootProperties> optional leaf
class Soot {
 public:
  Soot();
  void SetResults(PostProcessorCore* data);

  // True when kinetics.xml actually carries a <SootProperties> block.
  // Basically true for any CRECK-compiled model 2024+
  inline bool sootAvailable() const {
    return data_ != nullptr && data_->gasKinetics()->soot_available_;
  }

  inline unsigned int numberOfBins() const {
    return data_->gasKinetics()->soot_number_of_bins_;
  }

  // Properties as wrote into the kinetics.xml from BinProperties.txt
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
  // abscissa = additional[0]
  //
  //   particle_type = "all"     -> every Soot BIN counted as one particle each
  //                                (aggregates and singles alike)
  //                   "primary" -> the primary-particle size distribution (PPSD):
  //                                every BIN (single particles AND aggregates)
  //                                contributes the primary particles it is made of.
  //                                Only diameter_type = "dpp" is meaningful for 
  //                                this population; anything else throws.
  //   both restricted to BINs with Bin_section >= min_section.
  //   By default, min_section = 1 for PSD, min_section = 5 for PPSD
  //   diameter_type = "dmob"    -> mobility diameter, dm = Dpp * numPP^exponent;
  //                                the exact correlation (prefactor + exponent) is
  //                                picked by correlation_name (only used for this
  //                                diameter_type).
  //                                  "Kelesidis" -> dm = Dpp * numPP^0.45
  //                                    (Kelesidis et al. 2017, Carbon,
  //                                    10.1016/j.carbon.2017.06.004)
  //                                  "Sorensen"  -> dm = Dpp * numPP^0.465
  //                                    (Sorensen 2011, Aerosol Sci. Technol.,
  //                                    10.1080/02786826.2011.560909)
  //                                  "Rissler"   -> dm = 0.794 * Dpp * numPP^0.51
  //                                    (Rissler et al. 2013, Aerosol Sci. Technol.,
  //                                    10.1080/02786826.2013.791381)
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
  // N_i = omega_i * rho / Bin_mass_i  [#/m3] (particle_type == "all"), or that
  // times Bin_numpp_i (particle_type == "primary", see above).
  // Near-equal diameters are merged into fixed sections
  // (representative = geometric mean) before dN/dlog10(d) is formed. For
  // particle_type == "primary" the physically meaningful quantity is N_per_m3
  // itself (a primary-particle count).
  std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, double>>
  ParticleSizeDistribution(double local_value, const std::string& particle_type,
                           const std::string& diameter_type, int min_section,
                           const std::string& correlation_name, double merge_tol) const;

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
                    double merge_tol, const std::vector<double>& weight = {}) const;
};

#include "Soot.hpp"
#endif  // SOOT_H
