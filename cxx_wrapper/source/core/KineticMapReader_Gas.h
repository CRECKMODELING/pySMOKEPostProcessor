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

#ifndef KINETICMAPREADER_GAS_H
#define KINETICMAPREADER_GAS_H

#include "KineticMapReaderBase.h"

// Gas-phase kinetics reader: kinetics.xml + reaction_names.xml. Owns the
// OpenSMOKEpp gas thermo/kinetics maps (raw pointers, never deleted).
// Also owns the mechanism-wide <SootProperties> block
// and the ReactionsAssociatedToSpecies/isReactantProduct helpers,
// moved here from ProfilesDatabase since they only ever touched the gas kinetics map.
// Note: these are helpers related to the FluxAnalysis.
//       if the FluxAnalysis is extended to other phases, this has to be updated.
class KineticMapReader_Gas : public KineticMapReaderBase {
 public:
  KineticMapReader_Gas();
  ~KineticMapReader_Gas();

  bool Load(const std::string& folder_name);

  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMap() const { return thermodynamicsMap_; }
  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMap() const { return kineticsMap_; }

  double MW(unsigned int index) const override;

  // In SensitivityAnalysis, FallOff/CAB Reactions also undergo sensitivity on
  // kInf to get information on the fall-off behavior (e.g. if kInf is high in sensi,
  // the reaction is not in fall-off and vice-versa).
  // The index of these kInfs goes beyond the number of reactions, and therefore
  // requires its own dedicated parsing
  std::string FormattedReactionNameFromIndex(unsigned int index) const override;

  void ReactionsAssociatedToSpecies(const unsigned int index, std::vector<unsigned int>& indices) const;
  void IsReactantProduct(const unsigned int reaction_index, double& netStoichiometry) const;

  const Eigen::SparseMatrix<double>& StoichiometricMatrixReactants() const;
  const Eigen::SparseMatrix<double>& StoichiometricMatrixProducts() const;

  // <SootProperties> (kinetics.xml). Empty when the mechanism carries no soot bins.
  bool soot_available_;
  unsigned int soot_number_of_bins_;
  std::vector<unsigned int> soot_bin_index_;
  std::vector<int> soot_bin_section_;
  std::vector<double> soot_bin_nc_;
  std::vector<double> soot_bin_nh_;
  std::vector<double> soot_bin_no_;
  std::vector<double> soot_bin_htoc_;
  std::vector<double> soot_bin_mw_;
  std::vector<double> soot_bin_density_;
  std::vector<double> soot_bin_volume_;
  std::vector<double> soot_bin_mass_;
  std::vector<double> soot_bin_numpp_;
  std::vector<double> soot_bin_dsph_;
  std::vector<double> soot_bin_dcol_;
  std::vector<double> soot_bin_dpp_;
  std::vector<double> soot_bin_df_;

 private:
  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMap_;
  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMap_;

  // 1-based reaction numbers, parsed from <FallOff>/<CABR> (opensmoke.Kinetics.),
  // used only by FormattedReactionNameFromIndex's redirect above.
  std::vector<unsigned int> indices_of_falloff_reactions_;
  std::vector<unsigned int> indices_of_cabr_reactions_;

  void ReadSootProperties(const boost::property_tree::ptree& mechanism_ptree);
  void ReadFallOffAndCabrIndices(const boost::property_tree::ptree& mechanism_ptree);
};

#include "KineticMapReader_Gas.hpp"
#endif  // KINETICMAPREADER_GAS_H
