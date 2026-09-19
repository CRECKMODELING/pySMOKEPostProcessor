/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Authors: Timoteo Dinelli <timoteo.dinelli@polimi.it>                  |
|            Edoardo Ramalli <edoardo.ramalli@polimi.it>                  |
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

#ifndef ROPA_H
#define ROPA_H

#include "core/PostProcessorCore.h"
#include "core/Utilities.h"

class ROPA {
 public:
  ROPA();

  void SetResults(PostProcessorCore* data);

  void RateOfProductionAnalysis(const unsigned int number_of_reactions);

  void RateOfProductionAnalysis2D(const unsigned int number_of_reactions, const double local_x,
                                  const double local_z, const double region_low_x,
                                  const double region_up_x, const double region_low_z,
                                  const double region_up_z);

  void FluxAnalysis();

  void GetReactionRates(std::vector<unsigned int> reaction_indices, const bool sum_rates);

  void GetFormationRates(std::string specie, std::string units, std::string type);

  void SetKineticFolder(const std::string kineticFolder);

  void SetOutputFolder(const std::string outputFolder);

  void SetROPAType(const std::string type);

  void SetSpecies(const std::string species);

  void SetLocalValue(double localValue);

  void SetLowerBound(double lowerBound);

  void SetUpperBound(double upperBound);

  void SetElement(const std::string element);

  void SetThickness(const std::string thickness);

  void SetFluxAnalysisType(const std::string type);

  void SetWidth(const int width);

  void SetDepth(const int depth);

  void SetThreshold(const double threshold);

  void SetThicknessLogScale(bool thicknesslogscale);

  void SetLabelType(std::string type);

  inline const std::vector<unsigned int>& reactions() const { return reactions_; };

  inline const std::vector<double>& coefficients() const { return coefficients_; };

  inline const std::vector<int>& indexFirstName() const { return indexFirstName_; };

  inline const std::vector<int>& indexSecondName() const { return indexSecondName_; };

  inline const std::vector<double>& computedThickness() const { return computedThickness_; };

  inline const std::vector<double>& computedLabel() const { return computedLabel_; };

  inline const std::vector<double>& formationRates() const { return formationRates_; };

  inline const std::vector<std::vector<double>>& reactionRates() const { return reactionRates_; };

  inline const std::vector<double>& sumOfRates() const { return sumOfRates_; };

 protected:
  PostProcessorCore* data_;
  std::vector<unsigned int> indices_coarse_reactions_;
  std::vector<std::string> string_list_reactions;

  std::string ropaType_;
  std::string kineticFolder_;
  std::string outputFolder_;
  std::string species_;

  double localValue_;
  double upperBound_;
  double lowerBound_;
  bool speciesIsSelected;

  std::string element_;
  std::string thickness_;
  std::string flux_type_;
  int width_;
  int depth_;
  double threshold_;
  bool thicknesslogscale_;
  std::string label_type_;

  std::vector<unsigned int> reactions_;
  std::vector<double> coefficients_;

  std::vector<int> indexFirstName_;
  std::vector<int> indexSecondName_;
  std::vector<double> computedThickness_;
  std::vector<double> computedLabel_;

  std::vector<double> formationRates_;
  std::vector<std::vector<double>> reactionRates_;
  std::vector<double> sumOfRates_;
};

#include "ROPA.hpp"
#endif  // ROPA_H
