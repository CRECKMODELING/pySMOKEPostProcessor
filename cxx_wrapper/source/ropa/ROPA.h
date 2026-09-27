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

// RateOfProductionAnalysis/GetReactionRates/GetFormationRates take
// heterogeneous_reactions: bool = false and internally dispatch to whichever phase-2
// kinetics is actually loaded (data_->phase2Kind()) when true - Surface or Liquid today,
// Solid once it gets its own ROPA support. This replaces the former ROPA_Surface and
// ROPA_Liquid subclasses (which re-implemented these same three methods, in
// ROPA_Surface's case with a heterogeneous_reactions parameter that just branched to
// ROPA::<method>() when false) - same unification SensitivityReader::Prepare(Phase)
// already did for gas/surface sensitivity. A mechanism has at most one phase-2 kinetics
// loaded at a time (see PostProcessorCore.h), so a single bool is enough to mean
// "whichever one that is" - each dispatch site below still checks data_->phase2Kind()
// explicitly and throws a clear error rather than guessing. The *_Surface/*_Liquid
// bodies themselves live in ROPA_Surface.hpp/ROPA_Liquid.hpp, split out purely to keep
// this file a readable size.
class ROPA {
 public:
  ROPA();

  void SetResults(PostProcessorCore* data);

  void RateOfProductionAnalysis(const unsigned int number_of_reactions,
                                const bool heterogeneous_reactions = false);

  void RateOfProductionAnalysis2D(const unsigned int number_of_reactions,
                                  const double local_x,      const double local_z,
                                  const double region_low_x, const double region_up_x,
                                  const double region_low_z, const double region_up_z);

  void FluxAnalysis();

  void GetReactionRates(std::vector<unsigned int> reaction_indices, const bool sum_rates,
                        const bool heterogeneous_reactions = false);

  void GetFormationRates(std::string specie, std::string units, std::string type,
                         const bool heterogeneous_reactions = false);

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
  // ROPA on heterogeneous kinetics belongs to the same class, but is implemented
  // in different files (ROPA_Surface and ROPA_Liquid) to maintain clarity.
  // Calls to ROPA when heterogeneous_reactions=true are redirected.
  void RateOfProductionAnalysis_Surface(const unsigned int number_of_reactions);
  void GetReactionRates_Surface(const std::vector<unsigned int>& reaction_indices,
                                const bool sum_rates);
  void GetFormationRates_Surface(const std::string& specie, const std::string& units,
                                 const std::string& type);


  // To be tested, do not take this for granted
  void RateOfProductionAnalysis_Liquid(const unsigned int number_of_reactions);
  void GetReactionRates_Liquid(const std::vector<unsigned int>& reaction_indices,
                               const bool sum_rates);
  void GetFormationRates_Liquid(const std::string& specie, const std::string& units,
                                const std::string& type);

  PostProcessorCore* data_;
  std::vector<unsigned int> indices_coarse_reactions_;
  std::vector<std::string> string_list_reactions;

  std::string ropaType_;
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
#include "ROPA_Liquid.hpp"  // Class is shared, put in a different file for clarity
#include "ROPA_Surface.hpp" // Class is shared, put in a different file for clarity
#endif  // ROPA_H
