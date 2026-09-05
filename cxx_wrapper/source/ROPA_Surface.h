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

#ifndef ROPA_SURFACE_H
#define ROPA_SURFACE_H

#include "ProfilesDatabase.h"
#include "ROPA.h"  // Parent class

// [LG] DEV note: I am unsure what to do with the heterogeneous_reactions variable and
// SetROPAPhase function. On one side, setting it once and having it internally makes it
// easier to use. On the other side, using it as an input parameter is generally better
// for readability and requires no override for parent functions. Right now, the
// SetROPAPhase function is redundant as it is passed and updated within the other
// functions.

class ROPA_Surface : public ROPA {
 public:
  ROPA_Surface();

  void SetROPAPhase(const bool heterogeneous_reactions);

  void RateOfProductionAnalysis(const unsigned int number_of_reactions,
                                const bool heterogeneous_reactions);

  void GetReactionRates(std::vector<unsigned int> reaction_indices, const bool sum_rates,
                        const bool heterogeneous_reactions);

  void GetFormationRates(std::string specie, std::string units, std::string type,
                         const bool heterogeneous_reactions);

  // TODO: ROPA1D for UnsteadyPFR
  // void RateOfProductionAnalysis2D(const unsigned int number_of_reactions, const double
  // local_x,
  //                                 const double local_z, const double region_low_x,
  //                                 const double region_up_x, const double region_low_z,
  //                                 const double region_up_z);

  // inline const std::vector<int>& indexFirstName() const { return indexFirstName_; };

  // inline const std::vector<int>& indexSecondName() const { return indexSecondName_; };

  // inline const std::vector<double>& computedThickness() const { return
  // computedThickness_; };

  // inline const std::vector<double>& computedLabel() const { return computedLabel_; };

 protected:
  bool heterogeneous_reactions_;
};

#include "ROPA_Surface.hpp"
#endif  // ROPA_H
