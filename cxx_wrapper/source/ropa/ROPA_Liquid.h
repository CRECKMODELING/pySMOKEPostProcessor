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

#ifndef ROPA_LIQUID_H
#define ROPA_LIQUID_H

#include "core/PostProcessorCore.h"
#include "ROPA.h"  // Parent class

// Liquid-phase ROPA. Unlike ROPA_Surface (a genuinely different, composite
// gas+site+bulk state vector), a liquid mechanism's Output.xml is assumed to
// carry the same single-phase <mass-fractions>/<profiles> shape gas does -
// confirmed for kinetics.liquid.xml against a real mechanism, but no liquid
// Output.xml fixture exists anywhere to verify this against, so this class is
// structurally verified only (compiles, wired correctly), not exercised
// end-to-end. Because of that shape equivalence, RateOfProductionAnalysis/
// GetReactionRates/GetFormationRates below are ROPA's own single-phase logic
// with kineticsMapXML/thermodynamicsMapXML swapped for
// kineticsMapLiquidXML/thermodynamicsMapLiquidXML - not ROPA_Surface's
// composite-state logic, which doesn't apply here. FluxAnalysis stays
// inherited from ROPA unchanged (gas-only, same limitation ROPA_Surface
// already has - see its own header comment).
class ROPA_Liquid : public ROPA {
 public:
  ROPA_Liquid();

  // Also checks phase2Kind() == Liquid, so a mismatched PostProcessorCore
  // (e.g. surface, or gas-only) raises a clear error instead of this class's
  // methods segfaulting on a null kineticsMapLiquidXML.
  void SetResults(PostProcessorCore* data);

  void RateOfProductionAnalysis(const unsigned int number_of_reactions);

  void GetReactionRates(std::vector<unsigned int> reaction_indices, const bool sum_rates);

  void GetFormationRates(std::string specie, std::string units, std::string type);
};

#include "ROPA_Liquid.hpp"
#endif  // ROPA_LIQUID_H
