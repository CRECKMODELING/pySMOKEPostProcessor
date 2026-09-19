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

#ifndef POST_PROCESSOR_FLUX_MAP
#define POST_PROCESSOR_FLUX_MAP

namespace pySMOKEPostProcessor {
class PostProcessorFluxMap : public OpenSMOKE::FluxAnalysisMap {
 public:
  PostProcessorFluxMap(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML,
                       OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapXML);

  // Flux analysis for heterogeneous reactions requires addition of flux analysis map in OpenSMOKE previous to postprocessing.
  // As of right now, it is only developed for gas-phase reactions.

  void ComputeFluxAnalysis();
  void ComputeValues(const unsigned int index_j, std::vector<unsigned int>& local_indices,
                     std::vector<double>& local_thickness, std::vector<double>& local_normal_fluxes,
                     std::vector<double>& local_fluxes);

  std::vector<int> IndexFirstName;
  std::vector<int> IndexSecondName;
  std::vector<double> ComputedLabelValue;
  std::vector<double> ComputedThicknessValue;
};
}  // namespace pySMOKEPostProcessor

#include "PostProcessorFluxMap.hpp"
#endif  // POST_PROCESSOR_FLUX_MAP
