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

#ifndef KINETICMAPREADER_SOLID_H
#define KINETICMAPREADER_SOLID_H

#include "KineticMapReaderBase.h"

// Solid-phase kinetics reader: kinetics.solid.xml 
// Two hard asymmetries vs. every other phase:
//  - No reaction-names file exists for solid mechanisms, so reaction_names_
//    is synthesized ("Reaction #1", "Reaction #2", ...) instead of read from
//    a file - a missing *cosmetic* file degrades gracefully rather than
//    blocking construction, since nothing downstream needs the human name.
//  - OpenSMOKEpp's KineticsMap_Solid_CHEMKIN has no RateOfProductionAnalysis
//    primitive at all (a hard OpenSMOKEpp limitation, out of scope to patch),
//    so has_ropa_ is always false and there is deliberately no ROPA_Solid
//    widget anywhere in this codebase - not a runtime-guarded feature, an
//    absent one.
// @Riccardo a te l'onore
class KineticMapReader_Solid : public KineticMapReaderBase {
 public:
  KineticMapReader_Solid();
  ~KineticMapReader_Solid();

  bool Load(const std::string& folder_name);

  OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermodynamicsMap() const { return thermodynamicsMap_; }
  OpenSMOKE::KineticsMap_Solid_CHEMKIN* kineticsMap() const { return kineticsMap_; }

  double MW(unsigned int index) const override;

 private:
  OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermodynamicsMap_;
  OpenSMOKE::KineticsMap_Solid_CHEMKIN* kineticsMap_;
};

#include "KineticMapReader_Solid.hpp"
#endif  // KINETICMAPREADER_SOLID_H
