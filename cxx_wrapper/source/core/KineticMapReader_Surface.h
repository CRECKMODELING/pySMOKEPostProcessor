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

#ifndef KINETICMAPREADER_SURFACE_H
#define KINETICMAPREADER_SURFACE_H

#include "KineticMapReaderBase.h"

// Surface-phase kinetics reader: kinetics.surface.xml + surface_reaction_names.xml.
// <ReactionClasses> is supported, same as gas.
class KineticMapReader_Surface : public KineticMapReaderBase {
 public:
  KineticMapReader_Surface();
  ~KineticMapReader_Surface();

  bool Load(const std::string& folder_name);

  OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN* thermodynamicsMap() const { return thermodynamicsMap_; }
  OpenSMOKE::KineticsMap_Surface_CHEMKIN* kineticsMap() const { return kineticsMap_; }

  double MW(unsigned int index) const override;

  void ReactionsAssociatedToSpecies(const unsigned int index, std::vector<unsigned int>& indices) const;
  void IsReactantProduct(const unsigned int reaction_index, double& netStoichiometry) const;

  const Eigen::SparseMatrix<double>& StoichiometricMatrixReactants() const;
  const Eigen::SparseMatrix<double>& StoichiometricMatrixProducts() const;
  const Eigen::SparseMatrix<double>& ReactionOrdersMatrixReactants() const;
  const Eigen::SparseMatrix<double>& ReactionOrdersMatrixProducts() const;

 private:
  OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN* thermodynamicsMap_;
  OpenSMOKE::KineticsMap_Surface_CHEMKIN* kineticsMap_;
};

#include "KineticMapReader_Surface.hpp"
#endif  // KINETICMAPREADER_SURFACE_H
