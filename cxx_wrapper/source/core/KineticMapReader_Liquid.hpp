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

#include "KineticMapReaderUtils.h"

KineticMapReader_Liquid::KineticMapReader_Liquid() {
  phase_label_ = "liquid";
  has_ropa_ = true;
  thermodynamicsMap_ = nullptr;
  kineticsMap_ = nullptr;
}

KineticMapReader_Liquid::~KineticMapReader_Liquid() {}

bool KineticMapReader_Liquid::Load(const std::string& folder_name) {
  boost::filesystem::path folder(folder_name);
  boost::filesystem::path path_mechanism = folder / "kinetics.liquid.xml";

  if (!boost::filesystem::exists(path_mechanism)) {
    throw std::invalid_argument(
        "The folder of the kinetic mechanism does not contains the liquid kinetics.liquid.xml!");
  }

  boost::property_tree::ptree ptree;
  boost::property_tree::read_xml(path_mechanism.string(), ptree);

  // Same stdout-suppression precaution taken for the surface map constructor -
  // OpenSMOKEpp's non-gas map constructors tend to print a lot of setup noise.
  std::streambuf* old_buf = std::cout.rdbuf();
  std::ofstream null_stream("/dev/null");
  std::cout.rdbuf(null_stream.rdbuf());

  thermodynamicsMap_ = new OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN(ptree);
  // Material 1 (1-based): the single-non-gas-phase assumption this whole
  // design makes means there is always exactly one liquid material to pick.
  kineticsMap_ = new OpenSMOKE::KineticsMap_Liquid_CHEMKIN(*thermodynamicsMap_, ptree, 1u);

  std::cout.rdbuf(old_buf);

  species_names_ = thermodynamicsMap_->NamesOfSpecies();

  KineticMapReaderUtils::ReadReactionNamesFile(folder / "reaction_names.liquid.xml",
                                               kineticsMap_->NumberOfReactions(), reaction_names_);

  BuildNameIndexMaps();

  is_available_ = true;
  return true;
}

double KineticMapReader_Liquid::MW(unsigned int index) const { return thermodynamicsMap_->MW(index); }
