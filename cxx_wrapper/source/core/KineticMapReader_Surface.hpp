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

KineticMapReader_Surface::KineticMapReader_Surface() {
  phase_label_ = "surface";
  has_ropa_ = true;
  thermodynamicsMap_ = nullptr;
  kineticsMap_ = nullptr;
}

KineticMapReader_Surface::~KineticMapReader_Surface() {}

bool KineticMapReader_Surface::Load(const std::string& folder_name) {
  boost::filesystem::path folder(folder_name);
  boost::filesystem::path path_mechanism = folder / "kinetics.surface.xml";

  if (!boost::filesystem::exists(path_mechanism)) {
    throw std::invalid_argument(
        "The folder of the kinetic mechanism does not contains the heterogeneous kinetics.xml!");
  }

  boost::property_tree::ptree ptree;
  boost::property_tree::read_xml(path_mechanism.string(), ptree);

  // When constructing the thermo/kinetics surface maps, OpenSMOKEpp prints a
  // lot of (useless) stuff to stdout; send it to null for the duration of the call.
  std::streambuf* old_buf = std::cout.rdbuf();
  std::ofstream null_stream("/dev/null");
  std::cout.rdbuf(null_stream.rdbuf());

  thermodynamicsMap_ = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN(ptree);
  kineticsMap_ = new OpenSMOKE::KineticsMap_Surface_CHEMKIN(*thermodynamicsMap_, ptree);

  std::cout.rdbuf(old_buf);

  species_names_ = thermodynamicsMap_->NamesOfSpecies();

  has_reaction_classes_ = KineticMapReaderUtils::ReadReactionClassesBlock(
      ptree, kineticsMap_->NumberOfReactions(), reaction_main_class_, reaction_sub_class_);

  KineticMapReaderUtils::ReadReactionNamesFile(folder / "surface_reaction_names.xml",
                                               kineticsMap_->NumberOfReactions(),
                                               reaction_names_);

  BuildNameIndexMaps();

  is_available_ = true;
  return true;
}

double KineticMapReader_Surface::MW(unsigned int index) const { return thermodynamicsMap_->MW(index); }

void KineticMapReader_Surface::ReactionsAssociatedToSpecies(const unsigned int index,
                                                             std::vector<unsigned int>& indices) const {
  kineticsMap_->stoichiometry().BuildStoichiometricMatrix();

  for (int k = 0; k < kineticsMap_->stoichiometry().stoichiometric_matrix_reactants().outerSize();
       ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(
             kineticsMap_->stoichiometry().stoichiometric_matrix_reactants(), k);
         it; ++it)
      if (it.col() == index) indices.push_back(it.row());
  }
  for (int k = 0; k < kineticsMap_->stoichiometry().stoichiometric_matrix_products().outerSize();
       ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(
             kineticsMap_->stoichiometry().stoichiometric_matrix_products(), k);
         it; ++it)
      if (it.col() == index) indices.push_back(it.row());
  }

  std::sort(indices.begin(), indices.end());
}

void KineticMapReader_Surface::IsReactantProduct(const unsigned int reaction_index,
                                                 double& netStoichiometry) const {
  kineticsMap_->stoichiometry().BuildStoichiometricMatrix();

  std::vector<double> reactants_stoich;
  std::vector<double> products_stoich;
  std::vector<double> reactants_indices;
  std::vector<double> products_indices;
  std::vector<double> duplicate_species_indices;

  Eigen::SparseMatrix<double> reactants = kineticsMap_->stoichiometry().stoichiometric_matrix_reactants();
  Eigen::SparseMatrix<double> products = kineticsMap_->stoichiometry().stoichiometric_matrix_products();

  for (int k = 0; k < reactants.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(reactants, k); it; ++it) {
      if (it.row() == reaction_index) {
        reactants_stoich.push_back(it.value());
        reactants_indices.push_back(it.col());
      }
    }
  }

  for (int k = 0; k < products.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(products, k); it; ++it) {
      if (it.row() == reaction_index) {
        products_stoich.push_back(it.value());
        products_indices.push_back(it.col());
      }
    }
  }

  std::sort(reactants_indices.begin(), reactants_indices.end());
  std::sort(products_indices.begin(), products_indices.end());

  std::vector<double> common_species(reactants_indices.size() + products_indices.size());
  std::vector<double>::iterator it, end;

  end = std::set_intersection(reactants_indices.begin(), reactants_indices.end(),
                              products_indices.begin(), products_indices.end(),
                              common_species.begin());

  for (it = common_species.begin(); it != end; it++) duplicate_species_indices.push_back(*it);

  if (reactants_indices.size() != 1) {
    netStoichiometry = 1;
  } else if (duplicate_species_indices.size() == 1) {
    double idx = duplicate_species_indices[0];
    int pos_r = 0;
    int pos_p = 0;

    std::vector<double>::iterator it_r = std::find(reactants_indices.begin(), reactants_indices.end(), idx);
    if (it_r != reactants_indices.end()) pos_r = it_r - reactants_indices.begin();

    std::vector<double>::iterator it_p = std::find(products_indices.begin(), products_indices.end(), idx);
    if (it_r != products_indices.end()) pos_p = it_p - products_indices.begin();

    netStoichiometry = -reactants_stoich[pos_r] + products_stoich[pos_p];
  } else if (duplicate_species_indices.size() == 0) {
    netStoichiometry = 1;
  } else {
    std::string msg = "Something is wrong with the reaction you are asking for!";
    msg += "Reaction id: " + std::to_string(reaction_index);
    throw std::invalid_argument(msg);
  }
}

const Eigen::SparseMatrix<double>& KineticMapReader_Surface::StoichiometricMatrixReactants() const {
  kineticsMap_->stoichiometry().BuildStoichiometricMatrix();
  return kineticsMap_->stoichiometry().stoichiometric_matrix_reactants();
}

const Eigen::SparseMatrix<double>& KineticMapReader_Surface::StoichiometricMatrixProducts() const {
  kineticsMap_->stoichiometry().BuildStoichiometricMatrix();
  return kineticsMap_->stoichiometry().stoichiometric_matrix_products();
}
