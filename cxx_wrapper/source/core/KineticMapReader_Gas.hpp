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

KineticMapReader_Gas::KineticMapReader_Gas() {
  phase_label_ = "gas";
  has_ropa_ = true;
  soot_available_ = false;
  soot_number_of_bins_ = 0;
  thermodynamicsMap_ = nullptr;
  kineticsMap_ = nullptr;
}

KineticMapReader_Gas::~KineticMapReader_Gas() {}

bool KineticMapReader_Gas::Load(const std::string& folder_name) {
  boost::filesystem::path folder(folder_name);
  boost::filesystem::path path_mechanism = folder / "kinetics.xml";

  if (!boost::filesystem::exists(path_mechanism)) {
    throw std::invalid_argument(
        "The folder of the kinetic mechanism does not contains any kinetics.xml!");
  }

  boost::property_tree::ptree ptree;
  boost::property_tree::read_xml(path_mechanism.string(), ptree);

  thermodynamicsMap_ = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree, false);
  kineticsMap_ = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMap_, ptree, false);

  species_names_ = thermodynamicsMap_->NamesOfSpecies();

  KineticMapReaderUtils::ReadSpeciesClassesBlock(
      ptree, thermodynamicsMap_->NumberOfSpecies(), species_class_names_,
      species_class_members_, species_to_class_);
  has_species_classes_ = !species_class_names_.empty();

  has_reaction_classes_ = KineticMapReaderUtils::ReadReactionClassesBlock(
      ptree, kineticsMap_->NumberOfReactions(), reaction_main_class_,
      reaction_sub_class_);

  ReadSootProperties(ptree);
  ReadFallOffAndCabrIndices(ptree);

  KineticMapReaderUtils::ReadReactionNamesFile(
      folder / "reaction_names.xml", kineticsMap_->NumberOfReactions(), reaction_names_);

  BuildNameIndexMaps();

  is_available_ = true;
  return true;
}

double KineticMapReader_Gas::MW(unsigned int index) const { return thermodynamicsMap_->MW(index); }

std::string KineticMapReader_Gas::FormattedReactionNameFromIndex(unsigned int index) const {
  const unsigned int n_base = NumberOfReactions();
  const unsigned int n_falloff = static_cast<unsigned int>(indices_of_falloff_reactions_.size());
  const unsigned int n_cabr = static_cast<unsigned int>(indices_of_cabr_reactions_.size());

  if (index < n_base) {
    return "R" + std::to_string(index + 1) + ": " + ReactionNameFromIndex(index);
  }

  unsigned int global_index;  // 1-based, into the base reaction numbering
  if (index < n_base + n_falloff) {
    global_index = indices_of_falloff_reactions_[index - n_base];
  } else if (index < n_base + n_falloff + n_cabr) {
    global_index = indices_of_cabr_reactions_[index - n_base - n_falloff];
  } else {
    throw std::invalid_argument("Reaction index out of range: " + std::to_string(index));
  }
  return "R" + std::to_string(global_index) + "(inf): " + ReactionNameFromIndex(global_index - 1);
}

  // BuildStoichiometricMatrix()/BuildReactionOrdersMatrix() are called only once,
  // the first instance creates them and they exist aftwerwards.
void KineticMapReader_Gas::ReactionsAssociatedToSpecies(const unsigned int index,
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

void KineticMapReader_Gas::IsReactantProduct(const unsigned int reaction_index,
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

const Eigen::SparseMatrix<double>& KineticMapReader_Gas::StoichiometricMatrixReactants() const {
  kineticsMap_->stoichiometry().BuildStoichiometricMatrix();
  return kineticsMap_->stoichiometry().stoichiometric_matrix_reactants();
}

const Eigen::SparseMatrix<double>& KineticMapReader_Gas::StoichiometricMatrixProducts() const {
  kineticsMap_->stoichiometry().BuildStoichiometricMatrix();
  return kineticsMap_->stoichiometry().stoichiometric_matrix_products();
}

// --- optional <SootProperties> block ---------------------------------
// PolimiSoot BIN properties, one value per bin per field.
void KineticMapReader_Gas::ReadSootProperties(const boost::property_tree::ptree& mechanism_ptree) {
  soot_available_ = false;
  soot_number_of_bins_ = 0;

  boost::optional<const boost::property_tree::ptree&> block =
      mechanism_ptree.get_child_optional("opensmoke.SootProperties");
  if (!block) return;

  {
    std::stringstream s(block->get<std::string>("NumberOfBins"));
    s >> soot_number_of_bins_;
  }

  const auto read_uint = [&](const char* key, std::vector<unsigned int>& v) {
    std::stringstream s(block->get<std::string>(key));
    v.assign(soot_number_of_bins_, 0u);
    for (unsigned int k = 0; k < soot_number_of_bins_; k++) s >> v[k];
  };
  const auto read_int = [&](const char* key, std::vector<int>& v) {
    std::stringstream s(block->get<std::string>(key));
    v.assign(soot_number_of_bins_, 0);
    for (unsigned int k = 0; k < soot_number_of_bins_; k++) s >> v[k];
  };
  const auto read_double = [&](const char* key, std::vector<double>& v) {
    std::stringstream s(block->get<std::string>(key));
    v.assign(soot_number_of_bins_, 0.);
    for (unsigned int k = 0; k < soot_number_of_bins_; k++) s >> v[k];
  };

  read_uint("Bin_index", soot_bin_index_);
  read_int("Bin_section", soot_bin_section_);
  read_double("Bin_nc", soot_bin_nc_);
  read_double("Bin_nh", soot_bin_nh_);
  read_double("Bin_no", soot_bin_no_);
  read_double("Bin_htoc", soot_bin_htoc_);
  read_double("Bin_mw", soot_bin_mw_);
  read_double("Bin_density", soot_bin_density_);
  read_double("Bin_volume", soot_bin_volume_);
  read_double("Bin_mass", soot_bin_mass_);
  read_double("Bin_numpp", soot_bin_numpp_);
  read_double("Bin_dsph", soot_bin_dsph_);
  read_double("Bin_dcol", soot_bin_dcol_);
  read_double("Bin_dpp", soot_bin_dpp_);
  read_double("Bin_df", soot_bin_df_);

  soot_available_ = true;
}

void KineticMapReader_Gas::ReadFallOffAndCabrIndices(const boost::property_tree::ptree& mechanism_ptree) {
  const auto read_indices = [&](const char* key, std::vector<unsigned int>& v) {
    v.clear();
    boost::optional<std::string> text = mechanism_ptree.get_optional<std::string>(
        std::string("opensmoke.Kinetics.") + key);
    if (!text) return;
    std::stringstream s(*text);
    unsigned int n = 0;
    s >> n;
    v.assign(n, 0u);
    for (unsigned int k = 0; k < n; k++) s >> v[k];
  };

  read_indices("FallOff", indices_of_falloff_reactions_);
  read_indices("CABR", indices_of_cabr_reactions_);
}
