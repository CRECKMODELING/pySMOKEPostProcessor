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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

ReactionClass::ReactionClass() {
  data_ = nullptr;
  heterogeneous_ = false;
  classes_available_ = false;
}

void ReactionClass::SetHeterogeneous(const bool heterogeneous) {
  heterogeneous_ = heterogeneous;
}

void ReactionClass::SetDatabase(ProfilesDatabase* data) {
  data_ = data;
  ReadReactionClasses();
}

// The <ReactionClasses> block is parsed in 
// ProfilesDatabase::ReadKineticMechanism / ReadHeterogeneousKineticMechanism
void ReactionClass::ReadReactionClasses() {
  classes_available_ = false;
  main_class_.clear();
  sub_class_.clear();
  merge_group_.clear();

  if (data_ == nullptr) return;

  if (heterogeneous_ == true) {
    main_class_ = data_->reaction_main_class_heterogeneous_;
    sub_class_ = data_->reaction_sub_class_heterogeneous_;
    classes_available_ = data_->has_reaction_classes_heterogeneous_;
  } else {
    main_class_ = data_->reaction_main_class_;
    sub_class_ = data_->reaction_sub_class_;
    classes_available_ = data_->has_reaction_classes_;
  }

  ComputeMergeGroups();
}

// Serializes one reaction's reactant or product side, sorted by species index, from
// the mechanism stoichiometric matrix (species index + coefficient per
// reaction)
std::string ReactionClass::SerializeReactionSide(
    const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix, const int row) {
  std::vector<std::pair<int, double>> entries;
  for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(matrix, row); it;
       ++it)
    entries.emplace_back(static_cast<int>(it.col()), it.value());
  std::sort(entries.begin(), entries.end());

  std::ostringstream oss;
  for (const std::pair<int, double>& entry : entries)
    oss << entry.first << ':' << entry.second << ';';
  return oss.str();
}

// This only groups reactions by the mechanism stoichiometric matrix (static, computed once here); 
// it does not decide which member represents the group in a report - that choice depends 
// on the ROPA result being merged and is made per call, see MergeDuplicates().
void ReactionClass::ComputeMergeGroups() {
  const unsigned int nr = static_cast<unsigned int>(main_class_.size());

  // Eigen's stoichiometric_matrix_reactants()/products() are column-major; converted
  // once to row-major copies so InnerIterator(matrix, reaction_row) below actually
  // walks that reaction's species, not a column's.
  Eigen::SparseMatrix<double, Eigen::RowMajor> reactants;
  Eigen::SparseMatrix<double, Eigen::RowMajor> products;
  if (heterogeneous_ == true) {
    reactants =
        data_->kineticsMapSurfaceXML->stoichiometry().stoichiometric_matrix_reactants();
    products =
        data_->kineticsMapSurfaceXML->stoichiometry().stoichiometric_matrix_products();
  } else {
    reactants = data_->kineticsMapXML->stoichiometry().stoichiometric_matrix_reactants();
    products = data_->kineticsMapXML->stoichiometry().stoichiometric_matrix_products();
  }

  merge_group_.assign(nr, -1);
  std::unordered_map<std::string, int> canonical_to_group;
  for (unsigned int i = 0; i < nr; i++) {
    const std::string reactant_side =
        SerializeReactionSide(reactants, static_cast<int>(i));
    const std::string product_side = SerializeReactionSide(products, static_cast<int>(i));
    const std::string forward = reactant_side + "=" + product_side;
    const std::string backward = product_side + "=" + reactant_side;
    const std::string canonical = std::min(forward, backward);

    const auto inserted = canonical_to_group.emplace(
        canonical, static_cast<int>(canonical_to_group.size()));
    merge_group_[i] = inserted.first->second;
  }
}

void ReactionClass::MergeDuplicates(
    const std::vector<std::vector<int>>& species_indices,
    const std::vector<std::vector<double>>& species_coefficients,
    std::vector<int>& representative_indices,
    std::vector<std::vector<double>>& merged_coefficients) const {
  representative_indices.clear();
  const size_t n_species = species_indices.size();
  merged_coefficients.assign(n_species, {});

  // Every reaction index touched by any species this call, with its per-species
  // coefficient (0 where a given species doesn't cite it).
  std::unordered_map<int, std::vector<double>> idx_values;
  for (size_t s = 0; s < n_species; s++) {
    const std::vector<int>& idxs = species_indices[s];
    const std::vector<double>& coeffs = species_coefficients[s];
    for (size_t k = 0; k < idxs.size(); k++) {
      const int idx = idxs[k];
      if (idx < 0 || static_cast<size_t>(idx) >= merge_group_.size()) continue;
      auto it = idx_values.find(idx);
      if (it == idx_values.end())
        it = idx_values.emplace(idx, std::vector<double>(n_species, 0.0)).first;
      it->second[s] += coeffs[k];
    }
  }

  // Bucket the touched reactions by their (static) merge group.
  std::unordered_map<int, std::vector<int>> group_members;
  for (const std::pair<const int, std::vector<double>>& entry : idx_values)
    group_members[merge_group_[entry.first]].push_back(entry.first);

  for (const std::pair<const int, std::vector<int>>& group_entry : group_members) {
    const std::vector<int>& members = group_entry.second;

    // This query's representative: restrict to classified (non-UNSORTED) members
    // if any exist in the group, then take whichever has the largest |coefficient|
    // in any requested species column - i.e. whichever direction actually
    // dominates under the current local conditions. Ties go to the smallest
    // reaction index; order-independent, so unordered_map iteration order above
    // doesn't matter.
    bool any_classified = false;
    for (int idx : members)
      if (main_class_[idx] != "UNSORTED") {
        any_classified = true;
        break;
      }

    int best = -1;
    double best_score = -1.0;
    for (int idx : members) {
      if (any_classified && main_class_[idx] == "UNSORTED") continue;
      double score = 0.0;
      for (double v : idx_values.at(idx)) score = std::max(score, std::fabs(v));
      if (score > best_score || (score == best_score && idx < best)) {
        best_score = score;
        best = idx;
      }
    }

    const size_t slot = representative_indices.size();
    representative_indices.push_back(best);
    for (size_t s = 0; s < n_species; s++) merged_coefficients[s].push_back(0.0);

    for (int idx : members) {
      const std::vector<double>& values = idx_values.at(idx);
      for (size_t s = 0; s < n_species; s++) merged_coefficients[s][slot] += values[s];
    }
  }
}
