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

#ifndef REACTIONCLASSES_H
#define REACTIONCLASSES_H

#include "ProfilesDatabase.h"

// Post-processing of the mechanism's <ReactionClasses> block.
//
// The ReactionClasses block and the mechanism stoichiometric map are read once.
// Two things that reaction-class flux post-processing needs:
//  - which reactions are the *same* physical reaction written twice (an explicit
//    reversible forward/backward pair, or an outright duplicate declaration) -
//    decided from the mechanism's stoichiometric matrix (species index +
//    coefficient per reaction side)
//  - merging one or more species' ROPA (reaction_index, coefficient) pairs onto
//    one representative index per such group, chosen per call based on which
//    direction actually dominates the current query (see MergeDuplicates()).
//
// Grouping by class/subclass/reaction-type, thresholding, plotting remained in Python:
// it is cheap once duplicates are merged, and its shape varies with
// user-chosen sort criteria.
class ReactionClass {
 public:
  ReactionClass();

  void SetDatabase(ProfilesDatabase* data);
  void SetHeterogeneous( const bool heterogeneous);  
  // false (default) = gas, true = surface

  // Checks the <ReactionClasses> block was actually found.
  inline bool classesAvailable() const { return classes_available_; }

  // Per-reaction (0-based, sized to the selected phase's NumberOfReactions) raw
  // <ReactionClasses> labels; "UNSORTED" where a reaction has no entry - replaces
  // the python-side KineticMap.Classes()/rxnclass/rxnsubclass for this feature.
  inline const std::vector<std::string>& mainClass() const { return main_class_; }
  inline const std::vector<std::string>& subClass() const { return sub_class_; }

  // Merges one or more species' ROPA (reaction_index, coefficient) pairs
  // dealing with fw/bw reactions
  //
  // The representative is chosen per call, not fixed by the mechanism alone
  // (e.g. recombinations and decomposition being opposite of each other)
  // which direction should represent the merged flux is then a property of 
  // current local conditions, not of the mechanism only. Within each
  // group, the choice is: restrict to classified (non-UNSORTED) members if any
  // exist, then take whichever member's |coefficient| is largest in any of the
  // provided species columns (ties go to the smallest reaction index - very unlikely).
  //
  // representative_indices gets one entry per group touched by any species.
  // merged_coefficients gets one entry per species, each the same length as
  // representative_indices and in the same order - ready to build one pandas
  // column per species, all sharing the same row index.
  void MergeDuplicates(const std::vector<std::vector<int>>& species_indices,
                       const std::vector<std::vector<double>>& species_coefficients,
                       std::vector<int>& representative_indices,
                       std::vector<std::vector<double>>& merged_coefficients) const;

 protected:
  void ReadReactionClasses();
  void ComputeMergeGroups(); // [LG] Maybe we should extend this to all ROPAs/Sensitivities
  static std::string SerializeReactionSide(
      const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix, const int row);
  // The Eigen matrix is RowMajor to be able to use InnerIterator which is faster
  // when getting all reactions of the i-th species

  ProfilesDatabase* data_;
  bool heterogeneous_;
  bool classes_available_;

  std::vector<std::string> main_class_;  // per reaction, <MainClass name=...>
  std::vector<std::string> sub_class_;   // per reaction, <SubClass name=...>
  std::vector<int> merge_group_;  // per reaction, duplicate/reversible-pair group id
};

#include "ReactionClasses.hpp"
#endif  // REACTIONCLASSES_H
