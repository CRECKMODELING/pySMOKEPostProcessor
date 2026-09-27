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

#ifndef SPECIESCLASSES_H
#define SPECIESCLASSES_H

#include "core/PostProcessorCore.h"

// Post-processing of the <SpeciesClasses> block.
class SpeciesClass {
 public:
  SpeciesClass();
  void SetResults(PostProcessorCore* data);

  // Checks if the SpeciesClasses block is actually available
  bool speciesClassesAvailable() const { return species_classes_available_; }

  // --- class-wise elemental distribution along the whole abscissa ------------
  // For every class k and abscissa point i, accumulates the moles of the given
  // element carried by that class (per unit mass of mixture). normalize = true
  // rescales every abscissa column so the classes sum to one (it always should for gas-only).
  void ElementalDistribution(const std::string element, const bool normalize = false);

  // --- per-species elemental distribution ---
  // Same numeric core as ElementalDistribution, one row per species instead of
  // per class - does NOT require a <SpeciesClasses> block (every species in the
  // mechanism is covered, classified or not). ElementalDistribution calls this
  // internally and aggregates by class afterwards.
  // The location of this function is undoubtedly wrong, it should be on his own.
  void ElementMolesBySpecies(const std::string element);

  inline const std::vector<std::string>& speciesNames() const { return species_names_; }
  inline const std::vector<std::vector<double>>& elementMolesBySpecies() const {
    return element_moles_by_species_;
  }

  // --- element flux between species classes --------------------------------
  // FluxByClass() first builds the full directed class -> class element-flux
  // matrix (every reaction contributes, species filtering is only applied later)
  // Then, the full matrix class -> class is filtered at the end.
  void SetFluxPerClass(const bool flux_per_class);    // Choose target, either a speciesClass or a species
  void SetClassName(const std::string class_name);    // set target when flux_per_class == true
  void SetSpecies(const std::string species);         // set target when flux_per_class == false
  void SetElement(const std::string element);
  void SetFluxAnalysisType(const std::string type);  // production | destruction
  void SetLocalValue(const double localValue);
  void SetThickness(const std::string thickness);    // absolute | relative
  void SetLabelType(const std::string type);         // absolute | relative
  void SetThicknessLogScale(const bool thicknesslogscale);
  void SetWidth(const int width);                     // max links per class node (graph)
  void SetDepth(const int depth);                     // class-graph generations from the target
  void SetThreshold(const double threshold);          // min % of a node's flux to keep a link
  void SetCarbonWeighted(const bool carbon_weighted); // see FluxByClass() below (default true)
  void SetAutoPruneDiagonal(const bool auto_prune_diagonal);  // zero fluxMatrix()'s diagonal (default true)
  void FluxByClass();

  inline const std::vector<std::string>& classNames() const { return class_names_; }
  inline const std::vector<int>& speciesToClass() const { return species_to_class_; }

  inline const std::vector<double>& abscissa() const { return abscissa_; }
  inline const std::vector<std::vector<double>>& elementalFractions() const {
    return elemental_fractions_;
  }

  // full directed matrix: fluxMatrix()[i][j] = element flux from class i to class j.
  // The diagonal (i == j, intra-class flux) is zeroed when auto_prune_diagonal_ is
  // true (default): those entries are typically far larger than any cross-class
  // one and make a heatmap of the matrix unreadable; the graph never draws them
  // regardless (self-loops are excluded from the walk).
  inline const std::vector<std::vector<double>>& fluxMatrix() const { return flux_matrix_; }
  inline const std::vector<int>& indexFirstClass() const { return index_first_class_; }   // "Reactant" speciesClass
  inline const std::vector<int>& indexSecondClass() const { return index_second_class_; } // "Product"  speciesClass
  inline const std::vector<double>& computedThickness() const { return computed_thickness_; } // Thickness
  inline const std::vector<double>& computedLabel() const { return computed_label_; }         // Label

 protected:
  void ReadSpeciesClasses();
  unsigned int ElementIndex(const std::string& element) const;

  PostProcessorCore* data_;
  bool species_classes_available_; // Boolean flag

  std::vector<std::string> class_names_;                  // Class names (e.g. SOOT-X, AROMATICS, C2, etc)
  std::vector<std::vector<unsigned int>> class_species_;  // External vector: classes. Internal vector: indices of species inside that class.
  std::vector<int> species_to_class_;                     // For each species, index of the class it belongs to (-1 if none)

  std::string species_;
  std::string class_name_;
  std::string element_;
  std::string flux_type_;
  std::string thickness_;
  std::string label_type_;
  bool thicknesslogscale_;
  bool carbon_weighted_;
  bool flux_per_class_;
  bool auto_prune_diagonal_;
  int width_;
  int depth_;
  double threshold_;
  double localValue_;

  std::vector<double> abscissa_;                          // [point]
  std::vector<std::vector<double>> elemental_fractions_;  // [class][point]

  std::vector<std::string> species_names_;                     // [species], every species in the mechanism
  std::vector<std::vector<double>> element_moles_by_species_;  // [species][point]

  std::vector<std::vector<double>> flux_matrix_;       // [source class][target class], full directed
  std::vector<int> index_first_class_;                 // [edge] -> source class
  std::vector<int> index_second_class_;                // [edge] -> target class
  std::vector<double> computed_thickness_;             // [edge]
  std::vector<double> computed_label_;                 // [edge] -> flux (% of node, or absolute)
};

#include "SpeciesClasses.hpp"
#endif  // SPECIESCLASSES_H
