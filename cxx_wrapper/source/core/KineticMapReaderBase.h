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

#ifndef KINETICMAPREADERBASE_H
#define KINETICMAPREADERBASE_H

#include <unordered_map>

// Common read-only for every phase's kinetics reader exposes: 
// * reaction/species name<->index lookups
// * per-species MW, 
// * optional <SpeciesClasses>/<ReactionClasses> blocks. 
// reaction_names_, species_names_, the class blocks are filled by each 
// KineticMapReader_* subclass's own Load(); 
// this base is only a common interface but not a parser
// (that's KineticMapReaderUtils)
// MW(index) is the one thing every phase computes differently
// so it stays pure virtual.
class KineticMapReaderBase {
 public:
  KineticMapReaderBase();
  virtual ~KineticMapReaderBase();

  unsigned int NumberOfReactions() const;
  unsigned int NumberOfSpecies() const;

  std::string ReactionNameFromIndex(unsigned int index) const;  // 0-based
  std::string SpeciesNameFromIndex(unsigned int index) const;   // 0-based

  // "R{n}: {name}", 1-based n (same as previous implementations),
  // formerly reconstructed in Python by the KineticMap class. 
  // This default is enough for every phase whose widgets only ever report an
  // index in [0, NumberOfReactions).
  // The Gas phase requires more due to the presence of FallOff and CAB Reactions
  // and has therefore his own override.
  // Phase-specific reactions, like stick, coverage, and similar may require this
  // to change. The additional index is only required for sensitivity analyses.
  virtual std::string FormattedReactionNameFromIndex(unsigned int index) const;
  int ReactionIndexFromName(const std::string& name) const;     // throws if unknown
  int SpeciesIndexFromName(const std::string& name) const;      // throws if unknown

  const std::vector<std::string>& ReactionNames() const { return reaction_names_; }
  const std::vector<std::string>& SpeciesNames() const { return species_names_; }

  bool HasSpeciesClasses() const { return has_species_classes_; }
  const std::vector<std::string>& SpeciesClassNames() const { return species_class_names_; }
  const std::vector<std::vector<unsigned int>>& SpeciesClassMembers() const {
    return species_class_members_;
  }
  const std::vector<int>& SpeciesToClass() const { return species_to_class_; }
  int SpeciesClassOf(const std::string& species_name) const;  // -1 if unclassified/unknown

  bool HasReactionClasses() const { return has_reaction_classes_; }
  const std::vector<std::string>& ReactionMainClass() const { return reaction_main_class_; }
  const std::vector<std::string>& ReactionSubClass() const { return reaction_sub_class_; }

  const std::string& PhaseLabel() const { return phase_label_; }
  bool IsAvailable() const { return is_available_; }
  bool HasRopa() const { return has_ropa_; }

  virtual double MW(unsigned int index) const = 0;

 protected:
  std::string phase_label_;
  bool is_available_;
  bool has_ropa_;

  std::vector<std::string> reaction_names_;
  std::vector<std::string> species_names_;
  std::unordered_map<std::string, int> reaction_index_by_name_;
  std::unordered_map<std::string, int> species_index_by_name_;

  bool has_species_classes_;
  std::vector<std::string> species_class_names_;
  std::vector<std::vector<unsigned int>> species_class_members_;
  std::vector<int> species_to_class_;

  bool has_reaction_classes_;
  std::vector<std::string> reaction_main_class_;
  std::vector<std::string> reaction_sub_class_;

  // Rebuilds reaction_index_by_name_/species_index_by_name_ from
  // reaction_names_/species_names_. Every subclass's Load() calls this once
  // both vectors are final.
  void BuildNameIndexMaps();
};

#include "KineticMapReaderBase.hpp"
#endif  // KINETICMAPREADERBASE_H
