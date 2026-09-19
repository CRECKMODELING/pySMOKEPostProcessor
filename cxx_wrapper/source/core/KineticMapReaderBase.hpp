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

KineticMapReaderBase::KineticMapReaderBase() {
  is_available_ = false;
  has_ropa_ = false;
  has_species_classes_ = false;
  has_reaction_classes_ = false;
}

KineticMapReaderBase::~KineticMapReaderBase() {}

unsigned int KineticMapReaderBase::NumberOfReactions() const {
  return static_cast<unsigned int>(reaction_names_.size());
}

unsigned int KineticMapReaderBase::NumberOfSpecies() const {
  return static_cast<unsigned int>(species_names_.size());
}

std::string KineticMapReaderBase::ReactionNameFromIndex(unsigned int index) const {
  if (index >= reaction_names_.size()) {
    throw std::invalid_argument("Reaction index out of range: " + std::to_string(index));
  }
  return reaction_names_[index];
}

std::string KineticMapReaderBase::SpeciesNameFromIndex(unsigned int index) const {
  if (index >= species_names_.size()) {
    throw std::invalid_argument("Species index out of range: " + std::to_string(index));
  }
  return species_names_[index];
}

int KineticMapReaderBase::ReactionIndexFromName(const std::string& name) const {
  std::unordered_map<std::string, int>::const_iterator it = reaction_index_by_name_.find(name);
  if (it == reaction_index_by_name_.end()) {
    throw std::invalid_argument("Unknown reaction: " + name);
  }
  return it->second;
}

std::string KineticMapReaderBase::FormattedReactionNameFromIndex(unsigned int index) const {
  return "R" + std::to_string(index + 1) + ": " + ReactionNameFromIndex(index);
}

int KineticMapReaderBase::SpeciesIndexFromName(const std::string& name) const {
  std::unordered_map<std::string, int>::const_iterator it = species_index_by_name_.find(name);
  if (it == species_index_by_name_.end()) {
    throw std::invalid_argument("Unknown species: " + name);
  }
  return it->second;
}

int KineticMapReaderBase::SpeciesClassOf(const std::string& species_name) const {
  std::unordered_map<std::string, int>::const_iterator it = species_index_by_name_.find(species_name);
  if (it == species_index_by_name_.end()) {
    return -1;
  }
  const unsigned int idx = static_cast<unsigned int>(it->second);
  if (idx >= species_to_class_.size()) {
    return -1;
  }
  return species_to_class_[idx];
}

void KineticMapReaderBase::BuildNameIndexMaps() {
  // emplace (not operator[]): on a duplicate name, keep the first occurrence,
  // matching the previous Python implementation's list.index() semantics.
  // Mechanisms can legitimately declare the same reaction name twice (explicit
  // duplicates); silently letting the last one win would change which physical
  // reaction every by-name lookup resolves to.
  reaction_index_by_name_.clear();
  for (unsigned int j = 0; j < reaction_names_.size(); j++) {
    reaction_index_by_name_.emplace(reaction_names_[j], static_cast<int>(j));
  }

  species_index_by_name_.clear();
  for (unsigned int j = 0; j < species_names_.size(); j++) {
    species_index_by_name_.emplace(species_names_[j], static_cast<int>(j));
  }
}
