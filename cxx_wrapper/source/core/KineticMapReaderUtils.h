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

#ifndef KINETICMAPREADERUTILS_H
#define KINETICMAPREADERUTILS_H

// Shared ptree-scanning helpers used by every KineticMapReader_* 
// to parse the optional post-processing blocks out of a
// kinetics.xml/kinetics.surface.xml/... ptree already in memory (no second
// file parse). Lifted, unchanged in logic, from what used to be 
// ProfilesDatabase's private ReadSpeciesClassesBlock/ReadReactionClassesBlock/
// ReadSootProperties
namespace KineticMapReaderUtils {

// Locates an optional post-processing block in a pre-existing ptree.
// The block location is phase-specific: in gas-phase
// <SpeciesClasses>/<ReactionClasses> are under <opensmoke><Kinetics>;
// surface mechanism blocks are in <opensmoke><Kinetics><MaterialKinetics>.
inline boost::optional<const boost::property_tree::ptree&> FindMechanismBlock(
    const boost::property_tree::ptree& root, const std::string& name) {
  static const char* const parents[] = {"opensmoke.Kinetics.",
                                        "opensmoke.Kinetics.MaterialKinetics."};
  for (const char* parent : parents) {
    boost::optional<const boost::property_tree::ptree&> block =
        root.get_child_optional(std::string(parent) + name);
    if (block) return block;
  }
  return boost::none;
}

// --- optional <SpeciesClasses> block ------------------------------------
// Fills class_names / class_members and the per-species class index
// species_to_class (-1 = unclassified). The block is a flat list of
// <ClassSpecies name="..."> nodes whose text is whitespace-separated species indices.
inline bool ReadSpeciesClassesBlock(const boost::property_tree::ptree& mechanism_ptree,
                                    const unsigned int number_of_species,
                                    std::vector<std::string>& class_names,
                                    std::vector<std::vector<unsigned int>>& class_members,
                                    std::vector<int>& species_to_class) {
  class_names.clear();
  class_members.clear();
  species_to_class.assign(number_of_species, -1);

  boost::optional<const boost::property_tree::ptree&> block =
      FindMechanismBlock(mechanism_ptree, "SpeciesClasses");
  if (!block) return false;

  for (const boost::property_tree::ptree::value_type& child : *block) {
    if (child.first != "ClassSpecies") continue;

    const std::string name = child.second.get<std::string>(
        "<xmlattr>.name", "class_" + std::to_string(class_names.size()));

    std::vector<unsigned int> members;
    std::istringstream body(child.second.data());
    long idx;
    while (body >> idx) {
      if (idx < 0 || static_cast<unsigned int>(idx) >= number_of_species) continue;
      species_to_class[idx] = static_cast<int>(class_names.size());
      members.push_back(static_cast<unsigned int>(idx));
    }

    class_names.push_back(name);
    class_members.push_back(members);
  }

  return !class_names.empty();
}

// --- optional <ReactionClasses> block ----------------------------------
// Same idea for reactions. <MainClass name> -> <SubClass name> -> <ReactionIndices>
// (whitespace-separated). The label vectors are always sized to NR,
// "UNSORTED" where a reaction has no entry.
inline bool ReadReactionClassesBlock(const boost::property_tree::ptree& mechanism_ptree,
                                     const unsigned int nr,
                                     std::vector<std::string>& main_class,
                                     std::vector<std::string>& sub_class) {
  main_class.assign(nr, "UNSORTED");
  sub_class.assign(nr, "UNSORTED");

  boost::optional<const boost::property_tree::ptree&> block =
      FindMechanismBlock(mechanism_ptree, "ReactionClasses");
  if (!block) return false;

  for (const boost::property_tree::ptree::value_type& main_child : *block) {
    if (main_child.first != "MainClass") continue;
    const std::string main_name =
        main_child.second.get<std::string>("<xmlattr>.name", "");

    for (const boost::property_tree::ptree::value_type& sub_child : main_child.second) {
      if (sub_child.first != "SubClass") continue;
      const std::string sub_name =
          sub_child.second.get<std::string>("<xmlattr>.name", "");

      boost::optional<const boost::property_tree::ptree&> indices_node =
          sub_child.second.get_child_optional("ReactionIndices");
      if (!indices_node) continue;

      std::istringstream body(indices_node->data());
      long ridx;
      while (body >> ridx) {
        if (ridx < 0 || static_cast<unsigned int>(ridx) >= nr) continue;
        main_class[ridx] = main_name;
        sub_class[ridx] = sub_name;
      }
    }
  }

  return true;
}

// Reads a phase's reaction-names file (a flat <opensmoke.reaction-names> text
// leaf, one whitespace-separated token per reaction, exactly nr of them) into
// reaction_names. Throws if the file is missing - callers whose phase has no
// such file (solid) don't call this at all and synthesize names instead.
inline void ReadReactionNamesFile(const boost::filesystem::path& path,
                                  const unsigned int nr,
                                  std::vector<std::string>& reaction_names) {
  if (!boost::filesystem::exists(path)) {
    throw std::invalid_argument("Kinetic folder does not contain " + path.filename().string());
  }

  boost::property_tree::ptree ptree;
  boost::property_tree::read_xml(path.string(), ptree);

  std::stringstream stream;
  stream.str(ptree.get<std::string>("opensmoke.reaction-names"));

  reaction_names.clear();
  reaction_names.reserve(nr);
  for (unsigned int j = 0; j < nr; j++) {
    std::string reaction_string;
    stream >> reaction_string;
    reaction_names.push_back(reaction_string);
  }
}

}  // namespace KineticMapReaderUtils

#endif  // KINETICMAPREADERUTILS_H
