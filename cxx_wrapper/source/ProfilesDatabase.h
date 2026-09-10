/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Authors: Timoteo Dinelli <timoteo.dinelli@polimi.it>                  |
|            Edoardo Ramalli <edoardo.ramalli@polimi.it>                  |
|   CRECK Modeling Group <www.creckmodeling.polimi.it>                    |
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

#ifndef PROFILESDATABASE_H
#define PROFILESDATABASE_H

#include <Eigen/Sparse>

class ProfilesDatabase {
 public:
  ProfilesDatabase(void);
  ~ProfilesDatabase(void);

  bool ReadKineticMechanism(const std::string& folder_name);
  bool ReadHeterogeneousKineticMechanism(const std::string& folder_name, const std::string& phase_name); // Function added

  bool ReadFileResults(const std::string& folder_name, bool isHeterogeneous); // Modified number of input (added the isHeterogeneous variable)

  void Prepare();
  void PrepareHeterogeneous();  // Function added

  void SpeciesCoarsening(const double threshold);

  int number_of_abscissas_;
  int number_of_ordinates_;

  std::vector<int> column_index_of_massfractions_profiles;
  std::vector<std::string> string_list_additional;
  std::vector<int> list_of_conversion_species_;

  std::vector<std::string> string_list_massfractions_sorted;
  std::vector<int> sorted_index;
  std::vector<int> current_sorted_index;
  std::vector<double> sorted_max;

  std::vector<std::vector<double>> omega;
  std::vector<std::vector<double>> additional;

  unsigned int index_T;
  unsigned int index_P;
  unsigned int index_MW;
  unsigned int index_density;
  unsigned int index_velocity;
  unsigned int index_mass_flow_rate;
  unsigned int index_volume;

  unsigned int index_x_coord;
  unsigned int index_z_coord;

  // Surface-specific properties and quantities
  std::vector<std::vector<double>> Z;
  std::vector<std::vector<double>> massBulk;
  unsigned int index_area_over_volume;
  unsigned int index_surface_sites_concentration;
  std::vector<int> column_index_of_surfacefractions_profiles;
  std::vector<int> column_index_of_bulkmasses_profiles;
  std::vector<std::string> string_list_surfacefractions_sorted;
  std::vector<std::string> string_list_bulkmasses_sorted;
  std::vector<int> sorted_index_surface;
  std::vector<int> current_sorted_index_surface;
  std::vector<double> sorted_max_surface;
  std::vector<int> sorted_index_bulk;
  std::vector<int> current_sorted_index_bulk;
  std::vector<double> sorted_max_bulk;
  unsigned int number_of_gas_species;
  unsigned int number_of_surface_species;
  unsigned int number_of_bulk_species;

  std::vector<double> mw_species_;  // Note: for now, saving only the gas-phase molecular weights. 
      // Since surface species are already saved in molar fractions, and bulk species are in mass, 
      // but their activity is 1, no need to save their MW (no use for it)
      // Maybe, in different phases it is required to save those.

  boost::property_tree::ptree xml_main_input;

  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML;
  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML;
  OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN* thermodynamicsMapSurfaceXML;
  OpenSMOKE::KineticsMap_Surface_CHEMKIN* kineticsMapSurfaceXML;

  bool iSensitivityEnabled_;
  bool iROPAEnabled_;
  bool is_kinetics_available_;

  bool is_mechanism_heterogeneous_;
  bool iROPAHeterogeneousEnabled_;
  bool iSensitivityHeterogeneousEnabled_;
  bool is_heterogeneous_kinetics_available_;
  
  // Folder of kinetic mechanism and output are the same for heterogeneous mechanisms, no need to duplicate.
  boost::filesystem::path path_folder_results_;
  boost::filesystem::path path_folder_mechanism_;

  void ReactionsAssociatedToSpecies(const unsigned int index, std::vector<unsigned int>& indices);
  void isReactantProduct(const unsigned int reaction_index, double& netStoichiometry);
  
  void ReactionsAssociatedToSpecies_Surface(const unsigned int index, std::vector<unsigned int>& indices);  // Function added
  void isReactantProduct_Surface(const unsigned int reaction_index, double& netStoichiometry);              // Function added
  // In here I think its just adding a parameter kineticsMap and the functions are exactly the same, but templetization of the map is required.

  std::string name_reactions_;
  std::vector<std::string> reaction_strings_;

  std::string name_reactions_heterogeneous_;
  std::vector<std::string> reaction_strings_heterogeneous_;

  // ---- Optional mechanism blocks ----
  // Extracted inside ReadKineticMechanism / ReadHeterogeneousKineticMechanism while
  // the parsing full kinetics.xml / kinetics.surface.xml  
  // Every block is optional and the matching flag says whether it was found.

  // <SpeciesClasses> (kinetics.xml). species_class_members_[c] = species indices in
  // class c; species_to_class_[s] = the class of species s, or -1 if unclassified.
  bool has_species_classes_;
  std::vector<std::string> species_class_names_;
  std::vector<std::vector<unsigned int>> species_class_members_;
  std::vector<int> species_to_class_;

  // <ReactionClasses> (gas: kinetics.xml, surface: kinetics.surface.xml). Both label
  // vectors are always sized to the phase's NumberOfReactions, "UNSORTED" by default.
  bool has_reaction_classes_;
  std::vector<std::string> reaction_main_class_;
  std::vector<std::string> reaction_sub_class_;
  bool has_reaction_classes_heterogeneous_;
  std::vector<std::string> reaction_main_class_heterogeneous_;
  std::vector<std::string> reaction_sub_class_heterogeneous_;

  // <SootProperties> (kinetics.xml) Empty when the mechanism carries no soot bins.
  // Direct parser for newer mechanisms which contain the BinProperties on the XML
  bool soot_available_;
  unsigned int soot_number_of_bins_;
  std::vector<unsigned int> soot_bin_index_;
  std::vector<int> soot_bin_section_;
  std::vector<double> soot_bin_nc_;
  std::vector<double> soot_bin_nh_;
  std::vector<double> soot_bin_no_;
  std::vector<double> soot_bin_htoc_;
  std::vector<double> soot_bin_mw_;
  std::vector<double> soot_bin_density_;
  std::vector<double> soot_bin_volume_;
  std::vector<double> soot_bin_mass_;
  std::vector<double> soot_bin_numpp_;
  std::vector<double> soot_bin_dsph_;
  std::vector<double> soot_bin_dcol_;
  std::vector<double> soot_bin_dpp_;
  std::vector<double> soot_bin_df_;

 private:
  void ReadSpeciesClassesBlock(const boost::property_tree::ptree& mechanism_ptree);
  void ReadReactionClassesBlock(const boost::property_tree::ptree& mechanism_ptree,
                                bool heterogeneous);
  void ReadSootProperties(const boost::property_tree::ptree& mechanism_ptree);
};

#include "ProfilesDatabase.hpp"
#endif  // PROFILESDATABASE_H
