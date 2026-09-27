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

#ifndef OUTPUTREADER_H
#define OUTPUTREADER_H

#include <Eigen/Sparse>
#include <unordered_map>

// Reader for the file Output.xml: T/P/MW/species profiles, <additional>
// columns, and species profiles.
// ReadOutput/Prepare/PrepareHeterogeneous never touch anything from the
// kinetics*.xml, which is what allows UpdateOutput to exist.
// PostProcessors may be built using Output.xml only if profiles are the
// required variable. It is recommended to build with kinetics anyway.
class OutputReader {
 public:
  OutputReader(void);
  ~OutputReader(void);

  bool ReadOutput(const std::string& folder_name, bool isHeterogeneous);

  // Re-reads a different Output.xml, reusing the heterogeneity mode of the
  // last ReadOutput call. Requires ReadOutput to have already succeeded once.
  bool UpdateOutput(const std::string& folder_name);

  void Prepare();
  void PrepareHeterogeneous();

  void SpeciesCoarsening(const double threshold);

  bool HasOutput() const { return is_output_available_; }

  const std::vector<double>& mwSpecies() const { return mw_species_; }
  const std::vector<std::string>& speciesNames() const { return species_names_unsorted_; }
  int SpeciesIndexFromName(const std::string& name) const;

  // Returns (independent_variable, profile) for one species. basis is "mass" or "mole"
  // for mass fractions or mole fractions.
  std::pair<std::vector<double>, std::vector<double>> GetSpeciesProfile(
      const std::string& name, const std::string& basis) const;

  // additional[0] column, corresponding to either time or axial coordinate.
  // This is the column you work with when calling local/region ROPA/Sensitivity.
  const std::vector<double>& GetIndependentVariableProfile() const;

  // Looks up a named <additional> column.
  // Shared implementation for functions like GetTemperatureProfile,
  // GetPressureProfile, etc.
  const std::vector<double>& AdditionalProfile(const std::string& key) const;

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

  std::vector<double> mw_species_;  
                // Note: for now, saving only the gas-phase molecular weights.
                // Since surface species are already saved in molar fractions, and
                // bulk species are in mass, but their activity is 1, no need to save
                // their MW (no use for it) Maybe, in different phases it is required
                // to save those.

  boost::property_tree::ptree xml_main_input;

  bool iSensitivityEnabled_;
  bool iSensitivityHeterogeneousEnabled_;

  // Folder of kinetic mechanism and output are the same for heterogeneous mechanisms, 
  // no need to duplicate.
  boost::filesystem::path path_folder_results_;

  // Remembered from the last ReadOutput call so UpdateOutput doesn't need the
  // caller to repeat it, and so PostProcessorCore can pass it consistently.
  bool is_mechanism_heterogeneous_;

 private:
  bool is_output_available_;

  // Species names in original (unsorted) declaration order - i.e. the same
  // order as omega[]/mw_species_ rows. string_list_massfractions_sorted above
  // is the *sorted* view used for display; this is the lookup-by-index one.
  std::vector<std::string> species_names_unsorted_;
  std::unordered_map<std::string, int> species_index_by_name_;

  // Shared by Prepare()/PrepareHeterogeneous() - both start by reading the same
  // <t-p-mw> leaf.
  void ReadTPMW();

  // Species sorting, moved to a function for reusability across
  // Prepare and PrepareHeterogeneous
  void SortAndIndexSpecies(const std::vector<std::string>& unsorted,
                           std::vector<std::string>& sorted,
                           std::vector<int>& sorted_index) const;
};

#include "OutputReader.hpp"
#endif  // OUTPUTREADER_H
