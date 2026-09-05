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

Sensitivities_Database_Surface::Sensitivities_Database_Surface(
    bool heterogeneousSensitivity)
    : Sensitivities_Database(), heterogeneousSensitivity_(heterogeneousSensitivity) {}

void Sensitivities_Database_Surface::ReadParentFile() {
  boost::filesystem::path path_results;
  if (heterogeneousSensitivity_ ==
      true)  // We should look for this using the phase_name, not implemented here yet.
    path_results = data_->path_folder_results_ / "Sensitivities.Surface.xml";
  else
    path_results = data_->path_folder_results_ / "Sensitivities.xml";
  boost::property_tree::read_xml((path_results).string(), xml_main_input);

  {
    std::stringstream stream;
    stream.str(xml_main_input.get<std::string>("opensmoke.variables"));

    stream >> number_of_variables_;

    names_.resize(number_of_variables_);
    local_index_.resize(number_of_variables_);
    global_index_.resize(number_of_variables_);
    for (unsigned int j = 0; j < number_of_variables_; j++) {
      stream >> names_[j];
      stream >> local_index_[j];
      stream >> global_index_[j];
    }
  }

  {
    std::stringstream stream;
    stream.str(xml_main_input.get<std::string>("opensmoke.n-parameters"));

    stream >> number_of_parameters_;
    parameters_.resize(number_of_parameters_);
  }

  {
    std::stringstream stream;
    stream.str(xml_main_input.get<std::string>("opensmoke.constant-parameters"));

    for (unsigned int j = 0; j < number_of_parameters_; j++) stream >> parameters_[j];
  }

  {
    std::stringstream stream;
    stream.str(xml_main_input.get<std::string>("opensmoke.points"));

    stream >> number_of_points_;
  }

  // Memory allocation
  coefficients_.resize(number_of_parameters_);
  for (unsigned int j = 0; j < number_of_parameters_; j++)
    coefficients_[j].resize(number_of_points_);
  variable_.resize(number_of_points_);

  // Names of reactions
  {
    unsigned int NR;
    std::vector<std::string> reaction_strings;
    if (heterogeneousSensitivity_ == true) {
      NR = data_->kineticsMapSurfaceXML->NumberOfReactions();
      reaction_strings = data_->reaction_strings_heterogeneous_;
    } else {
      NR = data_->kineticsMapXML->NumberOfReactions();
      reaction_strings = data_->reaction_strings_;
    }

    string_list_reactions_.reserve(number_of_parameters_);

    for (unsigned int j = 0; j < number_of_parameters_; j++) {
      if (j + 1 <= NR) {
        std::stringstream index;
        index << j + 1;
        std::string tmp = "R" + index.str() + ": " + reaction_strings[j];
        string_list_reactions_.push_back(tmp);
      } else {
        unsigned int local_index = j + 1 - NR;

        if (heterogeneousSensitivity_ == false) {
          if (local_index <= data_->kineticsMapXML->NumberOfFallOffReactions()) {
            unsigned int global_index =
                data_->kineticsMapXML->IndicesOfFalloffReactions()[local_index - 1];
            std::stringstream index;
            index << global_index;
            std::string tmp = "R" + index.str() +
                              "(inf): " + data_->reaction_strings_[global_index - 1];
            string_list_reactions_.push_back(tmp);
          } else {
            unsigned int global_index =
                data_->kineticsMapXML->IndicesOfCabrReactions()
                    [local_index - data_->kineticsMapXML->NumberOfFallOffReactions() - 1];
            std::stringstream index;
            index << global_index;
            std::string tmp = "R" + index.str() +
                              "(inf): " + data_->reaction_strings_[global_index - 1];
            string_list_reactions_.push_back(tmp);
          }
        }
      }
    }
  }
}

void Sensitivities_Database_Surface::ReadFromChildFile(const std::string name) {
  std::string local_name;
  if (heterogeneousSensitivity_ == true)  
  // As above, we should look for the name using the phase_name variable
    local_name = "Sensitivities.Surface." + name + ".xml";
  else
    local_name = "Sensitivities." + name + ".xml";

  boost::property_tree::ptree ptree;
  try {
    boost::filesystem::path path_file = data_->path_folder_results_ / local_name;
    boost::property_tree::read_xml((path_file).string(), ptree);
  } catch (...) {
    std::cout << " The file named: " << local_name
              << " it is not present in the output folder select another specie."
              << std::endl;
    exit(-1);
  }

  boost::optional<boost::property_tree::ptree&> child =
      ptree.get_child_optional("opensmoke.coefficients");

  if (!child) {
    std::cout << "Corrupted xml file: missing the coefficients leaf" << std::endl;
    exit(-1);
  } else {
    std::stringstream stream;
    stream.str(ptree.get<std::string>("opensmoke.coefficients"));

    for (unsigned int i = 0; i < number_of_points_; i++)
      for (unsigned int j = 0; j < number_of_parameters_; j++)
        stream >> coefficients_[j][i];
  }

  for (unsigned int j = 0; j < number_of_variables_; j++)
    if (names_[j] == name) {
      current_local_index_ = j;
      break;
    }

  if (name == "temperature")
    variable_ = data_->additional[data_->index_T];
  else if (name == "density")
    variable_ = data_->additional[data_->index_density];
  else if (name == "velocity")
    variable_ = data_->additional[data_->index_velocity];
  else if (name == "mass-flow-rate")
    variable_ = data_->additional[data_->index_mass_flow_rate];
  else {
    // recognizing the global index of species
    unsigned int species_index;
    // IMPORTANT NOTE: we should put a try/catch here, right now I'm only performing
    // sensitivity on gas-phase species, I'm surely not interested in surface ones, but I
    // think I am actually interested in the bulk masses sensitivities. So - eventually -
    // we will deal with exceptions on the names
    for (unsigned int j = 0; j < data_->string_list_massfractions_sorted.size(); j++)
      if (name == data_->string_list_massfractions_sorted[j]) {
        species_index = data_->sorted_index[j];
        break;
      }

    // from mass to mole fractions (only performed on gas-phase species)
    for (unsigned int i = 0; i < number_of_points_; i++)
      variable_[i] = data_->omega[species_index][i] *
                     data_->additional[data_->index_MW][i] /
                     data_->mw_species_[species_index];
  }
}
