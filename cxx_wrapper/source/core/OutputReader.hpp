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

// Constructor
OutputReader::OutputReader(void) {
  index_density = -1;
  index_velocity = -1;
  index_mass_flow_rate = -1;

  iSensitivityEnabled_ = false;
  iSensitivityHeterogeneousEnabled_ = false;
  is_mechanism_heterogeneous_ = false;
  is_output_available_ = false;
}

// Destructor
OutputReader::~OutputReader(void) {}

bool OutputReader::ReadOutput(const std::string& folder_name, bool isHeterogeneous) {
  path_folder_results_ = folder_name;
  boost::filesystem::path path_results = path_folder_results_ / "Output.xml";
  // In all phases, the Output.xml contains the output for all phases
  // E.g. there is no Output.surface.xml like kinetics files do.

  if (!boost::filesystem::exists(path_results)) {
    throw std::invalid_argument("Output folder does not contain the Output.xml file");
  }

  // Cleared explicitly (not just default-constructed) because this same
  // object is reused across ReadOutput/UpdateOutput calls.
  xml_main_input.clear();
  boost::property_tree::read_xml((path_results).string(), xml_main_input);
  is_mechanism_heterogeneous_ = isHeterogeneous;

  if (is_mechanism_heterogeneous_ == true)
    PrepareHeterogeneous();
  else
    Prepare();

  boost::filesystem::path path_sensitivities = path_folder_results_ / "Sensitivities.xml";
  if (boost::filesystem::exists(path_sensitivities))
    iSensitivityEnabled_ = true;
  if (is_mechanism_heterogeneous_ == true)
  {
    boost::filesystem::path path_sensitivities_heterogeneous_ = path_folder_results_ / "Sensitivities.Surface.xml";
    if (boost::filesystem::exists(path_sensitivities_heterogeneous_))
      iSensitivityHeterogeneousEnabled_ = true;
  }

  is_output_available_ = true;
  return true;
}

bool OutputReader::UpdateOutput(const std::string& folder_name) {
  if (!is_output_available_) {
    throw std::invalid_argument(
        "UpdateOutput requires an output to already be loaded - call ReadOutput first");
  } // There is a point for not having this check. Let's see, it's ok for now.
  return ReadOutput(folder_name, is_mechanism_heterogeneous_);
}

void OutputReader::Prepare() {
  string_list_additional.clear();
  list_of_conversion_species_.clear();

  // Indices of T, P and MW
  {
    boost::optional<boost::property_tree::ptree&> child =
        xml_main_input.get_child_optional("opensmoke.t-p-mw");

    if (child) {
      std::stringstream stream;
      stream.str(xml_main_input.get<std::string>("opensmoke.t-p-mw"));
      stream >> index_T;
      stream >> index_P;
      stream >> index_MW;
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the t - p - mw leaf");
    }
  }

  {
    boost::optional<boost::property_tree::ptree&> child =
        xml_main_input.get_child_optional("opensmoke.additional");

    if (child) {
      std::stringstream stream;
      stream.str(xml_main_input.get<std::string>("opensmoke.additional"));

      unsigned int number_of_additional_profiles;
      stream >> number_of_additional_profiles;

      string_list_additional.reserve(number_of_additional_profiles);

      for (unsigned int j = 0; j < number_of_additional_profiles; j++) {
        std::string unit;
        std::string dummy;
        stream >> dummy;
        stream >> unit;
        string_list_additional.push_back(dummy + " " + unit);

        if (dummy == "density") index_density = j;
        if (dummy == "velocity") index_velocity = j;
        if (dummy == "mass-flow-rate") index_mass_flow_rate = j;
        if (dummy == "x-coord") index_x_coord = j;
        if (dummy == "z-coord") index_z_coord = j;
        if (dummy == "volume") index_volume = j;

        stream >> dummy;
      }
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the additional leaf");
    }
  }

  // Species (mass fractions)
  std::vector<std::string> string_list_massfractions_unsorted;
  {
    boost::optional<boost::property_tree::ptree&> child =
        xml_main_input.get_child_optional("opensmoke.mass-fractions");

    if (child) {
      std::stringstream stream;
      stream.str(xml_main_input.get<std::string>("opensmoke.mass-fractions"));

      unsigned int number_of_massfractions_profiles;
      stream >> number_of_massfractions_profiles;

      column_index_of_massfractions_profiles.resize(number_of_massfractions_profiles);
      string_list_massfractions_unsorted.reserve(number_of_massfractions_profiles);

      mw_species_.resize(number_of_massfractions_profiles);
      for (unsigned int j = 0; j < number_of_massfractions_profiles; j++) {
        std::string dummy;
        stream >> dummy;
        string_list_massfractions_unsorted.push_back(dummy);

        stream >> mw_species_[j];
        stream >> column_index_of_massfractions_profiles[j];
      }

      string_list_massfractions_sorted = string_list_massfractions_unsorted;

      std::sort(string_list_massfractions_sorted.begin(), string_list_massfractions_sorted.end());
      // string_list_massfractions_sorted.sort();

      sorted_index.resize(number_of_massfractions_profiles);
      for (unsigned int j = 0; j < number_of_massfractions_profiles; j++)
        for (unsigned int k = 0; k < number_of_massfractions_profiles; k++)
          if (string_list_massfractions_sorted[j] == string_list_massfractions_unsorted[k]) {
            sorted_index[j] = k;
            break;
          }
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the mass-fractions leaf");
    }
  }

  species_names_unsorted_ = string_list_massfractions_unsorted;
  species_index_by_name_.clear();
  for (unsigned int j = 0; j < species_names_unsorted_.size(); j++)
    species_index_by_name_[species_names_unsorted_[j]] = static_cast<int>(j);

  // Read profiles
  omega.resize(column_index_of_massfractions_profiles.size());
  additional.resize(string_list_additional.size());
  {
    {
      boost::optional<boost::property_tree::ptree&> child =
          xml_main_input.get_child_optional("opensmoke.profiles-size");

      if (child) {
        std::stringstream stream;
        stream.str(xml_main_input.get<std::string>("opensmoke.profiles-size"));
        stream >> number_of_abscissas_;
        stream >> number_of_ordinates_;
      } else {
        throw std::invalid_argument("Corrupted xml file: missing the profiles-size leaf");
      }
    }

    omega.resize(column_index_of_massfractions_profiles.size());
    for (unsigned int j = 0; j < column_index_of_massfractions_profiles.size(); j++)
      omega[j].resize(number_of_abscissas_);

    additional.resize(string_list_additional.size());
    for (unsigned int j = 0; j < string_list_additional.size(); j++)
      additional[j].resize(number_of_abscissas_);

    boost::optional<boost::property_tree::ptree&> child =
        xml_main_input.get_child_optional("opensmoke.profiles");
    if (child) {
      std::stringstream stream;
      stream.str(xml_main_input.get<std::string>("opensmoke.profiles"));

      for (unsigned int i = 0; i < number_of_abscissas_; i++) {
        for (unsigned int j = 0; j < string_list_additional.size(); j++) stream >> additional[j][i];
        for (unsigned int j = 0; j < column_index_of_massfractions_profiles.size(); j++)
          stream >> omega[j][i];
      }
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the profiles leaf");
    }

    sorted_max.resize(string_list_massfractions_sorted.size());
    for (unsigned int j = 0; j < column_index_of_massfractions_profiles.size(); j++) {
      sorted_max[j] = -1.e100;
      for (unsigned int i = 0; i < number_of_abscissas_; i++)
        if (omega[sorted_index[j]][i] > sorted_max[j]) sorted_max[j] = omega[sorted_index[j]][i];
    }
  }

  // Conversions
  {
    for (unsigned int j = 0; j < column_index_of_massfractions_profiles.size(); j++) {
      if (omega[j][0] > 1e-8) {
        list_of_conversion_species_.push_back(j);
        string_list_additional.push_back("conversion-" + string_list_massfractions_unsorted[j]);
        std::vector<double> tmp(number_of_abscissas_);
        for (unsigned int i = 0; i < number_of_abscissas_; i++)
          tmp[i] = (omega[j][0] - omega[j][i]) / omega[j][0];
        additional.push_back(tmp);
      }
    }
  }
}

void OutputReader::PrepareHeterogeneous() {
  string_list_additional.clear();
  list_of_conversion_species_.clear();

  // this part is the same for heterogeneous kinetics
  // Indices of T, P and MW
  {
    boost::optional<boost::property_tree::ptree&> child =
        xml_main_input.get_child_optional("opensmoke.t-p-mw");

    if (child) {
      std::stringstream stream;
      stream.str(xml_main_input.get<std::string>("opensmoke.t-p-mw"));
      stream >> index_T;
      stream >> index_P;
      stream >> index_MW;
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the t - p - mw leaf");
    }
  }

  // Additional
  // In heterogeneous phases, we need more/different additional parameters
  // For instance, in the surface phase we need A/V and the site density
  {
    boost::optional<boost::property_tree::ptree&> child =
        xml_main_input.get_child_optional("opensmoke.additional");

    if (child) {
      std::stringstream stream;
      stream.str(xml_main_input.get<std::string>("opensmoke.additional"));

      unsigned int number_of_additional_profiles;
      stream >> number_of_additional_profiles;

      string_list_additional.reserve(number_of_additional_profiles);

      for (unsigned int j = 0; j < number_of_additional_profiles; j++) {
        std::string unit;
        std::string dummy;
        stream >> dummy;
        stream >> unit;
        string_list_additional.push_back(dummy + " " + unit);

        if (dummy == "density") index_density = j;
        if (dummy == "velocity") index_velocity = j;
        if (dummy == "mass-flow-rate") index_mass_flow_rate = j;
        if (dummy == "x-coord") index_x_coord = j;
        if (dummy == "z-coord") index_z_coord = j;
        if (dummy == "volume") index_volume = j;
        if (dummy == "area-over-volume") index_area_over_volume = j;
        if (dummy == "CARBON") index_surface_sites_concentration = j;
        /*
        LG: About the "CARBON" name: I added its print in the OpenSMOKEpp library,
        and it is the sites concentration in kmol/m2.
        For C-deposition, it is ok, while for catalytic-type problems, it is not: 
        in the input, the surface would be called something like "Surface-NI", 
        and the printed name would then be just "NI". 
        It might be better to look not for the exact name but to print the "Surface-" 
        keyword and look for anything that starts with that.
        Further note: in OpenSMOKEpp, AC allowed to have multiple surface phases.
                      This is not considered here, and also not important at the moment.
        TODO: fix this better.
        */

        stream >> dummy;
      }
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the additional leaf");
    }
  }

  // Species (gas mass fractions)
  // Note the different name wrt homogeneous (mass-fractions vs gas-mass-fractions)
  std::vector<std::string> string_list_massfractions_unsorted;
  {
    boost::optional<boost::property_tree::ptree&> child =
        xml_main_input.get_child_optional("opensmoke.gas-mass-fractions");

    if (child) {
      std::stringstream stream;
      stream.str(xml_main_input.get<std::string>("opensmoke.gas-mass-fractions"));

      unsigned int number_of_massfractions_profiles;
      stream >> number_of_massfractions_profiles;
      number_of_gas_species = number_of_massfractions_profiles;

      column_index_of_massfractions_profiles.resize(number_of_massfractions_profiles);
      string_list_massfractions_unsorted.reserve(number_of_massfractions_profiles);

      mw_species_.resize(number_of_massfractions_profiles);
      for (unsigned int j = 0; j < number_of_massfractions_profiles; j++) {
        std::string dummy;
        stream >> dummy;
        string_list_massfractions_unsorted.push_back(dummy);

        stream >> mw_species_[j];
        stream >> column_index_of_massfractions_profiles[j];
      }

      string_list_massfractions_sorted = string_list_massfractions_unsorted;

      std::sort(string_list_massfractions_sorted.begin(), string_list_massfractions_sorted.end());
      // string_list_massfractions_sorted.sort();

      sorted_index.resize(number_of_massfractions_profiles);
      for (unsigned int j = 0; j < number_of_massfractions_profiles; j++)
        for (unsigned int k = 0; k < number_of_massfractions_profiles; k++)
          if (string_list_massfractions_sorted[j] == string_list_massfractions_unsorted[k]) {
            sorted_index[j] = k;
            break;
          }
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the gas-mass-fractions leaf");
    }
  }

  species_names_unsorted_ = string_list_massfractions_unsorted;
  species_index_by_name_.clear();
  for (unsigned int j = 0; j < species_names_unsorted_.size(); j++)
    species_index_by_name_[species_names_unsorted_[j]] = static_cast<int>(j);

  // Enhancement: create a sorting function and pass all of these "subsections" to that function sequentially. (Same story for the homogeneous call)
  // Species (surface fractions)
  std::vector<std::string> string_list_surfacefractions_unsorted;
  {
    boost::optional<boost::property_tree::ptree&> child =
        xml_main_input.get_child_optional("opensmoke.surface-moles-fractions");

    if (child) {
      std::stringstream stream;
      stream.str(xml_main_input.get<std::string>("opensmoke.surface-moles-fractions"));

      unsigned int number_of_surfacefractions_profiles;
      stream >> number_of_surfacefractions_profiles;
      number_of_surface_species = number_of_surfacefractions_profiles;

      column_index_of_surfacefractions_profiles.resize(number_of_surfacefractions_profiles);
      string_list_surfacefractions_unsorted.reserve(number_of_surfacefractions_profiles);

      for (unsigned int j = 0; j < number_of_surfacefractions_profiles; j++) {
        std::string dummy;
        stream >> dummy;
        string_list_surfacefractions_unsorted.push_back(dummy);

        stream >> dummy;
        stream >> column_index_of_surfacefractions_profiles[j];
      }
      string_list_surfacefractions_sorted = string_list_surfacefractions_unsorted;

      std::sort(string_list_surfacefractions_sorted.begin(), string_list_surfacefractions_sorted.end());

      sorted_index_surface.resize(number_of_surfacefractions_profiles);
      for (unsigned int j = 0; j < number_of_surfacefractions_profiles; j++)
        for (unsigned int k = 0; k < number_of_surfacefractions_profiles; k++)
          if (string_list_surfacefractions_sorted[j] == string_list_surfacefractions_unsorted[k]) {
            sorted_index_surface[j] = k;
            break;
        }
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the surface-moles-fractions leaf");
    }
  }
  // Species (bulk masses)
  std::vector<std::string> string_list_bulkmasses_unsorted;
  {
    boost::optional<boost::property_tree::ptree&> child =
        xml_main_input.get_child_optional("opensmoke.bulk-masses");

    if (child) {
      std::stringstream stream;
      stream.str(xml_main_input.get<std::string>("opensmoke.bulk-masses"));

      unsigned int number_of_bulkmasses_profiles;
      stream >> number_of_bulkmasses_profiles;
      number_of_bulk_species = number_of_bulkmasses_profiles;

      column_index_of_bulkmasses_profiles.resize(number_of_bulkmasses_profiles);
      string_list_bulkmasses_unsorted.reserve(number_of_bulkmasses_profiles);

      for (unsigned int j = 0; j < number_of_bulkmasses_profiles; j++) {
        std::string dummy;
        stream >> dummy;
        string_list_bulkmasses_unsorted.push_back(dummy);

        stream >> dummy;
        stream >> column_index_of_bulkmasses_profiles[j];
      }

      string_list_bulkmasses_sorted = string_list_bulkmasses_unsorted;

      std::sort(string_list_bulkmasses_sorted.begin(), string_list_bulkmasses_sorted.end());
      // string_list_massfractions_sorted.sort();

      sorted_index_bulk.resize(number_of_bulkmasses_profiles);
      for (unsigned int j = 0; j < number_of_bulkmasses_profiles; j++)
        for (unsigned int k = 0; k < number_of_bulkmasses_profiles; k++)
          if (string_list_bulkmasses_sorted[j] == string_list_bulkmasses_unsorted[k]) {
            sorted_index_bulk[j] = k;
            break;
          }
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the bulk-masses leaf");
    }
  }

  // Read profiles
  omega.resize(number_of_gas_species);
  Z.resize(number_of_surface_species);
  massBulk.resize(number_of_bulk_species);
  additional.resize(string_list_additional.size());
  {
    {
      boost::optional<boost::property_tree::ptree&> child =
          xml_main_input.get_child_optional("opensmoke.profiles-size");

      if (child) {
        std::stringstream stream;
        stream.str(xml_main_input.get<std::string>("opensmoke.profiles-size"));
        stream >> number_of_abscissas_;
        stream >> number_of_ordinates_;
      } else {
        throw std::invalid_argument("Corrupted xml file: missing the profiles-size leaf");
      }
    }

    omega.resize(number_of_gas_species);
    for (unsigned int j = 0; j < number_of_gas_species; j++)
      omega[j].resize(number_of_abscissas_);

    Z.resize(number_of_surface_species);
    for (unsigned int j = 0; j < number_of_surface_species; j++)
      Z[j].resize(number_of_abscissas_);

    massBulk.resize(number_of_bulk_species);
    for (unsigned int j = 0; j < number_of_bulk_species; j++)
      massBulk[j].resize(number_of_abscissas_);

    additional.resize(string_list_additional.size());
    for (unsigned int j = 0; j < string_list_additional.size(); j++)
      additional[j].resize(number_of_abscissas_);

    boost::optional<boost::property_tree::ptree&> child =
        xml_main_input.get_child_optional("opensmoke.profiles");
    if (child) {
      std::stringstream stream;
      stream.str(xml_main_input.get<std::string>("opensmoke.profiles"));

      for (unsigned int i = 0; i < number_of_abscissas_; i++) {
        for (unsigned int j = 0; j < string_list_additional.size(); j++) stream >> additional[j][i];
        for (unsigned int j = 0; j < column_index_of_massfractions_profiles.size(); j++)
          stream >> omega[j][i];
        for (unsigned int j = 0; j < column_index_of_surfacefractions_profiles.size(); j++)
          stream >> Z[j][i];
        for (unsigned int j = 0; j < column_index_of_bulkmasses_profiles.size(); j++)
          stream >> massBulk[j][i];
      }
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the profiles leaf");
    }

    sorted_max.resize(string_list_massfractions_sorted.size());
    for (unsigned int j = 0; j < column_index_of_massfractions_profiles.size(); j++) {
      sorted_max[j] = -1.e100;
      for (unsigned int i = 0; i < number_of_abscissas_; i++)
        if (omega[sorted_index[j]][i] > sorted_max[j]) sorted_max[j] = omega[sorted_index[j]][i];
    }

    sorted_max_surface.resize(string_list_surfacefractions_sorted.size());
    for (unsigned int j = 0; j < column_index_of_surfacefractions_profiles.size(); j++) {
      sorted_max_surface[j] = -1.e100;
      for (unsigned int i = 0; i < number_of_abscissas_; i++)
        if (Z[sorted_index_surface[j]][i] > sorted_max_surface[j]) sorted_max_surface[j] = Z[sorted_index_surface[j]][i];
    }

    sorted_max_bulk.resize(string_list_bulkmasses_sorted.size());
    for (unsigned int j = 0; j < column_index_of_bulkmasses_profiles.size(); j++) {
      sorted_max_bulk[j] = -1.e100;
      for (unsigned int i = 0; i < number_of_abscissas_; i++)
        if (massBulk[sorted_index_bulk[j]][i] > sorted_max_bulk[j]) sorted_max_bulk[j] = massBulk[sorted_index_bulk[j]][i];
    }
  }

  // Reactants conversion (probably not required in general)
  // LG This is wrong as there is no mass loss correction. 
  // TODO fix this (not important)
  {
    for (unsigned int j = 0; j < number_of_gas_species; j++) {
      if (omega[j][0] > 1e-8) {
        list_of_conversion_species_.push_back(j);
        string_list_additional.push_back("conversion-" + string_list_massfractions_unsorted[j]);
        std::vector<double> tmp(number_of_abscissas_);
        for (unsigned int i = 0; i < number_of_abscissas_; i++)
          tmp[i] = (omega[j][0] - omega[j][i]) / omega[j][0];
        additional.push_back(tmp);
      }
    }
  }
}

// LG Unused function, probably old, can we delete this?
void OutputReader::SpeciesCoarsening(const double threshold) {
  current_sorted_index.resize(0);
  for (unsigned int k = 0; k < string_list_massfractions_sorted.size(); k++)
    if (sorted_max[k] > threshold) current_sorted_index.push_back(k);
}

int OutputReader::SpeciesIndexFromName(const std::string& name) const {
  std::unordered_map<std::string, int>::const_iterator it = species_index_by_name_.find(name);
  if (it == species_index_by_name_.end()) {
    throw std::invalid_argument("Unknown species: " + name);
  }
  return it->second;
}

std::pair<std::vector<double>, std::vector<double>> OutputReader::GetSpeciesProfile(
    const std::string& name, const std::string& basis) const {
  const int idx = SpeciesIndexFromName(name);

  std::vector<double> x = GetIndependentVariableProfile();
  std::vector<double> y_mass = omega[idx];

  if (basis == "mass") {
    return std::make_pair(x, y_mass);
  }
  if (basis == "moles") {
    std::vector<double> y_mole(y_mass.size());
    const std::vector<double>& mw_mix = additional[index_MW];
    for (unsigned int i = 0; i < y_mass.size(); i++) {
      y_mole[i] = y_mass[i] * mw_mix[i] / mw_species_[idx];
    }
    return std::make_pair(x, y_mole);
  }
  throw std::invalid_argument("basis must be \"mass\" or \"moles\", got \"" + basis + "\"");
}

const std::vector<double>& OutputReader::GetIndependentVariableProfile() const {
  return additional[0];
}

const std::vector<double>& OutputReader::AdditionalProfile(const std::string& key) const {
  for (unsigned int j = 0; j < string_list_additional.size(); j++) {
    if (string_list_additional[j].find(key) != std::string::npos) {
      return additional[j];
    }
  }
  throw std::invalid_argument("No <additional> column matching \"" + key + "\"");
}
