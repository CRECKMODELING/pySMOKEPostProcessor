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

ProfilesDatabase::ProfilesDatabase(void) {
  iSensitivityEnabled_ = false;
  iROPAEnabled_ = false;
  is_kinetics_available_ = false;

  is_mechanism_heterogeneous_ = false;
  iROPAHeterogeneousEnabled_ = false;
  iSensitivityHeterogeneousEnabled_ = false;
  is_heterogeneous_kinetics_available_ = false;

  has_species_classes_ = false;
  has_reaction_classes_ = false;
  has_reaction_classes_heterogeneous_ = false;
  soot_available_ = false;
  soot_number_of_bins_ = 0;

  index_density = -1;
  index_velocity = -1;
  index_mass_flow_rate = -1;
}

ProfilesDatabase::~ProfilesDatabase(void) {}

bool ProfilesDatabase::ReadKineticMechanism(const std::string& folder_name) {
  path_folder_mechanism_ = folder_name;
  boost::filesystem::path path_mechanism = path_folder_mechanism_ / "kinetics.xml";

  if (!boost::filesystem::exists(path_mechanism)) {
    throw std::invalid_argument(
        "The folder of the kinetic mechanism does not contains any kinetics.xml!"
    );
  }

  // Read from file
  {
    boost::property_tree::ptree ptree;
    boost::property_tree::read_xml((path_mechanism).string(), ptree);

    thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree, false);
    kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, ptree, false);

    // Optional post-processing blocks
    ReadSpeciesClassesBlock(ptree);
    ReadReactionClassesBlock(ptree, false);
    ReadSootProperties(ptree);
  }

  if (thermodynamicsMapXML->NumberOfSpecies() == omega.size()) {
    iROPAEnabled_ = true;
  } else {
    throw std::invalid_argument(
        "Output.xml file contains only a subset of the total species in the kinetic mechanism"
    );
  }

  // Read the reaction strings
  {
    std::string local_name = "reaction_names.xml";
    boost::filesystem::path path_reaction_names = path_folder_mechanism_ / local_name;

    if (!boost::filesystem::exists(path_reaction_names)) {
      throw std::invalid_argument("Kinetic folder does not contain the reaction_names.xml file");
    }

    {
      boost::property_tree::ptree ptree;
      boost::property_tree::read_xml((path_reaction_names).string(), ptree);

      // Names of reactions
      {
        std::stringstream stream;
        stream.str(ptree.get<std::string>("opensmoke.reaction-names"));

        reaction_strings_.reserve(kineticsMapXML->NumberOfReactions());
        for (unsigned int j = 0; j < kineticsMapXML->NumberOfReactions(); j++) {
          std::string reaction_string;
          stream >> reaction_string;
          reaction_strings_.push_back(reaction_string);
        }
      }
    }
  }

  is_kinetics_available_ = true;

  return true;
}

bool ProfilesDatabase::ReadHeterogeneousKineticMechanism(const std::string& folder_name, const std::string& phase_name) {
  // phase name is not used right now, but it can be in the future (e.g. phase_name == Surface or Liquid etc)
  is_mechanism_heterogeneous_ = ReadKineticMechanism(folder_name); // The homogeneous kinetics.xml must be read regardless of heterogeneous.
  path_folder_mechanism_ = folder_name;

  // Using phase_name we can understand the name of the file. kinetics.?.xml
  boost::filesystem::path path_mechanism_het = path_folder_mechanism_ / "kinetics.surface.xml";

  if (!boost::filesystem::exists(path_mechanism_het)) {
    throw std::invalid_argument(
        "The folder of the kinetic mechanism does not contains the heterogeneous kinetics.xml!"
    );
  }

  // Read from file
  {
    boost::property_tree::ptree ptree_het;
    boost::property_tree::read_xml((path_mechanism_het).string(), ptree_het);

    // These names will be changed to heterogeneous maybe (notice in the surface map the bool is not required, this is general also to liquid and solid)
    
    std::streambuf* old_buf = std::cout.rdbuf();   // When calling the thermoSurfaceMap, the code prints a lot of (useless) stuff;
    std::ofstream null_stream("/dev/null");        // this is to deal with it.
    std::cout.rdbuf(null_stream.rdbuf());          // Basically: sends buffer to null, but only for a second

    thermodynamicsMapSurfaceXML = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN(ptree_het);
    kineticsMapSurfaceXML = new OpenSMOKE::KineticsMap_Surface_CHEMKIN(*thermodynamicsMapSurfaceXML, ptree_het);

    std::cout.rdbuf(old_buf);                      // Restores the buffer to the original (goes back to printing on screen)

    // Surface <ReactionClasses>
    ReadReactionClassesBlock(ptree_het, true);
  }

  // Disabling this check because number of species is no longer the same (e.g. in surface I have surface fractions etc)
  // In my case the correct check is NumberOfSpecies == omega + Z + massBulk but I don't know about the other phases
  // if (thermodynamicsSurfaceMapXML->NumberOfSpecies() == omega.size()) {
  //   iROPAEnabled_ = true;
  // } else {
  //   throw std::invalid_argument(
  //       "Output.xml file contains only a subset of the total species in the kinetic mechanism"
  //   );
  // }
  iROPAHeterogeneousEnabled_ = true;

  // Read the reaction strings
  { // The reaction_names file has a different name depending on phase. 
    // For me it's that, for liquid its reaction_names.liquid.xml, for solid no file exists (@Riccardo)
    // Again the phase_name arg should be used to define this variable properly.
    std::string local_name_het = "surface_reaction_names.xml"; 
    boost::filesystem::path path_reaction_names_het = path_folder_mechanism_ / local_name_het;
    if (!boost::filesystem::exists(path_reaction_names_het)) {
      throw std::invalid_argument("Kinetic folder does not contain the heterogeneous reaction_names.xml file");
    }

    {
      boost::property_tree::ptree ptree_het;
      boost::property_tree::read_xml((path_reaction_names_het).string(), ptree_het);

      // Names of reactions
      {
        std::stringstream stream;
        stream.str(ptree_het.get<std::string>("opensmoke.reaction-names"));

        reaction_strings_heterogeneous_.reserve(kineticsMapSurfaceXML->NumberOfReactions());
        for (unsigned int j = 0; j < kineticsMapSurfaceXML->NumberOfReactions(); j++) {
          std::string reaction_string;
          stream >> reaction_string;
          reaction_strings_heterogeneous_.push_back(reaction_string);
        }
      }
    }
  }

  is_heterogeneous_kinetics_available_ = true;
  return true;
}

// Locates an optional post-processing block in a pre-existing ptree.
// The block location is phase-specific: in gas-phase 
// <SpeciesClasses>/<ReactionClasses> are under <opensmoke><Kinetics>; 
// surface mechanism blocks are in  <opensmoke><Kinetics><MaterialKinetics>.
// <SootProperties> is a direct <opensmoke> child and is read separately.
static boost::optional<const boost::property_tree::ptree&> FindMechanismBlock(
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
// Fills species_class_names_ / species_class_members_ and the per-species class index
// species_to_class_ (-1 = unclassified). The block is a flat list of
// <ClassSpecies name="..."> nodes whose text is whitespace-separated species indices.
void ProfilesDatabase::ReadSpeciesClassesBlock(
    const boost::property_tree::ptree& mechanism_ptree) {
  // Init
  has_species_classes_ = false;
  species_class_names_.clear();
  species_class_members_.clear();
  const unsigned int ns = thermodynamicsMapXML->NumberOfSpecies();
  species_to_class_.assign(ns, -1);

  boost::optional<const boost::property_tree::ptree&> block =
      FindMechanismBlock(mechanism_ptree, "SpeciesClasses");
  if (!block) return; // Block not found: return

  for (const boost::property_tree::ptree::value_type& child : *block) {
    if (child.first != "ClassSpecies") continue;

    const std::string name = child.second.get<std::string>(
        "<xmlattr>.name", "class_" + std::to_string(species_class_names_.size()));

    std::vector<unsigned int> members;
    std::istringstream body(child.second.data());
    long idx;
    while (body >> idx) {
      if (idx < 0 || static_cast<unsigned int>(idx) >= ns) continue;
      species_to_class_[idx] = static_cast<int>(species_class_names_.size());
      members.push_back(static_cast<unsigned int>(idx));
    }

    species_class_names_.push_back(name);
    species_class_members_.push_back(members);
  }

  has_species_classes_ = !species_class_names_.empty();
}

// --- optional <ReactionClasses> block ----------------------------------
// Same idea for reactions. <MainClass name> -> <SubClass name> -> <ReactionIndices>
// (whitespace-separated). The label vectors are always sized to the phase's
// NumberOfReactions, "UNSORTED" where a reaction has no entry.
void ProfilesDatabase::ReadReactionClassesBlock(
    const boost::property_tree::ptree& mechanism_ptree, const bool heterogeneous) {
  std::vector<std::string>* main_class = &reaction_main_class_;
  std::vector<std::string>* sub_class = &reaction_sub_class_;
  bool* available = &has_reaction_classes_;
  unsigned int nr = kineticsMapXML->NumberOfReactions();
  if (heterogeneous) {
    main_class = &reaction_main_class_heterogeneous_;
    sub_class = &reaction_sub_class_heterogeneous_;
    available = &has_reaction_classes_heterogeneous_;
    nr = kineticsMapSurfaceXML->NumberOfReactions();
  }

  *available = false;
  main_class->assign(nr, "UNSORTED");
  sub_class->assign(nr, "UNSORTED");

  boost::optional<const boost::property_tree::ptree&> block =
      FindMechanismBlock(mechanism_ptree, "ReactionClasses");
  if (!block) return;

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
        (*main_class)[ridx] = main_name;
        (*sub_class)[ridx] = sub_name;
      }
    }
  }

  *available = true;
}

// --- optional <SootProperties> block ---------------------------------
// PolimiSoot BIN properties, one value per bin per field. This is just the flat
// array-read half of PolimiSoot_Analyzer::BinPropertiesFromXMLFile
// (utilities/soot/polimi), transcribed against the kinetics.xml tree already in
// memory: no re-parse, and none of that method's class-taxonomy-dependent bin
// grouping (which we don't expose). No block -> soot_available_ stays false and the
// vectors stay empty; the preprocessor writes the block whole, so a block that is
// present but missing a field is a corrupt mechanism and throws like the rest of
// ReadKineticMechanism.
void ProfilesDatabase::ReadSootProperties(
    const boost::property_tree::ptree& mechanism_ptree) {
  soot_available_ = false;
  soot_number_of_bins_ = 0;

  boost::optional<const boost::property_tree::ptree&> block =
      mechanism_ptree.get_child_optional("opensmoke.SootProperties");
  if (!block) return;

  {
    std::stringstream s(block->get<std::string>("NumberOfBins"));
    s >> soot_number_of_bins_;
  }

  const auto read_uint = [&](const char* key, std::vector<unsigned int>& v) {
    std::stringstream s(block->get<std::string>(key));
    v.assign(soot_number_of_bins_, 0u);
    for (unsigned int k = 0; k < soot_number_of_bins_; k++) s >> v[k];
  };
  const auto read_int = [&](const char* key, std::vector<int>& v) {
    std::stringstream s(block->get<std::string>(key));
    v.assign(soot_number_of_bins_, 0);
    for (unsigned int k = 0; k < soot_number_of_bins_; k++) s >> v[k];
  };
  const auto read_double = [&](const char* key, std::vector<double>& v) {
    std::stringstream s(block->get<std::string>(key));
    v.assign(soot_number_of_bins_, 0.);
    for (unsigned int k = 0; k < soot_number_of_bins_; k++) s >> v[k];
  };

  read_uint("Bin_index", soot_bin_index_);
  read_int("Bin_section", soot_bin_section_);
  read_double("Bin_nc", soot_bin_nc_);
  read_double("Bin_nh", soot_bin_nh_);
  read_double("Bin_no", soot_bin_no_);
  read_double("Bin_htoc", soot_bin_htoc_);
  read_double("Bin_mw", soot_bin_mw_);
  read_double("Bin_density", soot_bin_density_);
  read_double("Bin_volume", soot_bin_volume_);
  read_double("Bin_mass", soot_bin_mass_);
  read_double("Bin_numpp", soot_bin_numpp_);
  read_double("Bin_dsph", soot_bin_dsph_);
  read_double("Bin_dcol", soot_bin_dcol_);
  read_double("Bin_dpp", soot_bin_dpp_);
  read_double("Bin_df", soot_bin_df_);

  soot_available_ = true;
}

bool ProfilesDatabase::ReadFileResults(const std::string& folder_name, bool isHeterogeneous) {
  path_folder_results_ = folder_name;
  boost::filesystem::path path_results = path_folder_results_ / "Output.xml";   // Results are all in the same output.xml file (same as homogeneous)

  if (!boost::filesystem::exists(path_results)) {
    throw std::invalid_argument("Output folder does not contain the Output.xml file");
  }

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
  

  return true;
}

void ProfilesDatabase::Prepare() {
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

void ProfilesDatabase::PrepareHeterogeneous() {
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
        About the "CARBON" name: I added its print in the OpenSMOKEpp library,and it is the sites concentration in kmol/m2. 
        For C-deposition, it is ok, while for catalytic-type problems, it is not: in the input, the surface would be called
        something like "Surface-NI", and the printed name would then be just "NI". It might be better to look not for the
        exact name but to print the "Surface-" keyword and look for anything that starts with that.
        Further note: in OpenSMOKEpp, AC allowed to have multiple surface phases. 
                      This is not considered here, and also not important at the moment.
        */

        stream >> dummy;
      }
    } else {
      throw std::invalid_argument("Corrupted xml file: missing the additional leaf");
    }
  }

  // Species (gas mass fractions)
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
  // This is wrong as there is no mass loss correction. Not a priority right now but it has to be changed. I'll do it eventually.
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

// Unused function, probably old, can we delete this?
void ProfilesDatabase::SpeciesCoarsening(const double threshold) {
  current_sorted_index.resize(0);
  for (unsigned int k = 0; k < string_list_massfractions_sorted.size(); k++)
    if (sorted_max[k] > threshold) current_sorted_index.push_back(k);
}

// 0-based
void ProfilesDatabase::ReactionsAssociatedToSpecies(const unsigned int index,
                                                    std::vector<unsigned int>& indices) {
  kineticsMapXML->stoichiometry().BuildStoichiometricMatrix();

  for (int k = 0; k < kineticsMapXML->stoichiometry().stoichiometric_matrix_reactants().outerSize();
       ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(
             kineticsMapXML->stoichiometry().stoichiometric_matrix_reactants(), k);
         it; ++it)
      if (it.col() == index) indices.push_back(it.row());
  }
  for (int k = 0; k < kineticsMapXML->stoichiometry().stoichiometric_matrix_products().outerSize();
       ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(
             kineticsMapXML->stoichiometry().stoichiometric_matrix_products(), k);
         it; ++it)
      if (it.col() == index) indices.push_back(it.row());
  }

  std::sort(indices.begin(), indices.end());
}

// Function duplicated only cause it calls the kineticsMapSurface instead of the kineticsMap. (Same problem below)
// Map templatization and passing it as argument would remove duplicate code.
// As it is only required for the stoichiometric map, functions are exactly the same always.
void ProfilesDatabase::ReactionsAssociatedToSpecies_Surface(const unsigned int index,
                                                    std::vector<unsigned int>& indices) {
  kineticsMapSurfaceXML->stoichiometry().BuildStoichiometricMatrix();

  for (int k = 0; k < kineticsMapSurfaceXML->stoichiometry().stoichiometric_matrix_reactants().outerSize();
       ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(
             kineticsMapSurfaceXML->stoichiometry().stoichiometric_matrix_reactants(), k);
         it; ++it)
      if (it.col() == index) indices.push_back(it.row());
  }
  for (int k = 0; k < kineticsMapSurfaceXML->stoichiometry().stoichiometric_matrix_products().outerSize();
       ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(
             kineticsMapSurfaceXML->stoichiometry().stoichiometric_matrix_products(), k);
         it; ++it)
      if (it.col() == index) indices.push_back(it.row());
  }

  std::sort(indices.begin(), indices.end());
}

void ProfilesDatabase::isReactantProduct(const unsigned int reaction_index,
                                         double& netStoichiometry) {
  kineticsMapXML->stoichiometry().BuildStoichiometricMatrix();

  std::vector<double> reactants_stoich;
  std::vector<double> products_stoich;
  std::vector<double> reactants_indices;
  std::vector<double> products_indices;
  std::vector<double> duplicate_species_indices;

  Eigen::SparseMatrix<double> reactants =
      kineticsMapXML->stoichiometry().stoichiometric_matrix_reactants();
  Eigen::SparseMatrix<double> products =
      kineticsMapXML->stoichiometry().stoichiometric_matrix_products();

  // TODO
  // Here there is a large room for improvement in terms of computational
  // efficency keep in mind that the loops repeted can be easily condensed
  // into one but I don't have time now
  for (int k = 0; k < reactants.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(reactants, k); it; ++it) {
      if (it.row() == reaction_index) {
        reactants_stoich.push_back(it.value());
        reactants_indices.push_back(it.col());
      }
    }
  }

  for (int k = 0; k < products.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(products, k); it; ++it) {
      if (it.row() == reaction_index) {
        products_stoich.push_back(it.value());
        products_indices.push_back(it.col());
      }
    }
  }

  // Find if some of the species for the selected reactions appears
  // on the both side of the reaction

  // 1. sorting the vectors
  std::sort(reactants_indices.begin(), reactants_indices.end());
  std::sort(products_indices.begin(), products_indices.end());

  // 2. declaring result vector to store the common elements
  std::vector<double> common_species(reactants_indices.size() + products_indices.size());

  // 3. iterator to store return type
  std::vector<double>::iterator it, end;

  end = std::set_intersection(reactants_indices.begin(), reactants_indices.end(),
                              products_indices.begin(), products_indices.end(),
                              common_species.begin());

  for (it = common_species.begin(); it != end; it++) duplicate_species_indices.push_back(*it);

  if (reactants_indices.size() != 1) {
    netStoichiometry = 1;
  } else if (duplicate_species_indices.size() == 1) {
    double idx = duplicate_species_indices[0];
    int pos_r, pos_p;

    std::vector<double>::iterator it_r;
    std::vector<double>::iterator it_p;

    it_r = std::find(reactants_indices.begin(), reactants_indices.end(), idx);
    if (it_r != reactants_indices.end()) pos_r = it_r - reactants_indices.begin();

    it_p = std::find(products_indices.begin(), products_indices.end(), idx);
    if (it_r != products_indices.end()) pos_p = it_p - products_indices.begin();

    netStoichiometry = -reactants_stoich[pos_r] + products_stoich[pos_p];
  } else if (duplicate_species_indices.size() == 0) {
    netStoichiometry = 1;
  } else {
    std::string msg = "Something is wrong with the reaction you are asking for!";
    msg += "Reaction id: " + std::to_string(reaction_index);
    throw std::invalid_argument(msg);
  }
}

void ProfilesDatabase::isReactantProduct_Surface(const unsigned int reaction_index,
                                         double& netStoichiometry) {
  kineticsMapSurfaceXML->stoichiometry().BuildStoichiometricMatrix();

  std::vector<double> reactants_stoich;
  std::vector<double> products_stoich;
  std::vector<double> reactants_indices;
  std::vector<double> products_indices;
  std::vector<double> duplicate_species_indices;

  Eigen::SparseMatrix<double> reactants =
      kineticsMapSurfaceXML->stoichiometry().stoichiometric_matrix_reactants();
  Eigen::SparseMatrix<double> products =
      kineticsMapSurfaceXML->stoichiometry().stoichiometric_matrix_products();

  // TODO
  // Here there is a large room for improvement in terms of computational
  // efficency keep in mind that the loops repeted can be easily condensed
  // into one but I don't have time now
  for (int k = 0; k < reactants.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(reactants, k); it; ++it) {
      if (it.row() == reaction_index) {
        reactants_stoich.push_back(it.value());
        reactants_indices.push_back(it.col());
      }
    }
  }

  for (int k = 0; k < products.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(products, k); it; ++it) {
      if (it.row() == reaction_index) {
        products_stoich.push_back(it.value());
        products_indices.push_back(it.col());
      }
    }
  }

  // Find if some of the species for the selected reactions appears
  // on the both side of the reaction

  // 1. sorting the vectors
  std::sort(reactants_indices.begin(), reactants_indices.end());
  std::sort(products_indices.begin(), products_indices.end());

  // 2. declaring result vector to store the common elements
  std::vector<double> common_species(reactants_indices.size() + products_indices.size());

  // 3. iterator to store return type
  std::vector<double>::iterator it, end;

  end = std::set_intersection(reactants_indices.begin(), reactants_indices.end(),
                              products_indices.begin(), products_indices.end(),
                              common_species.begin());

  for (it = common_species.begin(); it != end; it++) duplicate_species_indices.push_back(*it);

  if (reactants_indices.size() != 1) {
    netStoichiometry = 1;
  } else if (duplicate_species_indices.size() == 1) {
    double idx = duplicate_species_indices[0];
    int pos_r, pos_p;

    std::vector<double>::iterator it_r;
    std::vector<double>::iterator it_p;

    it_r = std::find(reactants_indices.begin(), reactants_indices.end(), idx);
    if (it_r != reactants_indices.end()) pos_r = it_r - reactants_indices.begin();

    it_p = std::find(products_indices.begin(), products_indices.end(), idx);
    if (it_r != products_indices.end()) pos_p = it_p - products_indices.begin();

    netStoichiometry = -reactants_stoich[pos_r] + products_stoich[pos_p];
  } else if (duplicate_species_indices.size() == 0) {
    netStoichiometry = 1;
  } else {
    std::string msg = "Something is wrong with the reaction you are asking for!";
    msg += "Reaction id: " + std::to_string(reaction_index);
    throw std::invalid_argument(msg);
  }
}