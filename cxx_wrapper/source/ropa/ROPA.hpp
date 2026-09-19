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

#include <algorithm>

#include "PostProcessorFluxMap.h"
#include "math/OpenSMOKEUtilities.h"

ROPA::ROPA() {
  ropaType_ = "global";
  kineticFolder_ = "";
  outputFolder_ = "";
  species_ = "";

  localValue_ = 0;
  upperBound_ = 0;
  lowerBound_ = 0;
}

void ROPA::SetKineticFolder(const std::string kineticFolder) { kineticFolder_ = kineticFolder; }

void ROPA::SetOutputFolder(const std::string outputFolder) { outputFolder_ = outputFolder; }

void ROPA::SetROPAType(const std::string ropaType) {
  if (ropaType != "global" && ropaType != "local" && ropaType != "region")
    throw std::invalid_argument("Available ROPA types are: global | local | region");
  ropaType_ = ropaType;
}

void ROPA::SetSpecies(const std::string species) { species_ = species; }

void ROPA::SetLocalValue(double localValue) { localValue_ = localValue; }

void ROPA::SetLowerBound(double lowerBound) { lowerBound_ = lowerBound; }

void ROPA::SetUpperBound(double upperBound) { upperBound_ = upperBound; }

void ROPA::SetResults(PostProcessorCore* data) {
  if (!data->HasKinetics()) {
    throw std::invalid_argument("ROPA::SetResults: no kinetics mechanism loaded on this PostProcessorCore");
  }
  if (!data->HasOutput()) {
    throw std::invalid_argument("ROPA::SetResults: no simulation output loaded on this PostProcessorCore");
  }
  data_ = data;
}

void ROPA::SetElement(const std::string element) { element_ = element; }

void ROPA::SetThickness(const std::string thickness) {
  if (thickness != "absolute" && thickness != "relative")
    throw std::invalid_argument("Available thickness types are: absolute | relative(%)");
  thickness_ = thickness;
}

void ROPA::SetFluxAnalysisType(const std::string type) {
  if (type != "production" && type != "destruction")
    throw std::invalid_argument("Available flux analysis types are: production | destruction");
  flux_type_ = type;
}

void ROPA::SetWidth(const int width) { width_ = width; }

void ROPA::SetDepth(const int depth) { depth_ = depth; }

void ROPA::SetThreshold(const double threshold) { threshold_ = threshold; }

void ROPA::SetThicknessLogScale(bool thicknesslogscale) { thicknesslogscale_ = thicknesslogscale; }

void ROPA::SetLabelType(std::string type) {
  if (type != "absolute" && type != "relative")
    throw std::invalid_argument("Available label types are: absolute | relative(%)");
  label_type_ = type;
}

void ROPA::RateOfProductionAnalysis(const unsigned int number_of_reactions) {
  // Select y variables among the species
  if (std::find(data_->string_list_massfractions_sorted.begin(),
                data_->string_list_massfractions_sorted.end(),
                species_) != data_->string_list_massfractions_sorted.end()) {
    speciesIsSelected = true;
  } else {
    throw std::invalid_argument("Please select one of the available species!");
  }

  unsigned int index_of_species;
  for (unsigned int j = 0; j < data_->thermodynamicsMapXML->NumberOfSpecies(); j++) {
    if (speciesIsSelected == true) {
      if (species_ == data_->string_list_massfractions_sorted[j]) {
        index_of_species = data_->sorted_index[j];
        break;
      }
    }
  }

  OpenSMOKE::OpenSMOKEVectorDouble x(data_->thermodynamicsMapXML->NumberOfSpecies());
  OpenSMOKE::OpenSMOKEVectorDouble omega(data_->thermodynamicsMapXML->NumberOfSpecies());
  OpenSMOKE::OpenSMOKEVectorDouble c(data_->thermodynamicsMapXML->NumberOfSpecies());

  std::vector<int> reaction_indices;
  std::vector<double> reaction_coefficients;

  // Local Analysis
  // TODO: are we sure that we cannot put this as a function and call it multiple times with slightly different arguments?
  // It is basically the same thing twice on ROPA and twice on GetReactionRates (+ Surface, which has the same problem with its implementation)
  if (ropaType_ == "local") {
    unsigned int index = 0;
    for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
      if (data_->additional[0][j] >= localValue_) {
        index = j;
        break;
      }
    }
    // Recovers mass fractions
    for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++)
      omega[k + 1] = data_->omega[k][index];

    // Calculates mole fractions
    double MWmix;
    data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
                                                                  omega.GetHandle());

    // Calculates concentrations
    const double P_Pa = data_->additional[data_->index_P][index];
    const double T = data_->additional[data_->index_T][index];
    const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
    Product(cTot, x, &c);

    // Calculates formations rates
    data_->kineticsMapXML->SetTemperature(T);
    data_->kineticsMapXML->SetPressure(P_Pa);
    data_->thermodynamicsMapXML->SetTemperature(T);
    data_->thermodynamicsMapXML->SetPressure(P_Pa);

    data_->kineticsMapXML->KineticConstants();
    data_->kineticsMapXML->ReactionRates(c.GetHandle());

    // Ropa
    OpenSMOKE::ROPA_Data ropa;
    data_->kineticsMapXML->RateOfProductionAnalysis(ropa);

    MergePositiveAndNegativeBars(ropa.production_reaction_indices[index_of_species],
                                 ropa.destruction_reaction_indices[index_of_species],
                                 ropa.production_coefficients[index_of_species],
                                 ropa.destruction_coefficients[index_of_species], reaction_indices,
                                 reaction_coefficients);
  }  // Global | Region
  else {
    unsigned int index_min = 0;
    unsigned int index_max = data_->number_of_abscissas_ - 1;
    if (ropaType_ == "region") {
      for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
        if (data_->additional[0][j] >= lowerBound_) {
          index_min = j;
          break;
        }
      }
      for (unsigned int j = index_min; j < data_->number_of_abscissas_; j++) {
        if (data_->additional[0][j] >= upperBound_) {
          index_max = j;
          break;
        }
      }
      if (index_min == index_max) {
        if (index_max == data_->number_of_abscissas_ - 1)
          index_min = index_max - 1;
        else
          index_max = index_min + 1;
      }
    }

    const double delta = data_->additional[0][index_max] - data_->additional[0][index_min];

    std::vector<double> global_production_coefficients;
    std::vector<double> global_destruction_coefficients;
    std::vector<unsigned int> global_production_reaction_indices;
    std::vector<unsigned int> global_destruction_reaction_indices;

    for (unsigned int j = index_min; j < index_max - 1; j++) {
      // Recovers mass fractions
      for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++)
        omega[k + 1] = data_->omega[k][j];

      // Calculates mole fractions
      double MWmix;
      data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
                                                                    omega.GetHandle());

      // Calculates concentrations
      const double P_Pa = data_->additional[data_->index_P][j];
      const double T = data_->additional[data_->index_T][j];
      const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
      Product(cTot, x, &c);

      // Calculates formations rates
      data_->kineticsMapXML->SetTemperature(T);
      data_->kineticsMapXML->SetPressure(P_Pa);
      data_->thermodynamicsMapXML->SetTemperature(T);
      data_->thermodynamicsMapXML->SetPressure(P_Pa);

      data_->kineticsMapXML->KineticConstants();
      data_->kineticsMapXML->ReactionRates(c.GetHandle());

      // Ropa
      OpenSMOKE::ROPA_Data ropa;
      data_->kineticsMapXML->RateOfProductionAnalysis(ropa);

      if (ropa.production_coefficients[index_of_species].size() !=
              ropa.production_reaction_indices[index_of_species].size() ||
          ropa.destruction_coefficients[index_of_species].size() !=
              ropa.destruction_reaction_indices[index_of_species].size()) {
        throw std::invalid_argument("SSS");
      }

      if (j == index_min) {
        global_production_coefficients.resize(
            ropa.production_coefficients[index_of_species].size());
        global_destruction_coefficients.resize(
            ropa.destruction_coefficients[index_of_species].size());
        global_production_reaction_indices = ropa.production_reaction_indices[index_of_species];
        global_destruction_reaction_indices = ropa.destruction_reaction_indices[index_of_species];
      }

      const double dt = (data_->additional[0][j + 1] - data_->additional[0][j]) / delta;
      for (unsigned int k = 0; k < ropa.production_coefficients[index_of_species].size(); k++)
        global_production_coefficients[k] += dt * ropa.production_coefficients[index_of_species][k];

      for (unsigned int k = 0; k < ropa.destruction_coefficients[index_of_species].size(); k++)
        global_destruction_coefficients[k] +=
            dt * ropa.destruction_coefficients[index_of_species][k];
    }

    MergePositiveAndNegativeBars(global_production_reaction_indices,
                                 global_destruction_reaction_indices,
                                 global_production_coefficients, global_destruction_coefficients,
                                 reaction_indices, reaction_coefficients);
  }

  coefficients_.resize(std::min<int>(number_of_reactions, reaction_coefficients.size()));
  reactions_.resize(std::min<int>(number_of_reactions, reaction_coefficients.size()));

  for (int i = 0; i < std::min<int>(number_of_reactions, reaction_coefficients.size()); i++) {
    coefficients_[i] = reaction_coefficients[i];
    reactions_[i] = reaction_indices[i];
  }
}

void ROPA::FluxAnalysis() {
  // Select y variables among the species
  if (std::find(data_->string_list_massfractions_sorted.begin(),
                data_->string_list_massfractions_sorted.end(),
                species_) != data_->string_list_massfractions_sorted.end()) {
    speciesIsSelected = true;
  } else {
    throw std::invalid_argument("Please select one of the available species!");
  }
  unsigned int index_of_species;
  for (unsigned int j = 0; j < data_->thermodynamicsMapXML->NumberOfSpecies(); j++) {
    if (speciesIsSelected == true) {
      if (species_ == data_->string_list_massfractions_sorted[j]) {
        index_of_species = data_->sorted_index[j];
        break;
      }
    }
  }

  unsigned int index_element;  // = ui.comboBox_Elements->currentIndex();
  std::vector<std::string> elements_names = data_->thermodynamicsMapXML->elements();
  for (unsigned int k = 0; k < elements_names.size(); k++) {
    if (element_ == elements_names[k]) {
      index_element = k;
      break;
    }
  }
  const double n_elements =
      data_->thermodynamicsMapXML->atomic_composition()(index_of_species, index_element);

  if (n_elements == 0.) {
    throw std::invalid_argument("The selected species does not contain the selected element");
  }

  const int max_depth = depth_;
  const int max_width = width_;
  const double min_threshold_percentage = threshold_;

  // Local Analysis (Flux can be done only in when local ropa is available)
  unsigned int index = 0;
  for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
    if (data_->additional[0][j] >= localValue_) {
      index = j;
      break;
    }
  }
  OpenSMOKE::OpenSMOKEVectorDouble x(data_->thermodynamicsMapXML->NumberOfSpecies());
  OpenSMOKE::OpenSMOKEVectorDouble omega(data_->thermodynamicsMapXML->NumberOfSpecies());
  OpenSMOKE::OpenSMOKEVectorDouble c(data_->thermodynamicsMapXML->NumberOfSpecies());
  OpenSMOKE::OpenSMOKEVectorDouble r(data_->kineticsMapXML->NumberOfReactions());

  // Recovers mass fractions
  for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++) {
    omega[k + 1] = data_->omega[k][index];
  }
  // Calculates mole fractions
  double MWmix;
  data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
                                                                omega.GetHandle());

  // Calculates concentrations
  const double P_Pa = data_->additional[data_->index_P][index];
  const double T = data_->additional[data_->index_T][index];
  const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
  Product(cTot, x, &c);

  // Calculates formations rates
  data_->kineticsMapXML->SetTemperature(T);
  data_->kineticsMapXML->SetPressure(P_Pa);
  data_->thermodynamicsMapXML->SetTemperature(T);
  data_->thermodynamicsMapXML->SetPressure(P_Pa);

  data_->kineticsMapXML->KineticConstants();
  data_->kineticsMapXML->ReactionRates(c.GetHandle());
  data_->kineticsMapXML->GiveMeReactionRates(r.GetHandle());

  // OpenSMOKE::FluxAnalysisMap flux_analysis(*data_->thermodynamicsMapXML, *data_->kineticsMapXML);
  pySMOKEPostProcessor::PostProcessorFluxMap flux_analysis(*data_->thermodynamicsMapXML,
                                                           *data_->kineticsMapXML);

  bool destruction = false;
  bool relativethickness = false;
  bool labelrelative = false;

  if (flux_type_ == "destruction") destruction = true;

  if (thickness_ == "relative") relativethickness = true;

  if (label_type_ == "relative") labelrelative = true;

  flux_analysis.SetDestructionAnalysis(destruction);
  flux_analysis.SetNormalThickness(relativethickness);
  flux_analysis.SetNormalTags(labelrelative);
  flux_analysis.SetLogarithmicThickness(thicknesslogscale_);
  flux_analysis.SetMaxDepth(max_depth);
  flux_analysis.SetMaxWidth(max_width);
  flux_analysis.SetMinPercentageThreshold(min_threshold_percentage);
  flux_analysis.SetAtom(index_element);
  flux_analysis.SetReactionRates(r.Size(), r.GetHandle());

  std::vector<unsigned int> important_indices;
  important_indices.push_back(index_of_species);
  flux_analysis.GloballyAnalyze(important_indices, 0);
  flux_analysis.CalculateThickness();

  flux_analysis.ComputeFluxAnalysis();

  indexFirstName_ = flux_analysis.IndexFirstName;
  indexSecondName_ = flux_analysis.IndexSecondName;
  computedThickness_ = flux_analysis.ComputedThicknessValue;
  computedLabel_ = flux_analysis.ComputedLabelValue;
}

void ROPA::GetReactionRates(std::vector<unsigned int> reaction_indices, const bool sum_rates) {
  unsigned int numberOfReactions = reaction_indices.size();
  // Calculate the reaction rates
  {
    sumOfRates_.resize(data_->number_of_abscissas_);
    reactionRates_.resize(numberOfReactions, std::vector<double>(data_->number_of_abscissas_, 1));

    OpenSMOKE::OpenSMOKEVectorDouble x(data_->thermodynamicsMapXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble omega(data_->thermodynamicsMapXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble c(data_->thermodynamicsMapXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble r(data_->kineticsMapXML->NumberOfReactions());

    for (unsigned int i = 0; i < data_->number_of_abscissas_; i++) {
      // Recovers mass fractions
      for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++)
        omega[k + 1] = data_->omega[k][i];

      // Calculate mole fractions
      double MWmix;
      data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
                                                                    omega.GetHandle());

      // Calculate concentrations
      const double P_Pa = data_->additional[data_->index_P][i];
      const double T = data_->additional[data_->index_T][i];
      const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
      Product(cTot, x, &c);

      // Calculate reaction rates
      data_->kineticsMapXML->SetTemperature(T);
      data_->kineticsMapXML->SetPressure(P_Pa);
      data_->thermodynamicsMapXML->SetTemperature(T);
      data_->thermodynamicsMapXML->SetPressure(P_Pa);

      // data_->kineticsMapXML->KineticConstants();
      data_->kineticsMapXML->ReactionRates(c.GetHandle());
      data_->kineticsMapXML->GiveMeReactionRates(r.GetHandle());

      if (sum_rates) {
        double sum_rate = 0.;
        for (unsigned int k = 0; k < numberOfReactions; k++) {
          double multiplication_factor = 1;
          // data_->isReactantProduct(reaction_indices[k], multiplication_factor);
          const unsigned int j = reaction_indices[k] + 1;
          sum_rate += multiplication_factor * r[j];
        }
        sumOfRates_[i] = sum_rate;
      } else {
        for (unsigned int k = 0; k < numberOfReactions; k++) {
          double multiplication_factor = 1;
          // data_->isReactantProduct(reaction_indices[k], multiplication_factor);
          const unsigned int j = reaction_indices[k] + 1;
          // reactionRates_[k][i] = multiplication_factor * r[j];
          reactionRates_[k][i] = r[j];
        }
      }
    }
  }
}

void ROPA::GetFormationRates(std::string specie, std::string units, std::string type) {
  if (units != "mass" && units != "mole")
    throw std::invalid_argument("Available Formation Rates units are: mole | mass");

  // As for the Reaction Rates keep in mind that AC allow the plot of
  // several species at the same time now let's stay simple one at the time
  // Select y variables among the species
  OpenSMOKE::OpenSMOKEVector<unsigned int> formation_rates_to_plot;
  std::string selected_species = specie;

  {
    unsigned int n_selected_species = 1;
    ChangeDimensions(n_selected_species, &formation_rates_to_plot, true);
    for (unsigned int j = 0; j < n_selected_species; j++) {
      for (unsigned int k = 0; k < data_->string_list_massfractions_sorted.size(); k++) {
        if (selected_species == data_->string_list_massfractions_sorted[k]) {
          formation_rates_to_plot[j + 1] = k;
          break;
        }
      }
    }
  }

  // Calculate the formation rates
  {
    formationRates_.resize(data_->number_of_abscissas_);

    OpenSMOKE::OpenSMOKEVectorDouble P(data_->thermodynamicsMapXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble D(data_->thermodynamicsMapXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble x(data_->thermodynamicsMapXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble omega(data_->thermodynamicsMapXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble c(data_->thermodynamicsMapXML->NumberOfSpecies());

    for (unsigned int i = 0; i < data_->number_of_abscissas_; i++) {
      // Recovers mass fractions
      for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++)
        omega[k + 1] = data_->omega[k][i];

      // Calculates mole fractions
      double MWmix;
      data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
                                                                    omega.GetHandle());

      // Calculates concentrations
      const double P_Pa = data_->additional[data_->index_P][i];
      const double T = data_->additional[data_->index_T][i];
      const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
      Product(cTot, x, &c);

      // Calculates formations rates
      data_->kineticsMapXML->SetTemperature(T);
      data_->kineticsMapXML->SetPressure(P_Pa);
      data_->thermodynamicsMapXML->SetTemperature(T);
      data_->thermodynamicsMapXML->SetPressure(P_Pa);

      data_->kineticsMapXML->KineticConstants();
      data_->kineticsMapXML->ReactionRates(c.GetHandle());
      data_->kineticsMapXML->ProductionAndDestructionRates(P.GetHandle(),
                                                           D.GetHandle());  // kmol/m3/s

      if (type == "characteristic-time") {
        const unsigned k = data_->sorted_index[formation_rates_to_plot[1]] + 1;
        formationRates_[i] = c[k] / (D[k] + 1.e-32);
      } else {
        if (units == "mass") {
          OpenSMOKE::ElementByElementProduct(
              P.Size(), P.GetHandle(), data_->thermodynamicsMapXML->MWs().data(), P.GetHandle());
          OpenSMOKE::ElementByElementProduct(
              D.Size(), D.GetHandle(), data_->thermodynamicsMapXML->MWs().data(), D.GetHandle());
        }

        const unsigned k = data_->sorted_index[formation_rates_to_plot[1]] + 1;
        if (type == "net")
          formationRates_[i] = P[k] - D[k];
        else if (type == "production")
          formationRates_[i] = P[k];
        else if (type == "destruction")
          formationRates_[i] = D[k];
        else
          throw std::invalid_argument(
              "Available Formation Rates types are: net | production | destruction | "
              "characteristic-time");
      }
    }
  }
}

void ROPA::RateOfProductionAnalysis2D(const unsigned int number_of_reactions, double local_x,
                                      double local_z, double region_low_x, double region_up_x,
                                      double region_low_z, double region_up_z) {
  // This function is totally inefficent and by far the worst code I have evere written
  // Select y variables among the species
  if (std::find(data_->string_list_massfractions_sorted.begin(),
                data_->string_list_massfractions_sorted.end(),
                species_) != data_->string_list_massfractions_sorted.end()) {
    speciesIsSelected = true;
  } else {
    throw std::invalid_argument("Please select one of the available species!");
  }

  unsigned int index_of_species;
  for (unsigned int j = 0; j < data_->thermodynamicsMapXML->NumberOfSpecies(); j++) {
    if (speciesIsSelected == true) {
      if (species_ == data_->string_list_massfractions_sorted[j]) {
        index_of_species = data_->sorted_index[j];
        break;
      }
    }
  }

  OpenSMOKE::OpenSMOKEVectorDouble x(data_->thermodynamicsMapXML->NumberOfSpecies());
  OpenSMOKE::OpenSMOKEVectorDouble omega(data_->thermodynamicsMapXML->NumberOfSpecies());
  OpenSMOKE::OpenSMOKEVectorDouble c(data_->thermodynamicsMapXML->NumberOfSpecies());

  std::vector<int> reaction_indices;
  std::vector<double> reaction_coefficients;
  // Local Analysis
  if (ropaType_ == "local") {
    unsigned int index = 0;
    for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
      if (data_->additional[data_->index_x_coord][j] >= local_x &&
          data_->additional[data_->index_z_coord][j] >= local_z) {
        index = j;
        break;
      }
    }
    // Recovers mass fractions
    for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++)
      omega[k + 1] = data_->omega[k][index];

    // Calculates mole fractions
    double MWmix;
    data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
                                                                  omega.GetHandle());

    // Calculates concentrations
    const double P_Pa = data_->additional[data_->index_P][index];
    const double T = data_->additional[data_->index_T][index];
    const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
    Product(cTot, x, &c);

    // Calculates formations rates
    data_->kineticsMapXML->SetTemperature(T);
    data_->kineticsMapXML->SetPressure(P_Pa);
    data_->thermodynamicsMapXML->SetTemperature(T);
    data_->thermodynamicsMapXML->SetPressure(P_Pa);

    data_->kineticsMapXML->KineticConstants();
    data_->kineticsMapXML->ReactionRates(c.GetHandle());

    // Ropa
    OpenSMOKE::ROPA_Data ropa;
    data_->kineticsMapXML->RateOfProductionAnalysis(ropa);

    MergePositiveAndNegativeBars(ropa.production_reaction_indices[index_of_species],
                                 ropa.destruction_reaction_indices[index_of_species],
                                 ropa.production_coefficients[index_of_species],
                                 ropa.destruction_coefficients[index_of_species], reaction_indices,
                                 reaction_coefficients);
  }  // Region
  else if (ropaType_ == "region") {
    std::vector<double> global_production_coefficients;
    std::vector<double> global_destruction_coefficients;
    std::vector<unsigned int> global_production_reaction_indices;
    std::vector<unsigned int> global_destruction_reaction_indices;

    unsigned int index_ll;
    unsigned int index_lr;
    unsigned int index_ul;
    unsigned int index_ur;

    for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
      if (data_->additional[data_->index_x_coord][j] >= region_low_x &&
          data_->additional[data_->index_z_coord][j] >= region_low_z) {
        index_ll = j;
        break;
      }
    }

    for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
      if (data_->additional[data_->index_x_coord][j] >= region_low_x &&
          data_->additional[data_->index_z_coord][j] >= region_up_z) {
        index_ul = j;
        break;
      }
    }

    for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
      if (data_->additional[data_->index_x_coord][j] >= region_up_x &&
          data_->additional[data_->index_z_coord][j] >= region_low_z) {
        index_lr = j;
        break;
      }
    }

    for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
      if (data_->additional[data_->index_x_coord][j] >= region_up_x &&
          data_->additional[data_->index_z_coord][j] >= region_up_z) {
        index_ur = j;
        break;
      }
    }
    // TODO ADD CHECK on idices
    // if (index_min == index_max)
    // {
    //     if (index_max == data_->number_of_abscissas_ - 1)
    //         index_min = index_max - 1;
    //     else
    //         index_max = index_min + 1;
    // }
    //

    std::vector<double> volume = data_->additional[data_->index_volume];
    double total_volume;
    for (unsigned int i = 0; i < data_->number_of_abscissas_; i++) {
      if (data_->additional[data_->index_x_coord][i] >=
              data_->additional[data_->index_x_coord][index_ll] &&
          data_->additional[data_->index_x_coord][i] <=
              data_->additional[data_->index_x_coord][index_lr]) {
        if (data_->additional[data_->index_z_coord][i] >=
                data_->additional[data_->index_z_coord][index_ll] &&
            data_->additional[data_->index_z_coord][i] <=
                data_->additional[data_->index_z_coord][index_ul]) {
          total_volume += volume[i];
        }
      }
    }

    int counter = 0;

    for (unsigned int i = 0; i < data_->number_of_abscissas_; i++) {
      if (data_->additional[data_->index_x_coord][i] >=
              data_->additional[data_->index_x_coord][index_ll] &&
          data_->additional[data_->index_x_coord][i] <=
              data_->additional[data_->index_x_coord][index_lr]) {
        if (data_->additional[data_->index_z_coord][i] >=
                data_->additional[data_->index_z_coord][index_ll] &&
            data_->additional[data_->index_z_coord][i] <=
                data_->additional[data_->index_z_coord][index_ul]) {
          // Recovers mass fractions
          for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++)
            omega[k + 1] = data_->omega[k][i];

          // Calculates mole fractions
          double MWmix;
          data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
                                                                        omega.GetHandle());

          // Calculates concentrations
          const double P_Pa = data_->additional[data_->index_P][i];
          const double T = data_->additional[data_->index_T][i];
          const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
          Product(cTot, x, &c);

          // Calculates formations rates
          data_->kineticsMapXML->SetTemperature(T);
          data_->kineticsMapXML->SetPressure(P_Pa);
          data_->thermodynamicsMapXML->SetTemperature(T);
          data_->thermodynamicsMapXML->SetPressure(P_Pa);

          data_->kineticsMapXML->KineticConstants();
          data_->kineticsMapXML->ReactionRates(c.GetHandle());

          // Ropa
          OpenSMOKE::ROPA_Data ropa;
          data_->kineticsMapXML->RateOfProductionAnalysis(ropa);

          if (ropa.production_coefficients[index_of_species].size() !=
                  ropa.production_reaction_indices[index_of_species].size() ||
              ropa.destruction_coefficients[index_of_species].size() !=
                  ropa.destruction_reaction_indices[index_of_species].size()) {
            throw std::invalid_argument("SSS");
          }

          if (counter == 0) {
            global_production_coefficients.resize(
                ropa.production_coefficients[index_of_species].size());
            global_destruction_coefficients.resize(
                ropa.destruction_coefficients[index_of_species].size());
            global_production_reaction_indices = ropa.production_reaction_indices[index_of_species];
            global_destruction_reaction_indices =
                ropa.destruction_reaction_indices[index_of_species];
          }

          double integral_contribution = volume[i] / total_volume;

          for (unsigned int k = 0; k < ropa.production_coefficients[index_of_species].size(); k++)
            global_production_coefficients[k] +=
                integral_contribution * ropa.production_coefficients[index_of_species][k];

          for (unsigned int k = 0; k < ropa.destruction_coefficients[index_of_species].size(); k++)
            global_destruction_coefficients[k] +=
                integral_contribution * ropa.destruction_coefficients[index_of_species][k];

          counter += 1;
        }
      }
    }
    MergePositiveAndNegativeBars(global_production_reaction_indices,
                                 global_destruction_reaction_indices,
                                 global_production_coefficients, global_destruction_coefficients,
                                 reaction_indices, reaction_coefficients);
  } else  // Global
  {
    double tmp = 0;
    unsigned int index_min = 0;
    unsigned int index_max = data_->number_of_abscissas_;

    std::vector<double> volume = data_->additional[data_->index_volume];
    // double total_volume = std::accumulate(volume.begin(), volume.end(), 0.0d);
    double total_volume = 0;
    for (unsigned int i = 0; i < volume.size(); i++) {
      total_volume += volume[i];
    }

    std::vector<double> global_production_coefficients;
    std::vector<double> global_destruction_coefficients;
    std::vector<unsigned int> global_production_reaction_indices;
    std::vector<unsigned int> global_destruction_reaction_indices;

    for (unsigned int j = index_min; j < index_max; j++) {
      // Recovers mass fractions
      for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++)
        omega[k + 1] = data_->omega[k][j];

      // Calculates mole fractions
      double MWmix;
      data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
                                                                    omega.GetHandle());

      // Calculates concentrations
      const double P_Pa = data_->additional[data_->index_P][j];
      const double T = data_->additional[data_->index_T][j];
      const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
      Product(cTot, x, &c);

      // Calculates formations rates
      data_->kineticsMapXML->SetTemperature(T);
      data_->kineticsMapXML->SetPressure(P_Pa);
      data_->thermodynamicsMapXML->SetTemperature(T);
      data_->thermodynamicsMapXML->SetPressure(P_Pa);

      data_->kineticsMapXML->KineticConstants();
      data_->kineticsMapXML->ReactionRates(c.GetHandle());

      // Ropa
      OpenSMOKE::ROPA_Data ropa;
      data_->kineticsMapXML->RateOfProductionAnalysis(ropa);

      if (ropa.production_coefficients[index_of_species].size() !=
              ropa.production_reaction_indices[index_of_species].size() ||
          ropa.destruction_coefficients[index_of_species].size() !=
              ropa.destruction_reaction_indices[index_of_species].size()) {
        throw std::invalid_argument("SSS");
      }

      if (j == index_min) {
        global_production_coefficients.resize(
            ropa.production_coefficients[index_of_species].size());
        global_destruction_coefficients.resize(
            ropa.destruction_coefficients[index_of_species].size());
        global_production_reaction_indices = ropa.production_reaction_indices[index_of_species];
        global_destruction_reaction_indices = ropa.destruction_reaction_indices[index_of_species];
      }

      double integral_contribution = volume[j] / total_volume;

      for (unsigned int k = 0; k < ropa.production_coefficients[index_of_species].size(); k++)
        global_production_coefficients[k] +=
            integral_contribution * ropa.production_coefficients[index_of_species][k];

      for (unsigned int k = 0; k < ropa.destruction_coefficients[index_of_species].size(); k++)
        global_destruction_coefficients[k] +=
            integral_contribution * ropa.destruction_coefficients[index_of_species][k];
    }

    MergePositiveAndNegativeBars(global_production_reaction_indices,
                                 global_destruction_reaction_indices,
                                 global_production_coefficients, global_destruction_coefficients,
                                 reaction_indices, reaction_coefficients);
  }

  coefficients_.resize(std::min<int>(number_of_reactions, reaction_coefficients.size()));
  reactions_.resize(std::min<int>(number_of_reactions, reaction_coefficients.size()));

  for (int i = 0; i < std::min<int>(number_of_reactions, reaction_coefficients.size()); i++) {
    coefficients_[i] = reaction_coefficients[i];
    reactions_[i] = reaction_indices[i];
  }
}
