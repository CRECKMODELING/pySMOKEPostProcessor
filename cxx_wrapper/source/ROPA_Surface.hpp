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

#include "PostProcessorFluxMap.h"  // For gas-phase only
#include "math/OpenSMOKEUtilities.h"

ROPA_Surface::ROPA_Surface() {
  ROPA();  // Fallback on parent constructor
  heterogeneous_reactions_ = false;
}

void ROPA_Surface::SetROPAPhase(const bool heterogeneous_reactions) {
  heterogeneous_reactions_ = heterogeneous_reactions;
}

void ROPA_Surface::RateOfProductionAnalysis(const unsigned int number_of_reactions,
                                            const bool heterogeneous_reactions) {
  heterogeneous_reactions_ =
      heterogeneous_reactions;  // Update internal variable (maybe useless)
  if (heterogeneous_reactions_ == false) {
    // It might be correct to throw an exception here if het_rxns == false and the
    // selected species is from heterogeneous phase. Case homogeneous reactions: fallback
    // on homogeneous ROPA functions
    ROPA::RateOfProductionAnalysis(number_of_reactions);
  } else {
    // Case heterogeneous reactions: dedicated implementation
    if (std::find(data_->string_list_massfractions_sorted.begin(),
                  data_->string_list_massfractions_sorted.end(),
                  species_) != data_->string_list_massfractions_sorted.end()) {
      speciesIsSelected = true;
    } else if (std::find(data_->string_list_surfacefractions_sorted.begin(),
                         data_->string_list_surfacefractions_sorted.end(),
                         species_) != data_->string_list_surfacefractions_sorted.end()) {
      speciesIsSelected = true;
    } else if (std::find(data_->string_list_bulkmasses_sorted.begin(),
                         data_->string_list_bulkmasses_sorted.end(),
                         species_) != data_->string_list_bulkmasses_sorted.end()) {
      speciesIsSelected = true;
    } else {
      throw std::invalid_argument("Please select one of the available species!");
    }

    const unsigned int NSG = data_->thermodynamicsMapSurfaceXML->number_of_gas_species();
    const unsigned int NSS = data_->thermodynamicsMapSurfaceXML->number_of_site_species();
    const unsigned int NSB = data_->thermodynamicsMapSurfaceXML->number_of_bulk_species();

    unsigned int index_of_species;
    bool species_not_found = true;

    {  // Species search (it is probably not particularly efficient...)
      unsigned int dummy_index = 0;
      for (unsigned int j = 0; j < NSG; j++) {
        if (species_ == data_->string_list_massfractions_sorted[j]) {
          dummy_index = data_->sorted_index[j];
          index_of_species = dummy_index;
          species_not_found = false;
          break;
        }
      }

      if (species_not_found == true) {
        for (unsigned int j = 0; j < NSS; j++) {
          if (species_ == data_->string_list_surfacefractions_sorted[j]) {
            dummy_index = data_->sorted_index_surface[j];
            index_of_species = NSG + dummy_index;
            species_not_found = false;
            break;
          }
        }
      }
      if (species_not_found == true) {
        for (unsigned int j = 0; j < NSB; j++) {
          if (species_ == data_->string_list_bulkmasses_sorted[j]) {
            dummy_index = data_->sorted_index_bulk[j];
            index_of_species = NSG + NSS + dummy_index;
            species_not_found = false;
            break;
          }
        }
      }

      if (species_not_found == true) {
        throw std::invalid_argument("Species not found in any phase!");
      }
    }

    OpenSMOKE::OpenSMOKEVectorDouble x(NSG);
    OpenSMOKE::OpenSMOKEVectorDouble omega(NSG);
    OpenSMOKE::OpenSMOKEVectorDouble cGas(NSG);
    OpenSMOKE::OpenSMOKEVectorDouble Z(NSS);
    OpenSMOKE::OpenSMOKEVectorDouble cSurf(NSS);
    OpenSMOKE::OpenSMOKEVectorDouble aBulk(NSB);
    // Important note: in the XML file, and in the ProfilesDatabase file, the bulk
    // variable is massBulk which is the actual evolving variable For kinetics
    // evaluations, however, we consider the solid activity to be equal to 1, which is
    // also the parameter I'm passing to the ROPA function.

    std::vector<int> reaction_indices;
    std::vector<double> reaction_coefficients;

    if (ropaType_ == "local") {
      unsigned int index = 0;
      for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
        if (data_->additional[0][j] >= localValue_) {
          index = j;
          break;
        }
      }
      const double Gamma =
          data_->additional[data_->index_surface_sites_concentration]
                           [index];  // In OpenSMOKE, this is a OpenSMOKEVectorDouble
                                     // (more general for multiple surface phases).

      unsigned int k = 1;
      for (unsigned int j = 0; j < NSG; j++) {
        omega[j + 1] = data_->omega[j][index];
      }

      for (unsigned int j = 0; j < NSS; j++) {
        Z[j + 1] = data_->Z[j][index];
      }
      for (unsigned int j = 0; j < NSB; j++) {
        aBulk[j + 1] = 1.;  // If, in the future, a more complex model accounting for
                            // solid activity is introduced, this will have to be changed
      }

      // Calculates mole fractions
      double MWmix;
      data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
                                                                    omega.GetHandle());

      // Calculates gas-phase concentrations
      const double P_Pa = data_->additional[data_->index_P][index];
      const double T = data_->additional[data_->index_T][index];
      const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
      Product(cTot, x, &cGas);

      OpenSMOKE::ROPA_Data ropa;

      // Calculates formations rates
      data_->kineticsMapSurfaceXML->SetTemperature(T);
      data_->kineticsMapSurfaceXML->SetPressure(P_Pa);
      // Why are we setting T,P in thermo map? It does not seem to be used? 
      // (In standard ROPA it is the same)
      data_->thermodynamicsMapSurfaceXML->SetTemperature(T);
      data_->thermodynamicsMapSurfaceXML->SetPressure(P_Pa);

      data_->kineticsMapSurfaceXML->KineticConstants();
      data_->kineticsMapSurfaceXML->ReactionRates(cGas.GetHandle(), Z.GetHandle(),
                                                  aBulk.GetHandle(), &Gamma);
      // Performs ROPA
      data_->kineticsMapSurfaceXML->RateOfProductionAnalysis(ropa);

      ROPA::MergePositiveAndNegativeBars(
          ropa.production_reaction_indices[index_of_species],
          ropa.destruction_reaction_indices[index_of_species],
          ropa.production_coefficients[index_of_species],
          ropa.destruction_coefficients[index_of_species], reaction_indices,
          reaction_coefficients);

    } else {
      // Global/Region ROPA
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

      const double delta =
          data_->additional[0][index_max] - data_->additional[0][index_min];

      std::vector<double> global_production_coefficients;
      std::vector<double> global_destruction_coefficients;
      std::vector<unsigned int> global_production_reaction_indices;
      std::vector<unsigned int> global_destruction_reaction_indices;

      for (unsigned int idx = index_min; idx < index_max - 1; idx++) {
        const double Gamma =
            data_->additional[data_->index_surface_sites_concentration][idx];
        // In OpenSMOKE, this is a OpenSMOKEVectorDouble (more general for multiple
        // phases).

        for (unsigned int j = 0; j < NSG; j++) omega[j + 1] = data_->omega[j][idx];
        for (unsigned int j = 0; j < NSS; j++) Z[j + 1] = data_->Z[j][idx];
        for (unsigned int j = 0; j < NSB; j++) aBulk[j + 1] = 1.;
        // If, in the future, a more complex model accounting for solid activity is
        // introduced, this will have to be changed

        // Calculates mole fractions
        double MWmix;
        data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(
            x.GetHandle(), MWmix, omega.GetHandle());
        // Calculates concentrations
        const double P_Pa = data_->additional[data_->index_P][idx];
        const double T = data_->additional[data_->index_T][idx];
        const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
        Product(cTot, x, &cGas);

        OpenSMOKE::ROPA_Data ropa;
        // Calculates formations rates heterogeneous
        data_->kineticsMapSurfaceXML->SetTemperature(T);
        data_->kineticsMapSurfaceXML->SetPressure(P_Pa);
        // Why are we setting T,P in thermo map? It does not seem to be used? (In standard
        // ROPA it is the same)
        data_->thermodynamicsMapSurfaceXML->SetTemperature(T);
        data_->thermodynamicsMapSurfaceXML->SetPressure(P_Pa);

        data_->kineticsMapSurfaceXML->KineticConstants();
        data_->kineticsMapSurfaceXML->ReactionRates(cGas.GetHandle(), Z.GetHandle(),
                                                    aBulk.GetHandle(), &Gamma);
        // Performs ROPA
        data_->kineticsMapSurfaceXML->RateOfProductionAnalysis(ropa);

        if (ropa.production_coefficients[index_of_species].size() !=
                ropa.production_reaction_indices[index_of_species].size() ||
            ropa.destruction_coefficients[index_of_species].size() !=
                ropa.destruction_reaction_indices[index_of_species].size()) {
          throw std::invalid_argument("SSS");
        }

        if (idx == index_min) {
          global_production_coefficients.resize(
              ropa.production_coefficients[index_of_species].size());
          global_destruction_coefficients.resize(
              ropa.destruction_coefficients[index_of_species].size());
          global_production_reaction_indices =
              ropa.production_reaction_indices[index_of_species];
          global_destruction_reaction_indices =
              ropa.destruction_reaction_indices[index_of_species];
        }

        const double dt =
            (data_->additional[0][idx + 1] - data_->additional[0][idx]) / delta;

        for (unsigned int k = 0;
             k < ropa.production_coefficients[index_of_species].size(); k++)
          global_production_coefficients[k] +=
              dt * ropa.production_coefficients[index_of_species][k];

        for (unsigned int k = 0;
             k < ropa.destruction_coefficients[index_of_species].size(); k++)
          global_destruction_coefficients[k] +=
              dt * ropa.destruction_coefficients[index_of_species][k];
      }
      ROPA::MergePositiveAndNegativeBars(
          global_production_reaction_indices, global_destruction_reaction_indices,
          global_production_coefficients, global_destruction_coefficients,
          reaction_indices, reaction_coefficients);
    }

    coefficients_.resize(
        std::min<int>(number_of_reactions, reaction_coefficients.size()));
    reactions_.resize(std::min<int>(number_of_reactions, reaction_coefficients.size()));

    for (int i = 0; i < std::min<int>(number_of_reactions, reaction_coefficients.size());
         i++) {
      coefficients_[i] = reaction_coefficients[i];
      reactions_[i] = reaction_indices[i];
    }
  }
}

// void ROPA_Surface::FluxAnalysis() {
//   // Select y variables among the species
//   if (std::find(data_->string_list_massfractions_sorted.begin(),
//                 data_->string_list_massfractions_sorted.end(),
//                 species_) != data_->string_list_massfractions_sorted.end()) {
//     speciesIsSelected = true;
//   } else {
//     throw std::invalid_argument("Please select one of the available species!");
//   }
//   unsigned int index_of_species;
//   for (unsigned int j = 0; j < data_->thermodynamicsMapXML->NumberOfSpecies(); j++) {
//     if (speciesIsSelected == true) {
//       if (species_ == data_->string_list_massfractions_sorted[j]) {
//         index_of_species = data_->sorted_index[j];
//         break;
//       }
//     }
//   }

//   unsigned int index_element;  // = ui.comboBox_Elements->currentIndex();
//   std::vector<std::string> elements_names = data_->thermodynamicsMapXML->elements();
//   for (unsigned int k = 0; k < elements_names.size(); k++) {
//     if (element_ == elements_names[k]) {
//       index_element = k;
//       break;
//     }
//   }
//   const double n_elements =
//       data_->thermodynamicsMapXML->atomic_composition()(index_of_species,
//       index_element);

//   if (n_elements == 0.) {
//     throw std::invalid_argument("The selected species does not contain the selected
//     element");
//   }

//   const int max_depth = depth_;
//   const int max_width = width_;
//   const double min_threshold_percentage = threshold_;

//   // Local Analysis (Flux can be done only in when local ropa is available)
//   unsigned int index = 0;
//   for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
//     if (data_->additional[0][j] >= localValue_) {
//       index = j;
//       break;
//     }
//   }
//   OpenSMOKE::OpenSMOKEVectorDouble x(data_->thermodynamicsMapXML->NumberOfSpecies());
//   OpenSMOKE::OpenSMOKEVectorDouble
//   omega(data_->thermodynamicsMapXML->NumberOfSpecies());
//   OpenSMOKE::OpenSMOKEVectorDouble c(data_->thermodynamicsMapXML->NumberOfSpecies());
//   OpenSMOKE::OpenSMOKEVectorDouble r(data_->kineticsMapXML->NumberOfReactions());

//   // Recovers mass fractions
//   for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++) {
//     omega[k + 1] = data_->omega[k][index];
//   }
//   // Calculates mole fractions
//   double MWmix;
//   data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
//                                                                 omega.GetHandle());

//   // Calculates concentrations
//   const double P_Pa = data_->additional[data_->index_P][index];
//   const double T = data_->additional[data_->index_T][index];
//   const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
//   Product(cTot, x, &c);

//   // Calculates formations rates
//   data_->kineticsMapXML->SetTemperature(T);
//   data_->kineticsMapXML->SetPressure(P_Pa);
//   data_->thermodynamicsMapXML->SetTemperature(T);
//   data_->thermodynamicsMapXML->SetPressure(P_Pa);

//   data_->kineticsMapXML->KineticConstants();
//   data_->kineticsMapXML->ReactionRates(c.GetHandle());
//   data_->kineticsMapXML->GiveMeReactionRates(r.GetHandle());

//   // OpenSMOKE::FluxAnalysisMap flux_analysis(*data_->thermodynamicsMapXML,
//   *data_->kineticsMapXML); pySMOKEPostProcessor::PostProcessorFluxMap
//   flux_analysis(*data_->thermodynamicsMapXML,
//                                                            *data_->kineticsMapXML);

//   bool destruction = false;
//   bool relativethickness = false;
//   bool labelrelative = false;

//   if (flux_type_ == "destruction") destruction = true;

//   if (thickness_ == "relative") relativethickness = true;

//   if (label_type_ == "relative") labelrelative = true;

//   flux_analysis.SetDestructionAnalysis(destruction);
//   flux_analysis.SetNormalThickness(relativethickness);
//   flux_analysis.SetNormalTags(labelrelative);
//   flux_analysis.SetLogarithmicThickness(thicknesslogscale_);
//   flux_analysis.SetMaxDepth(max_depth);
//   flux_analysis.SetMaxWidth(max_width);
//   flux_analysis.SetMinPercentageThreshold(min_threshold_percentage);
//   flux_analysis.SetAtom(index_element);
//   flux_analysis.SetReactionRates(r.Size(), r.GetHandle());

//   std::vector<unsigned int> important_indices;
//   important_indices.push_back(index_of_species);
//   flux_analysis.GloballyAnalyze(important_indices, 0);
//   flux_analysis.CalculateThickness();

//   flux_analysis.ComputeFluxAnalysis();

//   indexFirstName_ = flux_analysis.IndexFirstName;
//   indexSecondName_ = flux_analysis.IndexSecondName;
//   computedThickness_ = flux_analysis.ComputedThicknessValue;
//   computedLabel_ = flux_analysis.ComputedLabelValue;
// }

void ROPA_Surface::GetReactionRates(std::vector<unsigned int> reaction_indices,
                                    const bool sum_rates,
                                    const bool heterogeneous_reactions) {
  heterogeneous_reactions_ = heterogeneous_reactions;
  if (heterogeneous_reactions_ ==
      false) {  // For homogeneous reactions, fallback on standard function
    ROPA::GetReactionRates(reaction_indices, sum_rates);
  } else {  // Dedicated implementation
    const unsigned int NR = data_->kineticsMapSurfaceXML->NumberOfReactions();
    // Calculate the reaction rates
    {
      sumOfRates_.resize(data_->number_of_abscissas_);
      reactionRates_.resize(NR, std::vector<double>(data_->number_of_abscissas_, 1));

      const unsigned int NSG =
          data_->thermodynamicsMapSurfaceXML->number_of_gas_species();
      const unsigned int NSS =
          data_->thermodynamicsMapSurfaceXML->number_of_site_species();
      const unsigned int NSB =
          data_->thermodynamicsMapSurfaceXML->number_of_bulk_species();

      OpenSMOKE::OpenSMOKEVectorDouble x(NSG);
      OpenSMOKE::OpenSMOKEVectorDouble omega(NSG);
      OpenSMOKE::OpenSMOKEVectorDouble cGas(NSG);

      OpenSMOKE::OpenSMOKEVectorDouble Z(NSS);
      OpenSMOKE::OpenSMOKEVectorDouble cSurf(NSS);
      OpenSMOKE::OpenSMOKEVectorDouble aBulk(NSB);

      OpenSMOKE::OpenSMOKEVectorDouble r(NR);

      for (unsigned int idx = 0; idx < data_->number_of_abscissas_; idx++) {
        const double Gamma =
            data_->additional[data_->index_surface_sites_concentration][idx];
        // In OpenSMOKE, this is a OpenSMOKEVectorDouble (more general for multiple
        // phases).

        for (unsigned int j = 0; j < NSG; j++) omega[j + 1] = data_->omega[j][idx];
        for (unsigned int j = 0; j < NSS; j++) Z[j + 1] = data_->Z[j][idx];
        for (unsigned int j = 0; j < NSB; j++) aBulk[j + 1] = 1.;
        // If, in the future, a more complex model accounting for solid activity is
        // introduced, this will have to be changed

        // Calculates mole fractions
        double MWmix;
        data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(
            x.GetHandle(), MWmix, omega.GetHandle());
        // Calculates concentrations
        const double P_Pa = data_->additional[data_->index_P][idx];
        const double T = data_->additional[data_->index_T][idx];
        const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
        Product(cTot, x, &cGas);

        OpenSMOKE::ROPA_Data ropa;
        // Calculates formations rates heterogeneous
        data_->kineticsMapSurfaceXML->SetTemperature(T);
        data_->kineticsMapSurfaceXML->SetPressure(P_Pa);
        data_->thermodynamicsMapSurfaceXML->SetTemperature(T);
        data_->thermodynamicsMapSurfaceXML->SetPressure(P_Pa);

        data_->kineticsMapSurfaceXML->ReactionRates(cGas.GetHandle(), Z.GetHandle(),
                                                    aBulk.GetHandle(), &Gamma);
        data_->kineticsMapSurfaceXML->GiveMeReactionRates(r.GetHandle());

        if (sum_rates) {
          double sum_rate = 0.;
          for (unsigned int k = 0; k < NR; k++) {
            double multiplication_factor = 1;
            // data_->isReactantProduct(reaction_indices[k], multiplication_factor);
            const unsigned int j = reaction_indices[k] + 1;
            sum_rate += multiplication_factor * r[j];
          }
          sumOfRates_[idx] = sum_rate;
        } else {
          for (unsigned int k = 0; k < NR; k++) {
            double multiplication_factor = 1;
            // data_->isReactantProduct(reaction_indices[k], multiplication_factor);
            const unsigned int j = reaction_indices[k] + 1;
            // reactionRates_[k][i] = multiplication_factor * r[j];
            reactionRates_[k][idx] = r[j];
          }
        }
      }
    }
  }
}

void ROPA_Surface::GetFormationRates(std::string specie, std::string units,
                                     std::string type,
                                     const bool heterogeneous_reactions) {
  heterogeneous_reactions_ = heterogeneous_reactions;
  if (heterogeneous_reactions_ ==
      false) {  // For homogeneous reactions, fallback on standard function
    ROPA::GetFormationRates(specie, units, type);
  } else {  // Dedicated implementation

    if (units != "mass" && units != "mole")
      throw std::invalid_argument("Available Formation Rates units are: mole | mass");

    // As for the Reaction Rates keep in mind that AC allow the plot of
    // several species at the same time now let's stay simple one at the time
    // Select y variables among the species
    OpenSMOKE::OpenSMOKEVector<unsigned int> formation_rates_to_plot;
    std::string selected_species = specie;

    unsigned int index_of_species;
    bool species_not_found = true;
    unsigned int phase_of_species = 7;  // 0 for gas, 1 for surface, 2 for bulk

    const unsigned int NSG = data_->thermodynamicsMapSurfaceXML->number_of_gas_species();
    const unsigned int NSS = data_->thermodynamicsMapSurfaceXML->number_of_site_species();
    const unsigned int NSB = data_->thermodynamicsMapSurfaceXML->number_of_bulk_species();

    {
      for (unsigned int j = 0; j < NSG; j++) {
        if (species_ == data_->string_list_massfractions_sorted[j]) {
          species_not_found = false;
          phase_of_species = 0;
          break;
        }
      }

      if (species_not_found == true) {
        for (unsigned int j = 0; j < NSS; j++) {
          if (species_ == data_->string_list_surfacefractions_sorted[j]) {
            species_not_found = false;
            phase_of_species = 1;
            break;
          }
        }
      }
      if (species_not_found == true) {
        for (unsigned int j = 0; j < NSB; j++) {
          if (species_ == data_->string_list_bulkmasses_sorted[j]) {
            species_not_found = false;
            phase_of_species = 2;
            break;
          }
        }
      }

      if (species_not_found == true) {
        throw std::invalid_argument("Species not found in any phase!");
      }

      unsigned int n_selected_species = 1;
      ChangeDimensions(n_selected_species, &formation_rates_to_plot, true);

      if (phase_of_species == 0) {
        for (unsigned int j = 0; j < n_selected_species; j++) {
          for (unsigned int k = 0; k < data_->string_list_massfractions_sorted.size();
               k++) {
            if (selected_species == data_->string_list_massfractions_sorted[k]) {
              formation_rates_to_plot[j + 1] = k;
              break;
            }
          }
        }
      } else if (phase_of_species == 1) {
        for (unsigned int j = 0; j < n_selected_species; j++) {
          for (unsigned int k = 0; k < data_->string_list_surfacefractions_sorted.size();
               k++) {
            if (selected_species == data_->string_list_surfacefractions_sorted[k]) {
              formation_rates_to_plot[j + 1] = k + NSG;
              break;
            }
          }
        }
      } else {
        for (unsigned int j = 0; j < n_selected_species; j++) {
          for (unsigned int k = 0; k < data_->string_list_bulkmasses_sorted.size(); k++) {
            if (selected_species == data_->string_list_bulkmasses_sorted[k]) {
              formation_rates_to_plot[j + 1] = k + NSG + NSS;
              break;
            }
          }
        }
      }  // base case (exception) with phase_of_species should not be a problem 
         // as it is caught by previous checks
    }
    OpenSMOKE::OpenSMOKEVectorDouble P(NSG + NSS + NSB);
    OpenSMOKE::OpenSMOKEVectorDouble D(NSG + NSS + NSB);

    OpenSMOKE::OpenSMOKEVectorDouble x(NSG);
    OpenSMOKE::OpenSMOKEVectorDouble omega(NSG);
    OpenSMOKE::OpenSMOKEVectorDouble cGas(NSG);
    OpenSMOKE::OpenSMOKEVectorDouble Z(NSS);
    OpenSMOKE::OpenSMOKEVectorDouble cSurf(NSS);
    OpenSMOKE::OpenSMOKEVectorDouble aBulk(NSB);

    // Calculate the formation rates
    {
      formationRates_.resize(data_->number_of_abscissas_);

      for (unsigned int index = 0; index < data_->number_of_abscissas_; index++) {
        for (unsigned int j = 0; j < NSG; j++) omega[j + 1] = data_->omega[j][index];

        for (unsigned int j = 0; j < NSS; j++) Z[j + 1] = data_->Z[j][index];

        for (unsigned int j = 0; j < NSB; j++)
          aBulk[j + 1] = 1.;  
        // If, in the future, a more complex model accounting for solid
        // activity is introduced, this will have to be changed

        const double Gamma =
            data_->additional[data_->index_surface_sites_concentration][index];
        // Calculates mole fractions
        double MWmix;
        data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(
            x.GetHandle(), MWmix, omega.GetHandle());

        // Calculates gas-phase concentrations
        const double P_Pa = data_->additional[data_->index_P][index];
        const double T = data_->additional[data_->index_T][index];
        const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
        Product(cTot, x, &cGas);

        // Calculates formations rates
        data_->kineticsMapSurfaceXML->SetTemperature(T);
        data_->kineticsMapSurfaceXML->SetPressure(P_Pa);

        data_->thermodynamicsMapSurfaceXML->SetTemperature(T);
        data_->thermodynamicsMapSurfaceXML->SetPressure(P_Pa);

        data_->kineticsMapSurfaceXML->KineticConstants();
        data_->kineticsMapSurfaceXML->ReactionRates(cGas.GetHandle(), Z.GetHandle(),
                                                    aBulk.GetHandle(), &Gamma);
        data_->kineticsMapSurfaceXML->ProductionAndDestructionRates(P.GetHandle(),
                                                                    D.GetHandle());

        // LG: disabling characteristic times for now, maybe re-add later (cGas, Z,
        // aBulk..., tau(Surf) = Z*Gamma/D I suppose?) if (type == "characteristic-time")
        // {
        //   const unsigned k = data_->sorted_index[formation_rates_to_plot[1]] + 1;
        //   formationRates_[index] = c[k] / (D[k] + 1.e-32);
        // } else {
        if (units == "mass") {
          OpenSMOKE::ElementByElementProduct(
              P.Size(), P.GetHandle(), data_->thermodynamicsMapSurfaceXML->MWs().data(),
              P.GetHandle());
          OpenSMOKE::ElementByElementProduct(
              D.Size(), D.GetHandle(), data_->thermodynamicsMapSurfaceXML->MWs().data(),
              D.GetHandle());
        }

        const unsigned k = data_->sorted_index[formation_rates_to_plot[1]] + 1;
        if (type == "net")
          formationRates_[index] = P[k] - D[k];
        else if (type == "production")
          formationRates_[index] = P[k];
        else if (type == "destruction")
          formationRates_[index] = D[k];
        else
          throw std::invalid_argument(
              "Available Heterogeneous Formation Rates types are: net | production | "
              "destruction");  // | "
                               //"characteristic-time");
        //} // End of characteristic time if/else
      }
    }
  }
}

// void ROPA_Surface::RateOfProductionAnalysis2D(const unsigned int number_of_reactions,
// double local_x,
//                                       double local_z, double region_low_x, double
//                                       region_up_x, double region_low_z, double
//                                       region_up_z) {
//   // This function is totally inefficent and by far the worst code I have evere written
//   // Select y variables among the species
//   if (std::find(data_->string_list_massfractions_sorted.begin(),
//                 data_->string_list_massfractions_sorted.end(),
//                 species_) != data_->string_list_massfractions_sorted.end()) {
//     speciesIsSelected = true;
//   } else {
//     throw std::invalid_argument("Please select one of the available species!");
//   }

//   unsigned int index_of_species;
//   for (unsigned int j = 0; j < data_->thermodynamicsMapXML->NumberOfSpecies(); j++) {
//     if (speciesIsSelected == true) {
//       if (species_ == data_->string_list_massfractions_sorted[j]) {
//         index_of_species = data_->sorted_index[j];
//         break;
//       }
//     }
//   }

//   OpenSMOKE::OpenSMOKEVectorDouble x(data_->thermodynamicsMapXML->NumberOfSpecies());
//   OpenSMOKE::OpenSMOKEVectorDouble
//   omega(data_->thermodynamicsMapXML->NumberOfSpecies());
//   OpenSMOKE::OpenSMOKEVectorDouble c(data_->thermodynamicsMapXML->NumberOfSpecies());

//   std::vector<int> reaction_indices;
//   std::vector<double> reaction_coefficients;
//   // Local Analysis
//   if (ropaType_ == "local") {
//     unsigned int index = 0;
//     for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
//       if (data_->additional[data_->index_x_coord][j] >= local_x &&
//           data_->additional[data_->index_z_coord][j] >= local_z) {
//         index = j;
//         break;
//       }
//     }
//     // Recovers mass fractions
//     for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++)
//       omega[k + 1] = data_->omega[k][index];

//     // Calculates mole fractions
//     double MWmix;
//     data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MWmix,
//                                                                   omega.GetHandle());

//     // Calculates concentrations
//     const double P_Pa = data_->additional[data_->index_P][index];
//     const double T = data_->additional[data_->index_T][index];
//     const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
//     Product(cTot, x, &c);

//     // Calculates formations rates
//     data_->kineticsMapXML->SetTemperature(T);
//     data_->kineticsMapXML->SetPressure(P_Pa);
//     data_->thermodynamicsMapXML->SetTemperature(T);
//     data_->thermodynamicsMapXML->SetPressure(P_Pa);

//     data_->kineticsMapXML->KineticConstants();
//     data_->kineticsMapXML->ReactionRates(c.GetHandle());

//     // Ropa
//     OpenSMOKE::ROPA_Data ropa;
//     data_->kineticsMapXML->RateOfProductionAnalysis(ropa);

//     MergePositiveAndNegativeBars(ropa.production_reaction_indices[index_of_species],
//                                  ropa.destruction_reaction_indices[index_of_species],
//                                  ropa.production_coefficients[index_of_species],
//                                  ropa.destruction_coefficients[index_of_species],
//                                  reaction_indices, reaction_coefficients);
//   }  // Region
//   else if (ropaType_ == "region") {
//     std::vector<double> global_production_coefficients;
//     std::vector<double> global_destruction_coefficients;
//     std::vector<unsigned int> global_production_reaction_indices;
//     std::vector<unsigned int> global_destruction_reaction_indices;

//     unsigned int index_ll;
//     unsigned int index_lr;
//     unsigned int index_ul;
//     unsigned int index_ur;

//     for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
//       if (data_->additional[data_->index_x_coord][j] >= region_low_x &&
//           data_->additional[data_->index_z_coord][j] >= region_low_z) {
//         index_ll = j;
//         break;
//       }
//     }

//     for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
//       if (data_->additional[data_->index_x_coord][j] >= region_low_x &&
//           data_->additional[data_->index_z_coord][j] >= region_up_z) {
//         index_ul = j;
//         break;
//       }
//     }

//     for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
//       if (data_->additional[data_->index_x_coord][j] >= region_up_x &&
//           data_->additional[data_->index_z_coord][j] >= region_low_z) {
//         index_lr = j;
//         break;
//       }
//     }

//     for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
//       if (data_->additional[data_->index_x_coord][j] >= region_up_x &&
//           data_->additional[data_->index_z_coord][j] >= region_up_z) {
//         index_ur = j;
//         break;
//       }
//     }
//     // TODO ADD CHECK on idices
//     // if (index_min == index_max)
//     // {
//     //     if (index_max == data_->number_of_abscissas_ - 1)
//     //         index_min = index_max - 1;
//     //     else
//     //         index_max = index_min + 1;
//     // }
//     //

//     std::vector<double> volume = data_->additional[data_->index_volume];
//     double total_volume;
//     for (unsigned int i = 0; i < data_->number_of_abscissas_; i++) {
//       if (data_->additional[data_->index_x_coord][i] >=
//               data_->additional[data_->index_x_coord][index_ll] &&
//           data_->additional[data_->index_x_coord][i] <=
//               data_->additional[data_->index_x_coord][index_lr]) {
//         if (data_->additional[data_->index_z_coord][i] >=
//                 data_->additional[data_->index_z_coord][index_ll] &&
//             data_->additional[data_->index_z_coord][i] <=
//                 data_->additional[data_->index_z_coord][index_ul]) {
//           total_volume += volume[i];
//         }
//       }
//     }

//     int counter = 0;

//     for (unsigned int i = 0; i < data_->number_of_abscissas_; i++) {
//       if (data_->additional[data_->index_x_coord][i] >=
//               data_->additional[data_->index_x_coord][index_ll] &&
//           data_->additional[data_->index_x_coord][i] <=
//               data_->additional[data_->index_x_coord][index_lr]) {
//         if (data_->additional[data_->index_z_coord][i] >=
//                 data_->additional[data_->index_z_coord][index_ll] &&
//             data_->additional[data_->index_z_coord][i] <=
//                 data_->additional[data_->index_z_coord][index_ul]) {
//           // Recovers mass fractions
//           for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies();
//           k++)
//             omega[k + 1] = data_->omega[k][i];

//           // Calculates mole fractions
//           double MWmix;
//           data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(),
//           MWmix,
//                                                                         omega.GetHandle());

//           // Calculates concentrations
//           const double P_Pa = data_->additional[data_->index_P][i];
//           const double T = data_->additional[data_->index_T][i];
//           const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
//           Product(cTot, x, &c);

//           // Calculates formations rates
//           data_->kineticsMapXML->SetTemperature(T);
//           data_->kineticsMapXML->SetPressure(P_Pa);
//           data_->thermodynamicsMapXML->SetTemperature(T);
//           data_->thermodynamicsMapXML->SetPressure(P_Pa);

//           data_->kineticsMapXML->KineticConstants();
//           data_->kineticsMapXML->ReactionRates(c.GetHandle());

//           // Ropa
//           OpenSMOKE::ROPA_Data ropa;
//           data_->kineticsMapXML->RateOfProductionAnalysis(ropa);

//           if (ropa.production_coefficients[index_of_species].size() !=
//                   ropa.production_reaction_indices[index_of_species].size() ||
//               ropa.destruction_coefficients[index_of_species].size() !=
//                   ropa.destruction_reaction_indices[index_of_species].size()) {
//             throw std::invalid_argument("SSS");
//           }

//           if (counter == 0) {
//             global_production_coefficients.resize(
//                 ropa.production_coefficients[index_of_species].size());
//             global_destruction_coefficients.resize(
//                 ropa.destruction_coefficients[index_of_species].size());
//             global_production_reaction_indices =
//             ropa.production_reaction_indices[index_of_species];
//             global_destruction_reaction_indices =
//                 ropa.destruction_reaction_indices[index_of_species];
//           }

//           double integral_contribution = volume[i] / total_volume;

//           for (unsigned int k = 0; k <
//           ropa.production_coefficients[index_of_species].size(); k++)
//             global_production_coefficients[k] +=
//                 integral_contribution *
//                 ropa.production_coefficients[index_of_species][k];

//           for (unsigned int k = 0; k <
//           ropa.destruction_coefficients[index_of_species].size(); k++)
//             global_destruction_coefficients[k] +=
//                 integral_contribution *
//                 ropa.destruction_coefficients[index_of_species][k];

//           counter += 1;
//         }
//       }
//     }
//     MergePositiveAndNegativeBars(global_production_reaction_indices,
//                                  global_destruction_reaction_indices,
//                                  global_production_coefficients,
//                                  global_destruction_coefficients, reaction_indices,
//                                  reaction_coefficients);
//   } else  // Global
//   {
//     double tmp = 0;
//     unsigned int index_min = 0;
//     unsigned int index_max = data_->number_of_abscissas_;

//     std::vector<double> volume = data_->additional[data_->index_volume];
//     // double total_volume = std::accumulate(volume.begin(), volume.end(), 0.0d);
//     double total_volume = 0;
//     for (unsigned int i = 0; i < volume.size(); i++) {
//       total_volume += volume[i];
//     }

//     std::vector<double> global_production_coefficients;
//     std::vector<double> global_destruction_coefficients;
//     std::vector<unsigned int> global_production_reaction_indices;
//     std::vector<unsigned int> global_destruction_reaction_indices;

//     for (unsigned int j = index_min; j < index_max; j++) {
//       // Recovers mass fractions
//       for (unsigned int k = 0; k < data_->thermodynamicsMapXML->NumberOfSpecies(); k++)
//         omega[k + 1] = data_->omega[k][j];

//       // Calculates mole fractions
//       double MWmix;
//       data_->thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(),
//       MWmix,
//                                                                     omega.GetHandle());

//       // Calculates concentrations
//       const double P_Pa = data_->additional[data_->index_P][j];
//       const double T = data_->additional[data_->index_T][j];
//       const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
//       Product(cTot, x, &c);

//       // Calculates formations rates
//       data_->kineticsMapXML->SetTemperature(T);
//       data_->kineticsMapXML->SetPressure(P_Pa);
//       data_->thermodynamicsMapXML->SetTemperature(T);
//       data_->thermodynamicsMapXML->SetPressure(P_Pa);

//       data_->kineticsMapXML->KineticConstants();
//       data_->kineticsMapXML->ReactionRates(c.GetHandle());

//       // Ropa
//       OpenSMOKE::ROPA_Data ropa;
//       data_->kineticsMapXML->RateOfProductionAnalysis(ropa);

//       if (ropa.production_coefficients[index_of_species].size() !=
//               ropa.production_reaction_indices[index_of_species].size() ||
//           ropa.destruction_coefficients[index_of_species].size() !=
//               ropa.destruction_reaction_indices[index_of_species].size()) {
//         throw std::invalid_argument("SSS");
//       }

//       if (j == index_min) {
//         global_production_coefficients.resize(
//             ropa.production_coefficients[index_of_species].size());
//         global_destruction_coefficients.resize(
//             ropa.destruction_coefficients[index_of_species].size());
//         global_production_reaction_indices =
//         ropa.production_reaction_indices[index_of_species];
//         global_destruction_reaction_indices =
//         ropa.destruction_reaction_indices[index_of_species];
//       }

//       double integral_contribution = volume[j] / total_volume;

//       for (unsigned int k = 0; k <
//       ropa.production_coefficients[index_of_species].size(); k++)
//         global_production_coefficients[k] +=
//             integral_contribution * ropa.production_coefficients[index_of_species][k];

//       for (unsigned int k = 0; k <
//       ropa.destruction_coefficients[index_of_species].size(); k++)
//         global_destruction_coefficients[k] +=
//             integral_contribution * ropa.destruction_coefficients[index_of_species][k];
//     }

//     MergePositiveAndNegativeBars(global_production_reaction_indices,
//                                  global_destruction_reaction_indices,
//                                  global_production_coefficients,
//                                  global_destruction_coefficients, reaction_indices,
//                                  reaction_coefficients);
//   }

//   coefficients_.resize(std::min<int>(number_of_reactions,
//   reaction_coefficients.size())); reactions_.resize(std::min<int>(number_of_reactions,
//   reaction_coefficients.size()));

//   for (int i = 0; i < std::min<int>(number_of_reactions, reaction_coefficients.size());
//   i++) {
//     coefficients_[i] = reaction_coefficients[i];
//     reactions_[i] = reaction_indices[i];
//   }
// }
