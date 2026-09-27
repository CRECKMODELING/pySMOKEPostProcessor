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

// Out-of-line ROPA::*_Liquid method definitions - RateOfProductionAnalysis_Liquid,
// GetReactionRates_Liquid, GetFormationRates_Liquid. Not a separate class: these are the
// liquid-kinetics bodies
// ROPA::RateOfProductionAnalysis/GetReactionRates/GetFormationRates dispatch to when
// heterogeneous_reactions=true and data_->phase2Kind() == Phase2Kind::Liquid (see ROPA.h
// and ROPA.hpp) - this used to be the entire ROPA_Liquid subclass, folded into ROPA the
// same way ROPA_Surface was (see ROPA_Surface.hpp).
//
// Unlike Surface (a genuinely different, composite gas+site+bulk state vector), a liquid
// mechanism's Output.xml is assumed to carry the same single-phase <mass-fractions>/
// <profiles> shape gas does - confirmed for kinetics.liquid.xml against a real mechanism,
// but no liquid Output.xml fixture exists anywhere to verify this against, so this file
// is structurally verified only (compiles, wired correctly), not exercised end-to-end.
// Because of that shape equivalence, the three methods below are ROPA's own single-phase
// logic with kineticsMapXML/thermodynamicsMapXML swapped for
// kineticsMapLiquidXML/thermodynamicsMapLiquidXML - not ROPA_Surface's composite-state
// logic, which doesn't apply here. FluxAnalysis stays inherited from ROPA unchanged
// (gas-only, same limitation Surface has - see ROPA_Surface.hpp).

#include <algorithm>

#include "math/OpenSMOKEUtilities.h"

void ROPA::RateOfProductionAnalysis_Liquid(const unsigned int number_of_reactions) {
  if (data_->phase2Kind() != Phase2Kind::Liquid) {
    throw std::invalid_argument(
        "ROPA::RateOfProductionAnalysis: heterogeneous_reactions=true requires a "
        "PostProcessorCore with liquid kinetics loaded");
  }

  // Select y variables among the species
  if (std::find(data_->string_list_massfractions_sorted.begin(),
                data_->string_list_massfractions_sorted.end(),
                species_) != data_->string_list_massfractions_sorted.end()) {
    speciesIsSelected = true;
  } else {
    throw std::invalid_argument("Please select one of the available species!");
  }

  unsigned int index_of_species;
  for (unsigned int j = 0; j < data_->thermodynamicsMapLiquidXML->NumberOfSpecies();
       j++) {
    if (speciesIsSelected == true) {
      if (species_ == data_->string_list_massfractions_sorted[j]) {
        index_of_species = data_->sorted_index[j];
        break;
      }
    }
  }

  OpenSMOKE::OpenSMOKEVectorDouble x(
      data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
  OpenSMOKE::OpenSMOKEVectorDouble omega(
      data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
  OpenSMOKE::OpenSMOKEVectorDouble c(
      data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
  // KineticsMap_Liquid_CHEMKIN::ReactionRates also takes the gas-phase
  // concentrations (for gas<->liquid interface reactions). No gas-phase
  // concentration profile is tracked in a liquid Output.xml (assumed
  // liquid-species-only, same as gas Output.xml is gas-species-only), so this
  // is zero throughout - no gas-phase species are treated as present.
  OpenSMOKE::OpenSMOKEVectorDouble c_gas(data_->thermodynamicsMapXML->NumberOfSpecies());
  c_gas = 0.;

  std::vector<int> reaction_indices;
  std::vector<double> reaction_coefficients;

  // Local Analysis
  if (ropaType_ == "local") {
    unsigned int index = 0;
    for (unsigned int j = 0; j < data_->number_of_abscissas_; j++) {
      if (data_->additional[0][j] >= localValue_) {
        index = j;
        break;
      }
    }
    // Recovers mass fractions
    for (unsigned int k = 0; k < data_->thermodynamicsMapLiquidXML->NumberOfSpecies();
         k++)
      omega[k + 1] = data_->omega[k][index];

    // Calculates mole fractions
    double MWmix;
    data_->thermodynamicsMapLiquidXML->MoleFractions_From_MassFractions(
        x.GetHandle(), MWmix, omega.GetHandle());

    // Calculates concentrations
    const double P_Pa = data_->additional[data_->index_P][index];
    const double T = data_->additional[data_->index_T][index];
    const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
    Product(cTot, x, &c);

    // Calculates formations rates
    data_->kineticsMapLiquidXML->SetTemperature(T);
    data_->kineticsMapLiquidXML->SetPressure(P_Pa);
    data_->thermodynamicsMapLiquidXML->SetTemperature(T);
    data_->thermodynamicsMapLiquidXML->SetPressure(P_Pa);

    data_->kineticsMapLiquidXML->KineticConstants();
    data_->kineticsMapLiquidXML->ReactionRates(c_gas.GetHandle(), c.GetHandle());

    // Ropa
    OpenSMOKE::ROPA_Data ropa;
    data_->kineticsMapLiquidXML->RateOfProductionAnalysis(ropa);

    MergePositiveAndNegativeBars(ropa.production_reaction_indices[index_of_species],
                                 ropa.destruction_reaction_indices[index_of_species],
                                 ropa.production_coefficients[index_of_species],
                                 ropa.destruction_coefficients[index_of_species],
                                 reaction_indices, reaction_coefficients);
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

    const double delta =
        data_->additional[0][index_max] - data_->additional[0][index_min];

    std::vector<double> global_production_coefficients;
    std::vector<double> global_destruction_coefficients;
    std::vector<unsigned int> global_production_reaction_indices;
    std::vector<unsigned int> global_destruction_reaction_indices;

    for (unsigned int j = index_min; j < index_max - 1; j++) {
      // Recovers mass fractions
      for (unsigned int k = 0; k < data_->thermodynamicsMapLiquidXML->NumberOfSpecies();
           k++)
        omega[k + 1] = data_->omega[k][j];

      // Calculates mole fractions
      double MWmix;
      data_->thermodynamicsMapLiquidXML->MoleFractions_From_MassFractions(
          x.GetHandle(), MWmix, omega.GetHandle());

      // Calculates concentrations
      const double P_Pa = data_->additional[data_->index_P][j];
      const double T = data_->additional[data_->index_T][j];
      const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
      Product(cTot, x, &c);

      // Calculates formations rates
      data_->kineticsMapLiquidXML->SetTemperature(T);
      data_->kineticsMapLiquidXML->SetPressure(P_Pa);
      data_->thermodynamicsMapLiquidXML->SetTemperature(T);
      data_->thermodynamicsMapLiquidXML->SetPressure(P_Pa);

      data_->kineticsMapLiquidXML->KineticConstants();
      data_->kineticsMapLiquidXML->ReactionRates(c_gas.GetHandle(), c.GetHandle());

      // Ropa
      OpenSMOKE::ROPA_Data ropa;
      data_->kineticsMapLiquidXML->RateOfProductionAnalysis(ropa);

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
        global_production_reaction_indices =
            ropa.production_reaction_indices[index_of_species];
        global_destruction_reaction_indices =
            ropa.destruction_reaction_indices[index_of_species];
      }

      const double dt = (data_->additional[0][j + 1] - data_->additional[0][j]) / delta;
      for (unsigned int k = 0; k < ropa.production_coefficients[index_of_species].size();
           k++)
        global_production_coefficients[k] +=
            dt * ropa.production_coefficients[index_of_species][k];

      for (unsigned int k = 0; k < ropa.destruction_coefficients[index_of_species].size();
           k++)
        global_destruction_coefficients[k] +=
            dt * ropa.destruction_coefficients[index_of_species][k];
    }

    MergePositiveAndNegativeBars(
        global_production_reaction_indices, global_destruction_reaction_indices,
        global_production_coefficients, global_destruction_coefficients, reaction_indices,
        reaction_coefficients);
  }

  coefficients_.resize(std::min<int>(number_of_reactions, reaction_coefficients.size()));
  reactions_.resize(std::min<int>(number_of_reactions, reaction_coefficients.size()));

  for (int i = 0; i < std::min<int>(number_of_reactions, reaction_coefficients.size());
       i++) {
    coefficients_[i] = reaction_coefficients[i];
    reactions_[i] = reaction_indices[i];
  }
}

void ROPA::GetReactionRates_Liquid(const std::vector<unsigned int>& reaction_indices,
                                   const bool sum_rates) {
  if (data_->phase2Kind() != Phase2Kind::Liquid) {
    throw std::invalid_argument(
        "ROPA::GetReactionRates: heterogeneous_reactions=true requires a "
        "PostProcessorCore with liquid kinetics loaded");
  }

  unsigned int numberOfReactions = reaction_indices.size();
  // Calculate the reaction rates
  {
    sumOfRates_.resize(data_->number_of_abscissas_);
    reactionRates_.resize(numberOfReactions,
                          std::vector<double>(data_->number_of_abscissas_, 1));

    OpenSMOKE::OpenSMOKEVectorDouble x(
        data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble omega(
        data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble c(
        data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble r(data_->kineticsMapLiquidXML->NumberOfReactions());
    // See RateOfProductionAnalysis_Liquid() above for why this is zero throughout.
    OpenSMOKE::OpenSMOKEVectorDouble c_gas(
        data_->thermodynamicsMapXML->NumberOfSpecies());
    c_gas = 0.;

    for (unsigned int i = 0; i < data_->number_of_abscissas_; i++) {
      // Recovers mass fractions
      for (unsigned int k = 0; k < data_->thermodynamicsMapLiquidXML->NumberOfSpecies();
           k++)
        omega[k + 1] = data_->omega[k][i];

      // Calculate mole fractions
      double MWmix;
      data_->thermodynamicsMapLiquidXML->MoleFractions_From_MassFractions(
          x.GetHandle(), MWmix, omega.GetHandle());

      // Calculate concentrations
      const double P_Pa = data_->additional[data_->index_P][i];
      const double T = data_->additional[data_->index_T][i];
      const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
      Product(cTot, x, &c);

      // Calculate reaction rates
      data_->kineticsMapLiquidXML->SetTemperature(T);
      data_->kineticsMapLiquidXML->SetPressure(P_Pa);
      data_->thermodynamicsMapLiquidXML->SetTemperature(T);
      data_->thermodynamicsMapLiquidXML->SetPressure(P_Pa);

      data_->kineticsMapLiquidXML->ReactionRates(c_gas.GetHandle(), c.GetHandle());
      data_->kineticsMapLiquidXML->GiveMeReactionRates(r.GetHandle());

      if (sum_rates) {
        double sum_rate = 0.;
        for (unsigned int k = 0; k < numberOfReactions; k++) {
          const unsigned int j = reaction_indices[k] + 1;
          sum_rate += r[j];
        }
        sumOfRates_[i] = sum_rate;
      } else {
        for (unsigned int k = 0; k < numberOfReactions; k++) {
          const unsigned int j = reaction_indices[k] + 1;
          reactionRates_[k][i] = r[j];
        }
      }
    }
  }
}

void ROPA::GetFormationRates_Liquid(const std::string& specie, const std::string& units,
                                    const std::string& type) {
  if (data_->phase2Kind() != Phase2Kind::Liquid) {
    throw std::invalid_argument(
        "ROPA::GetFormationRates: heterogeneous_reactions=true requires a "
        "PostProcessorCore with liquid kinetics loaded");
  }
  if (units != "mass" && units != "mole")
    throw std::invalid_argument("Available Formation Rates units are: mole | mass");

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

    OpenSMOKE::OpenSMOKEVectorDouble P(
        data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble D(
        data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble x(
        data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble omega(
        data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
    OpenSMOKE::OpenSMOKEVectorDouble c(
        data_->thermodynamicsMapLiquidXML->NumberOfSpecies());
    // See RateOfProductionAnalysis_Liquid() above for why this is zero throughout.
    OpenSMOKE::OpenSMOKEVectorDouble c_gas(
        data_->thermodynamicsMapXML->NumberOfSpecies());
    c_gas = 0.;

    for (unsigned int i = 0; i < data_->number_of_abscissas_; i++) {
      // Recovers mass fractions
      for (unsigned int k = 0; k < data_->thermodynamicsMapLiquidXML->NumberOfSpecies();
           k++)
        omega[k + 1] = data_->omega[k][i];

      // Calculates mole fractions
      double MWmix;
      data_->thermodynamicsMapLiquidXML->MoleFractions_From_MassFractions(
          x.GetHandle(), MWmix, omega.GetHandle());

      // Calculates concentrations
      const double P_Pa = data_->additional[data_->index_P][i];
      const double T = data_->additional[data_->index_T][i];
      const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
      Product(cTot, x, &c);

      // Calculates formations rates
      data_->kineticsMapLiquidXML->SetTemperature(T);
      data_->kineticsMapLiquidXML->SetPressure(P_Pa);
      data_->thermodynamicsMapLiquidXML->SetTemperature(T);
      data_->thermodynamicsMapLiquidXML->SetPressure(P_Pa);

      data_->kineticsMapLiquidXML->KineticConstants();
      data_->kineticsMapLiquidXML->ReactionRates(c_gas.GetHandle(), c.GetHandle());
      data_->kineticsMapLiquidXML->ProductionAndDestructionRates(
          P.GetHandle(),
          D.GetHandle());  // kmol/m3/s

      if (type == "characteristic-time") {
        const unsigned k = data_->sorted_index[formation_rates_to_plot[1]] + 1;
        formationRates_[i] = c[k] / (D[k] + 1.e-32);
      } else {
        if (units == "mass") {
          OpenSMOKE::ElementByElementProduct(
              P.Size(), P.GetHandle(), data_->thermodynamicsMapLiquidXML->MWs().data(),
              P.GetHandle());
          OpenSMOKE::ElementByElementProduct(
              D.Size(), D.GetHandle(), data_->thermodynamicsMapLiquidXML->MWs().data(),
              D.GetHandle());
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
