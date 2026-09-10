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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <queue>
#include <sstream>
#include <iostream>

SpeciesClass::SpeciesClass() {
  data_ = nullptr;
  species_classes_available_ = false;

  species_ = "";
  class_name_ = "";
  element_ = "C";
  flux_type_ = "destruction";
  thickness_ = "relative";
  label_type_ = "relative";
  thicknesslogscale_ = true;
  carbon_weighted_ = true;
  flux_per_class_ = true;
  auto_prune_diagonal_ = true;
  width_ = 3;
  depth_ = 2;
  threshold_ = 0.;
  localValue_ = 0.;
}

void SpeciesClass::SetDatabase(ProfilesDatabase* data) {
  data_ = data;
  ReadSpeciesClasses();
}

void SpeciesClass::SetFluxPerClass(const bool flux_per_class) { flux_per_class_ = flux_per_class; }

void SpeciesClass::SetClassName(const std::string class_name) { class_name_ = class_name; }

void SpeciesClass::SetSpecies(const std::string species) { species_ = species; }

void SpeciesClass::SetElement(const std::string element) { element_ = element; }

void SpeciesClass::SetFluxAnalysisType(const std::string type) {
  if (type != "production" && type != "destruction")
    throw std::invalid_argument("Available flux analysis types are: production | destruction");
  flux_type_ = type;
}

void SpeciesClass::SetLocalValue(const double localValue) { localValue_ = localValue; }

void SpeciesClass::SetThickness(const std::string thickness) {
  if (thickness != "absolute" && thickness != "relative")
    throw std::invalid_argument("Available thickness types are: absolute | relative(%)");
  thickness_ = thickness;
}

void SpeciesClass::SetLabelType(const std::string type) {
  if (type != "absolute" && type != "relative")
    throw std::invalid_argument("Available label types are: absolute | relative(%)");
  label_type_ = type;
}

void SpeciesClass::SetThicknessLogScale(const bool thicknesslogscale) {
  thicknesslogscale_ = thicknesslogscale;
}

void SpeciesClass::SetWidth(const int width) { width_ = width; }

void SpeciesClass::SetDepth(const int depth) { depth_ = depth; }

void SpeciesClass::SetThreshold(const double threshold) { threshold_ = threshold; }

void SpeciesClass::SetCarbonWeighted(const bool carbon_weighted) { carbon_weighted_ = carbon_weighted; }

void SpeciesClass::SetAutoPruneDiagonal(const bool auto_prune_diagonal) {
  auto_prune_diagonal_ = auto_prune_diagonal;
}

// The <SpeciesClasses> block is parsed in ProfilesDatabase::ReadKineticMechanism
void SpeciesClass::ReadSpeciesClasses() {
  species_classes_available_ = false;
  class_names_.clear();
  class_species_.clear();
  species_to_class_.clear();

  if (data_ == nullptr) return;

  class_names_ = data_->species_class_names_;
  class_species_ = data_->species_class_members_;
  species_to_class_ = data_->species_to_class_;
  species_classes_available_ = data_->has_species_classes_;
}

unsigned int SpeciesClass::ElementIndex(const std::string& element) const {
  const std::vector<std::string> elements_names = data_->thermodynamicsMapXML->elements();
  for (unsigned int k = 0; k < elements_names.size(); k++)
    if (element == elements_names[k]) return k;
  throw std::invalid_argument("The element '" + element + "' is not present in the mechanism");
}

void SpeciesClass::ElementalDistribution(const std::string element, const bool normalize) {
  if (!species_classes_available_)
    throw std::invalid_argument("The kinetic mechanism does not contain a <SpeciesClasses> block");

  const unsigned int index_element = ElementIndex(element);
  const unsigned int npts = data_->number_of_abscissas_;
  const std::vector<double>& mws = data_->thermodynamicsMapXML->MWs();
  const Eigen::MatrixXd& elemental_composition = data_->thermodynamicsMapXML->atomic_composition();

  abscissa_ = data_->additional[0];

  const unsigned int nclasses = class_names_.size();
  elemental_fractions_.assign(nclasses, std::vector<double>(npts, 0.));

  for (unsigned int k = 0; k < nclasses; k++) {
    for (unsigned int s_i = 0; s_i < class_species_[k].size(); s_i++) {
      const unsigned int s = class_species_[k][s_i];
      const double n_atoms = elemental_composition(s, index_element);
      if (n_atoms == 0.) continue;
      const double moles_of_element_per_kg = n_atoms / mws[s];
      for (unsigned int i = 0; i < npts; i++)
        elemental_fractions_[k][i] += data_->omega[s][i] * moles_of_element_per_kg;
    }
  }

  if (normalize) {
    for (unsigned int i = 0; i < npts; i++) {
      double total = 0.;
      for (unsigned int k = 0; k < nclasses; k++) total += elemental_fractions_[k][i];
      if (total > 0.)
        for (unsigned int k = 0; k < nclasses; k++) elemental_fractions_[k][i] /= total;
    }
  }
}

// Class -> class element-flux matrix, computed straight from every per-reaction flux. 
// Two weightings are available via carbon_weighted_: 
// the per-reaction "carbon share" routing, and the OpenSMOKE carbon-atom throughput of
// maps/FluxAnalysisMap.hpp::AnalyzeNetFluxes. 
// The matrix is always available via fluxMatrix(); 
// width_/depth_/threshold_ only are used to simplify the class graph, 
// coming straight from the seed class.
void SpeciesClass::FluxByClass() {
  if (!species_classes_available_)
    throw std::invalid_argument("The kinetic mechanism does not contain a <SpeciesClasses> block");

  OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermo = *data_->thermodynamicsMapXML;
  OpenSMOKE::KineticsMap_CHEMKIN& kinetics = *data_->kineticsMapXML;

  const unsigned int ns = thermo.NumberOfSpecies();
  const unsigned int nr = kinetics.NumberOfReactions();
  const unsigned int nclasses = class_names_.size();
  const unsigned int index_element = ElementIndex(element_);
  const Eigen::MatrixXd& elemental_composition = thermo.atomic_composition();

  std::vector<double> n_atom(ns);
  for (unsigned int s = 0; s < ns; s++) n_atom[s] = elemental_composition(s, index_element);

  // --- thermochemical state at the requested point of the independent variable ---
  unsigned int index = 0;
  for (unsigned int j = 0; j < static_cast<unsigned int>(data_->number_of_abscissas_); j++)
    if (data_->additional[0][j] >= localValue_) {
      index = j;
      break;
    }

  OpenSMOKE::OpenSMOKEVectorDouble x(ns);
  OpenSMOKE::OpenSMOKEVectorDouble omega(ns);
  OpenSMOKE::OpenSMOKEVectorDouble c(ns);
  OpenSMOKE::OpenSMOKEVectorDouble r(nr);

  for (unsigned int k = 0; k < ns; k++) omega[k + 1] = data_->omega[k][index];

  double MWmix;
  thermo.MoleFractions_From_MassFractions(x.GetHandle(), MWmix, omega.GetHandle());

  const double P_Pa = data_->additional[data_->index_P][index];
  const double T = data_->additional[data_->index_T][index];
  const double cTot = P_Pa / PhysicalConstants::R_J_kmol / T;
  Product(cTot, x, &c);

  kinetics.SetTemperature(T);
  kinetics.SetPressure(P_Pa);
  thermo.SetTemperature(T);
  thermo.SetPressure(P_Pa);
  kinetics.KineticConstants();
  kinetics.ReactionRates(c.GetHandle());
  kinetics.GiveMeReactionRates(r.GetHandle());  // signed net rates, 1-based

  // --- reactant / product species carrying `element`, per reaction ---
  kinetics.stoichiometry().BuildStoichiometricMatrix();
  const Eigen::SparseMatrix<double>& SR = kinetics.stoichiometry().stoichiometric_matrix_reactants();
  const Eigen::SparseMatrix<double>& SP = kinetics.stoichiometry().stoichiometric_matrix_products();

  std::vector<std::vector<unsigned int>> reactants_of(nr), products_of(nr);
  std::vector<double> element_atoms_reactants(nr, 0.);  // OpenSMOKE convention: no stoich multiplicity
  std::vector<double> element_atoms_products(nr, 0.);
  for (int k = 0; k < SR.outerSize(); ++k)              // k = species (column)
    for (Eigen::SparseMatrix<double>::InnerIterator it(SR, k); it; ++it) {
      const unsigned int rr = it.row();                 // reaction
      element_atoms_reactants[rr] += n_atom[k];
      if (n_atom[k] != 0.) reactants_of[rr].push_back(static_cast<unsigned int>(k));
    }
  for (int k = 0; k < SP.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(SP, k); it; ++it) {
      const unsigned int rr = it.row();
      element_atoms_products[rr] += n_atom[k];
      if (n_atom[k] != 0.) products_of[rr].push_back(static_cast<unsigned int>(k));
    }

  // class -> class element-flux matrix
  //   phi_ij, the contribution of reaction r to the flux from the class of
  //   reactant i to the class of product j, splits reaction r among its
  //   (C-reactant, C-product) pairs by carbon share w_i = n_i / N_reactants,
  //   w_j = n_j / N_products:
  //     carbon_weighted_ == true (default, OpenSMOKE): phi = |R_r| * n_i * n_j / N_reactants
  //         -> sum over ij ~= |R_r| * N_C,r ; the carbon-atom throughput. Faithful
  //            to maps/FluxAnalysisMap.hpp::AnalyzeNetFluxes. Note (bin C count)^2
  //            makes soot-soot terms dominate a mechanism with soot.
  //     carbon_weighted_ == false: phi = |R_r| * w_i * w_j
  //         -> sum over ij == |R_r| ; a per-reaction molar rate routed by carbon
  //            share. Bounded, so BINs stay comparable to small species.
  flux_matrix_.assign(nclasses, std::vector<double>(nclasses, 0.));
  for (unsigned int rr = 0; rr < nr; ++rr) {
    const double R = r[rr + 1];
    const double Nr = element_atoms_reactants[rr];
    const double Np = element_atoms_products[rr];
    if (R == 0. || Nr == 0. || Np == 0.) continue;
    double norm;
    if (carbon_weighted_) {
      norm = std::fabs(R) / Nr;
    } else {
      norm = std::fabs(R) / (Nr * Np); // stupid soot.
    }
    for (unsigned int a = 0; a < reactants_of[rr].size(); ++a) {
      const unsigned int i = reactants_of[rr][a];
      const int ci = species_to_class_[i];
      if (ci < 0) continue;
      for (unsigned int b = 0; b < products_of[rr].size(); ++b) {
        const unsigned int j = products_of[rr][b];
        if (j == i) continue;  // species unchanged by the reaction: no flux (as in OpenSMOKE)
        const int cj = species_to_class_[j];
        if (cj < 0) continue;
        const double phi = n_atom[i] * n_atom[j] * norm;
        if (R > 0.)
          flux_matrix_[ci][cj] += phi;  // Forward reaction
        else
          flux_matrix_[cj][ci] += phi;  // Reversible reaction going backwards
      }
    }
  }

  // Intra-class flux (i == j) is typically very large (e.g. HABS reactions),
  // it is therefore set to zero for visualization purposes with heatmaps
  if (auto_prune_diagonal_) {
    for (unsigned int k = 0; k < nclasses; k++) flux_matrix_[k][k] = 0.;
  }

  index_first_class_.clear();
  index_second_class_.clear();
  computed_label_.clear();
  computed_thickness_.clear();

  int seed_class = -1;
  if (flux_per_class_) {
    for (unsigned int k = 0; k < class_names_.size(); k++)
      if (class_names_[k] == class_name_) {
        seed_class = static_cast<int>(k);
        break;
      }
    if (seed_class < 0)
      throw std::invalid_argument("The species class '" + class_name_ +
                                  "' is not defined in this mechanism");
  } else if (!species_.empty()) {
    const std::vector<std::string>& names = thermo.NamesOfSpecies();
    for (unsigned int s = 0; s < names.size(); s++)
      if (names[s] == species_) {
        seed_class = species_to_class_[s];
        break;
      }
    if (seed_class < 0)
      throw std::invalid_argument("The species '" + species_ +
                                  "' is not in the mechanism or has no species class");
  }
  if (seed_class < 0) return;  // no seed: the full matrix is still available via fluxMatrix()

  const bool follow_incoming = (flux_type_ == "production");  // where the element comes from

  // Used for different production/destruction thickness evaluation
  std::vector<int> edge_source, edge_target;
  std::vector<double> edge_flux, edge_node_total;

  std::vector<char> visited(nclasses, 0);
  std::queue<std::pair<int, int>> frontier;  // (class, depth)
  frontier.push(std::make_pair(seed_class, 0));
  visited[seed_class] = 1;

  while (!frontier.empty()) {
    const int cur = frontier.front().first;
    const int d = frontier.front().second;
    frontier.pop();
    if (d >= depth_) continue;

    std::vector<std::pair<double, int>> cand;  // (flux, neighbour class)
    double node_total = 0.;
    for (unsigned int k = 0; k < nclasses; k++) {
      if (static_cast<int>(k) == cur) continue;
      double f;
      if (follow_incoming) {
        f = flux_matrix_[k][cur];  // carbon flowing INTO cur from k
      } else {
        f = flux_matrix_[cur][k];  // carbon flowing OUT OF cur into k
      }
      if (f <= 0.) continue;
      cand.push_back(std::make_pair(f, static_cast<int>(k)));
      node_total += f;
    }
    if (node_total == 0.) continue;
    std::sort(cand.begin(), cand.end(),
              [](const std::pair<double, int>& p, const std::pair<double, int>& q) {
                return p.first > q.first;
              });

    int kept = 0;
    for (unsigned int m = 0; m < cand.size() && kept < width_; ++m) {
      const double f = cand[m].first;
      if (f / node_total * 100. < threshold_) break;  // candidates are sorted descending
      const int nb = cand[m].second;
      ++kept;

      edge_source.push_back(cur);
      edge_target.push_back(nb);
      edge_flux.push_back(f);
      edge_node_total.push_back(node_total);

      if (!visited[nb]) {
        visited[nb] = 1;
        frontier.push(std::make_pair(nb, d + 1));
      }
    }
  }

  //   label_type_ == "relative": % of the SOURCE class's own production/destruction
  //                  "absolute": the raw flux value.
  //   thickness_  == "relative": the arrow width follows that very same per-node percentage, 
  //                              so a 96% edge is drawn (almost) as thick as physically possible
  //                  "absolute": the arrow width follows the edge's flux relative to the single 
  //                              largest edge actually present in this graph (not the
  //                              whole class x class matrix, most of which isn't even drawn).
  double local_max_flux = 0.;
  if (thickness_ == "absolute") {
    for (unsigned int e = 0; e < edge_flux.size(); e++)
      local_max_flux = std::max(local_max_flux, edge_flux[e]);
    if (local_max_flux == 0.) local_max_flux = 1.;
  }

  for (unsigned int e = 0; e < edge_flux.size(); e++) {
    index_first_class_.push_back(edge_source[e]);
    index_second_class_.push_back(edge_target[e]);

    double label_value;
    if (label_type_ == "relative") {
      label_value = edge_flux[e] / edge_node_total[e] * 100.;
    } else {
      label_value = edge_flux[e];
    }
    computed_label_.push_back(label_value);

    double t;
    if (thickness_ == "relative") {
      t = edge_flux[e] / edge_node_total[e];
    } else {
      t = edge_flux[e] / local_max_flux;
    }

    double thickness_value;
    if (thicknesslogscale_) {
      thickness_value = 0.5 + 5. * std::max(0., 3. + std::log10(t + 1.e-30));
    } else {
      thickness_value = 0.5 + 15. * t;
    }
    computed_thickness_.push_back(thickness_value);
  }
}
