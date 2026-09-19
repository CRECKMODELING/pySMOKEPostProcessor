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

PostProcessorCore::PostProcessorCore(void) {
  thermodynamicsMapXML = nullptr;
  kineticsMapXML = nullptr;
  thermodynamicsMapSurfaceXML = nullptr;
  kineticsMapSurfaceXML = nullptr;
  thermodynamicsMapLiquidXML = nullptr;
  kineticsMapLiquidXML = nullptr;

  iROPAEnabled_ = false;
  is_kinetics_available_ = false;

  iROPAHeterogeneousEnabled_ = false;
  is_heterogeneous_kinetics_available_ = false;

  phase2_kind_ = Phase2Kind::None;
}

PostProcessorCore::~PostProcessorCore(void) {}

bool PostProcessorCore::LoadKinetics(const std::string& folder_name) {
  path_folder_mechanism_ = folder_name;

  kinetics_gas_.reset(new KineticMapReader_Gas());
  kinetics_gas_->Load(folder_name);

  thermodynamicsMapXML = kinetics_gas_->thermodynamicsMap();
  kineticsMapXML = kinetics_gas_->kineticsMap();

  reaction_strings_ = kinetics_gas_->ReactionNames();

  is_kinetics_available_ = true;
  iROPAEnabled_ = true;

  LoadPhase2(folder_name);

  ValidateConsistency();

  return true;
}

void PostProcessorCore::LoadPhase2(const std::string& folder_name) {
  boost::filesystem::path folder(folder_name);

  if (boost::filesystem::exists(folder / "kinetics.surface.xml")) {
    phase2_kind_ = Phase2Kind::Surface;

    kinetics_surface_.reset(new KineticMapReader_Surface());
    kinetics_surface_->Load(folder_name);

    thermodynamicsMapSurfaceXML = kinetics_surface_->thermodynamicsMap();
    kineticsMapSurfaceXML = kinetics_surface_->kineticsMap();

    reaction_strings_heterogeneous_ = kinetics_surface_->ReactionNames();

    iROPAHeterogeneousEnabled_ = true;
    is_heterogeneous_kinetics_available_ = true;
  } else if (boost::filesystem::exists(folder / "kinetics.liquid.xml")) {
    phase2_kind_ = Phase2Kind::Liquid;

    kinetics_liquid_.reset(new KineticMapReader_Liquid());
    kinetics_liquid_->Load(folder_name);

    thermodynamicsMapLiquidXML = kinetics_liquid_->thermodynamicsMap();
    kineticsMapLiquidXML = kinetics_liquid_->kineticsMap();

    reaction_strings_heterogeneous_ = kinetics_liquid_->ReactionNames();

    iROPAHeterogeneousEnabled_ = true;
    is_heterogeneous_kinetics_available_ = true;
  } else if (boost::filesystem::exists(folder / "kinetics.solid.xml")) {
    phase2_kind_ = Phase2Kind::Solid;

    kinetics_solid_.reset(new KineticMapReader_Solid());
    kinetics_solid_->Load(folder_name);

    reaction_strings_heterogeneous_ = kinetics_solid_->ReactionNames();

    // No ROPA for solid - OpenSMOKEpp has no primitive to wrap - so the ROPA
    // enablement flags stay false even though kinetics did load successfully.
    is_heterogeneous_kinetics_available_ = true;
  } else {
    phase2_kind_ = Phase2Kind::None;
  }
}

KineticMapReaderBase* PostProcessorCore::phase2Kinetics() const {
  if (phase2_kind_ == Phase2Kind::Surface) return kinetics_surface_.get();
  if (phase2_kind_ == Phase2Kind::Liquid) return kinetics_liquid_.get();
  if (phase2_kind_ == Phase2Kind::Solid) return kinetics_solid_.get();
  return nullptr;
}

bool PostProcessorCore::ReadOutput(const std::string& folder_name, bool isHeterogeneous) {
  const bool ok = OutputReader::ReadOutput(folder_name, isHeterogeneous);
  ValidateConsistency();
  return ok;
}

bool PostProcessorCore::UpdateOutput(const std::string& folder_name) {
  const bool ok = OutputReader::UpdateOutput(folder_name);
  ValidateConsistency();
  return ok;
}

void PostProcessorCore::ValidateConsistency() {
  if (!HasKinetics() || !HasOutput()) return;

  if (kinetics_gas_->NumberOfSpecies() != omega.size()) {
    throw std::invalid_argument(
        "Output.xml file contains only a subset of the total species in the kinetic mechanism"
    );
  }
}

// 0-based. Forwards to kinetics_gas_/kinetics_surface_.
void PostProcessorCore::ReactionsAssociatedToSpecies(const unsigned int index,
                                                      std::vector<unsigned int>& indices) {
  kinetics_gas_->ReactionsAssociatedToSpecies(index, indices);
}

void PostProcessorCore::ReactionsAssociatedToSpecies_Surface(const unsigned int index,
                                                              std::vector<unsigned int>& indices) {
  kinetics_surface_->ReactionsAssociatedToSpecies(index, indices);
}

void PostProcessorCore::isReactantProduct(const unsigned int reaction_index,
                                          double& netStoichiometry) {
  kinetics_gas_->IsReactantProduct(reaction_index, netStoichiometry);
}

void PostProcessorCore::isReactantProduct_Surface(const unsigned int reaction_index,
                                                  double& netStoichiometry) {
  kinetics_surface_->IsReactantProduct(reaction_index, netStoichiometry);
}
