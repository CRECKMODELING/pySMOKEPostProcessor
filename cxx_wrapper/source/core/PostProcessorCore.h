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

#ifndef POSTPROCESSORCORE_H
#define POSTPROCESSORCORE_H

#include <Eigen/Sparse>
#include <memory>

#include "KineticMapReader_Gas.h"
#include "KineticMapReader_Liquid.h"
#include "KineticMapReader_Solid.h"
#include "KineticMapReader_Surface.h"
#include "OutputReader.h"

// A mechanism has gas + at most one heterogeneous phase kinetics,
// OpenSMOKE (v0.24) gives an error at mechanism compilation time
// if two heterogeneous phases are selected.
// LoadKinetics auto-detects this by which kinetics file is present 
// (kinetics.surface.xml / kinetics.liquid.xml / kinetics.solid.xml)
// and records which one here.
enum class Phase2Kind { None, Surface, Liquid, Solid };

// This is the main working class in the C++ part of the PostProcessor.
// All fields are publicly accessible
// (data_->kineticsMapXML, data_->omega, data_->gasKinetics()->soot_bin_*, ...) -
// inherited from OutputReader for the output half; most kinetics-half fields
// are copied out of the owned KineticMapReader_* objects, except phase-specific
// data (e.g. soot, gas-only) that are read straight off gasKinetics()/etc.
// instead of duplicating it here.
//
// Construction is willingly partial: LoadKinetics()
// and ReadOutput()/UpdateOutput() (inherited from OutputReader) are
// independent calls, and a caller may use only one of them. 
// This is intended for possible uses where the user wants to generate
// the PostProcessor once using a fixed kinetic mechanism and then update
// the Output multiple times, e.g. for parsing results from a ParametricAnalysis.
// HasKinetics()/HasOutput()/HasPhase2Kinetics() report what is actually loaded.

// Why this exists: codes even in the examples were rebuilding the pp object from scratch. 
// This is not efficient, and this modification lead to a rough /4 reduction in time 
// for the file ReactionClasses_FDI (ca 8-10 min to 2)
// Another case: the DepositionPlot requires a long series of calls to the pp at different
// locations, and the required memory was building up over time, reaching above 8 GB of RAM.
class PostProcessorCore : public OutputReader {
 public:
  PostProcessorCore(void);
  ~PostProcessorCore(void);

  bool LoadKinetics(const std::string& folder_name);

  // Hide the base versions so calls through a PostProcessorCore also run
  // ValidateConsistency() once both sides are loaded. 
  // Every call passes through a PostProcessorCore*, never an OutputReader*
  bool ReadOutput(const std::string& folder_name, bool isHeterogeneous);
  bool UpdateOutput(const std::string& folder_name);

  bool HasKinetics() const { return is_kinetics_available_; }
  bool HasPhase2Kinetics() const { return phase2_kind_ != Phase2Kind::None; }
  Phase2Kind phase2Kind() const { return phase2_kind_; }

  KineticMapReader_Gas* gasKinetics() const { return kinetics_gas_.get(); }

  // Parses the heterogeneous kinetics, whichever phase is available
  // The check on which type of phase that is is performed automatically.
  KineticMapReaderBase* phase2Kinetics() const;

  void ValidateConsistency();

  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML;
  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML;
  OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN* thermodynamicsMapSurfaceXML;
  OpenSMOKE::KineticsMap_Surface_CHEMKIN* kineticsMapSurfaceXML;
  OpenSMOKE::ThermodynamicsMap_Liquid_CHEMKIN* thermodynamicsMapLiquidXML;
  OpenSMOKE::KineticsMap_Liquid_CHEMKIN* kineticsMapLiquidXML;

  bool iROPAEnabled_;
  bool is_kinetics_available_;

  bool iROPAHeterogeneousEnabled_;
  bool is_heterogeneous_kinetics_available_;

  boost::filesystem::path path_folder_mechanism_;

  void ReactionsAssociatedToSpecies(const unsigned int index, std::vector<unsigned int>& indices);
  void isReactantProduct(const unsigned int reaction_index, double& netStoichiometry);

  void ReactionsAssociatedToSpecies_Surface(const unsigned int index, std::vector<unsigned int>& indices);
  void isReactantProduct_Surface(const unsigned int reaction_index, double& netStoichiometry);

  std::string name_reactions_;
  std::vector<std::string> reaction_strings_;

  std::string name_reactions_heterogeneous_;
  std::vector<std::string> reaction_strings_heterogeneous_;

  // ---- Optional mechanism blocks
  // Species-class (<SpeciesClasses>) and reaction-class (<ReactionClasses>)
  // metadata are NOT mirrored here:
  // HasSpeciesClasses()/SpeciesClassNames()/SpeciesClassMembers()/SpeciesToClass()
  // and HasReactionClasses()/ReactionMainClass()/ReactionSubClass() are already
  // public reference-returning accessors on KineticMapReaderBase (the common
  // base of every phase reader), so SpeciesClass/ReactionClass read them
  // straight off gasKinetics() (gas side) / phase2Kinetics() (heterogeneous
  // side) instead of PostProcessorCore holding its own copy of each.

  // Soot BIN properties (<SootProperties>, generated at mechanism compile-time
  // with the keyword @SootProperties) are NOT mirrored here: soot is always
  // gas-phase (KineticMapReader_Gas.h), and its fields are already public
  // there, so Soot.hpp reads them straight off gasKinetics()->soot_bin_*
  // instead of paying for a second copy.

  // Owns the actual kinetics.xml/kinetics.{surface,liquid,solid}.xml parse.
  // The raw pointers/vectors above are copies of what these already hold -
  // kept so widgets reading data_->kineticsMapXML etc. directly don't need
  // to change. At most one of kinetics_surface_/kinetics_liquid_/kinetics_solid_
  // is ever non-null, matching phase2_kind_.
  std::unique_ptr<KineticMapReader_Gas> kinetics_gas_;
  std::unique_ptr<KineticMapReader_Surface> kinetics_surface_;
  std::unique_ptr<KineticMapReader_Liquid> kinetics_liquid_;
  std::unique_ptr<KineticMapReader_Solid> kinetics_solid_;

 private:
  Phase2Kind phase2_kind_;

  void LoadPhase2(const std::string& folder_name);
};

#include "PostProcessorCore.hpp"
#endif  // POSTPROCESSORCORE_H
