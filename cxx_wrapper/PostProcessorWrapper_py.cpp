#include "PostProcessorWrapper_py.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "source/classes/ReactionClasses.h"
#include "source/classes/Soot.h"
#include "source/classes/SpeciesClasses.h"
#include "source/core/KineticMapReaderBase.h"
#include "source/core/KineticMapReader_Gas.h"
#include "source/core/KineticMapReader_Liquid.h"
#include "source/core/KineticMapReader_Solid.h"
#include "source/core/KineticMapReader_Surface.h"
#include "source/core/PostProcessorCore.h"
#include "source/ropa/ROPA.h"
#include "source/sensitivity/SensitivityCalculator.h"

namespace py = pybind11;
constexpr auto byref = py::return_value_policy::reference_internal;

PYBIND11_MODULE(pySMOKEPostProcessor, m) {
  m.doc() = "Python Interface to the OpenSMOKEpp Graphical Post Processor";

  py::enum_<Phase2Kind>(m, "Phase2Kind")
      .value("None_", Phase2Kind::None)
      .value("Surface", Phase2Kind::Surface)
      .value("Liquid", Phase2Kind::Liquid)
      .value("Solid", Phase2Kind::Solid);

  // Common read-only surface shared by every phase's kinetics reader (name<->index
  // lookups, real per-species MW, the optional SpeciesClasses/ReactionClasses
  // blocks). No py::init<>() - Python never constructs these directly, only
  // reaches them via PostProcessorCore.gasKinetics()/phase2Kinetics(). Registering
  // the base (it has a virtual destructor) is what lets pybind11 automatically
  // downcast phase2Kinetics()'s KineticMapReaderBase* to whichever concrete
  // Surface/Liquid/Solid type is actually loaded.
  py::class_<KineticMapReaderBase>(m, "KineticMapReaderBase")
      .def("numberOfReactions", &KineticMapReaderBase::NumberOfReactions,
           py::call_guard<py::gil_scoped_release>())
      .def("numberOfSpecies", &KineticMapReaderBase::NumberOfSpecies,
           py::call_guard<py::gil_scoped_release>())
      .def("reactionNameFromIndex", &KineticMapReaderBase::ReactionNameFromIndex,
           py::arg("index"), py::call_guard<py::gil_scoped_release>())
      .def("formattedReactionNameFromIndex",
           &KineticMapReaderBase::FormattedReactionNameFromIndex, py::arg("index"),
           py::call_guard<py::gil_scoped_release>())
      .def("speciesNameFromIndex", &KineticMapReaderBase::SpeciesNameFromIndex,
           py::arg("index"), py::call_guard<py::gil_scoped_release>())
      .def("reactionIndexFromName", &KineticMapReaderBase::ReactionIndexFromName,
           py::arg("name"), py::call_guard<py::gil_scoped_release>())
      .def("speciesIndexFromName", &KineticMapReaderBase::SpeciesIndexFromName,
           py::arg("name"), py::call_guard<py::gil_scoped_release>())
      .def("speciesNames", &KineticMapReaderBase::SpeciesNames,
           py::call_guard<py::gil_scoped_release>())
      .def("hasSpeciesClasses", &KineticMapReaderBase::HasSpeciesClasses,
           py::call_guard<py::gil_scoped_release>())
      .def("speciesClassNames", &KineticMapReaderBase::SpeciesClassNames,
           py::call_guard<py::gil_scoped_release>())
      .def("speciesClassMembers", &KineticMapReaderBase::SpeciesClassMembers,
           py::call_guard<py::gil_scoped_release>())
      .def("mw", &KineticMapReaderBase::MW, py::arg("index"),
           py::call_guard<py::gil_scoped_release>());

  // stoichiometricMatrixReactants/Products return Eigen::SparseMatrix<double>,
  // converted to scipy.sparse automatically by pybind11/eigen.h - replaces the
  // from-scratch pySMOKEPostProcessor/maps/StoichiometricMap.py. Unused, but exists as
  // stoichiometric matrix interface.
  py::class_<KineticMapReader_Gas, KineticMapReaderBase>(m, "KineticMapReader_Gas")
      .def("stoichiometricMatrixReactants",
           &KineticMapReader_Gas::StoichiometricMatrixReactants,
           py::call_guard<py::gil_scoped_release>())
      .def("stoichiometricMatrixProducts",
           &KineticMapReader_Gas::StoichiometricMatrixProducts,
           py::call_guard<py::gil_scoped_release>());

  py::class_<KineticMapReader_Surface, KineticMapReaderBase>(m,
                                                             "KineticMapReader_Surface")
      .def("stoichiometricMatrixReactants",
           &KineticMapReader_Surface::StoichiometricMatrixReactants,
           py::call_guard<py::gil_scoped_release>())
      .def("stoichiometricMatrixProducts",
           &KineticMapReader_Surface::StoichiometricMatrixProducts,
           py::call_guard<py::gil_scoped_release>());

  // Liquid/Solid: structurally verified only (no example data)
  py::class_<KineticMapReader_Liquid, KineticMapReaderBase>(m, "KineticMapReader_Liquid");
  py::class_<KineticMapReader_Solid, KineticMapReaderBase>(m, "KineticMapReader_Solid");

  py::class_<PostProcessorCore>(m, "PostProcessorCore")
      .def(py::init<>())
      .def("loadKinetics", &PostProcessorCore::LoadKinetics,
           py::arg("folder_name") = "kinetics", py::call_guard<py::gil_scoped_release>())
      .def("readFileResults", &PostProcessorCore::ReadOutput,
           py::arg("folder_name") = "Output", py::arg("isHeterogeneous") = false,
           py::call_guard<py::gil_scoped_release>())
      .def("updateOutput", &PostProcessorCore::UpdateOutput, py::arg("folder_name"),
           py::call_guard<py::gil_scoped_release>())
      .def("hasKinetics", &PostProcessorCore::HasKinetics,
           py::call_guard<py::gil_scoped_release>())
      .def("hasOutput", &PostProcessorCore::HasOutput,
           py::call_guard<py::gil_scoped_release>())
      .def("hasPhase2Kinetics", &PostProcessorCore::HasPhase2Kinetics,
           py::call_guard<py::gil_scoped_release>())
      .def("phase2Kind", &PostProcessorCore::phase2Kind,
           py::call_guard<py::gil_scoped_release>())
      .def("gasKinetics", &PostProcessorCore::gasKinetics, byref,
           py::call_guard<py::gil_scoped_release>())
      .def("phase2Kinetics", &PostProcessorCore::phase2Kinetics, byref,
           py::call_guard<py::gil_scoped_release>())
      // --- OutputReader's read surface (kinetics-independent, see OutputReader.h) ---
      .def("getSpeciesProfile", &PostProcessorCore::GetSpeciesProfile, py::arg("name"),
           py::arg("basis") = "mole", py::call_guard<py::gil_scoped_release>())
      .def("getIndependentVariableProfile",
           &PostProcessorCore::GetIndependentVariableProfile,
           py::call_guard<py::gil_scoped_release>())
      .def("additionalProfile", &PostProcessorCore::AdditionalProfile, py::arg("key"),
           py::call_guard<py::gil_scoped_release>());

  py::class_<ROPA>(m, "ROPA")
      .def(py::init<>())
      .def("setResults", &ROPA::SetResults, py::arg("data"),
           py::call_guard<py::gil_scoped_release>())
      .def("rateOfProductionAnalysis", &ROPA::RateOfProductionAnalysis,
           py::arg("number_of_reactions") = 10,
           py::arg("heterogeneous_reactions") = false,
           py::call_guard<py::gil_scoped_release>())
      .def("RateOfProductionAnalysis2D", &ROPA::RateOfProductionAnalysis2D,
           py::call_guard<py::gil_scoped_release>())  // TODO keyword arguments
      .def("fluxAnalysis", &ROPA::FluxAnalysis,
           py::call_guard<py::gil_scoped_release>())  // No arguments
      .def("getReactionRates", &ROPA::GetReactionRates, py::arg("reaction_indices"),
           py::arg("sum_rates") = false, py::arg("heterogeneous_reactions") = false,
           py::call_guard<py::gil_scoped_release>())
      .def("getFormationRates", &ROPA::GetFormationRates, py::arg("specie"),
           py::arg("units"), py::arg("type"), py::arg("heterogeneous_reactions") = false,
           py::call_guard<py::gil_scoped_release>())
      .def("setROPAType", &ROPA::SetROPAType, py::arg("type") = "global",
           py::call_guard<py::gil_scoped_release>())
      .def("setSpecies", &ROPA::SetSpecies, py::arg("species"),
           py::call_guard<py::gil_scoped_release>())
      .def("setLocalValue", &ROPA::SetLocalValue, py::arg("localValue"),
           py::call_guard<py::gil_scoped_release>())
      .def("setLowerBound", &ROPA::SetLowerBound, py::arg("lowerBound"),
           py::call_guard<py::gil_scoped_release>())
      .def("setUpperBound", &ROPA::SetUpperBound, py::arg("upperBound"),
           py::call_guard<py::gil_scoped_release>())
      // TODO keyword args for flux analysis
      .def("setElement", &ROPA::SetElement, py::call_guard<py::gil_scoped_release>())
      .def("setThickness", &ROPA::SetThickness, py::call_guard<py::gil_scoped_release>())
      .def("setFluxAnalysisType", &ROPA::SetFluxAnalysisType,
           py::call_guard<py::gil_scoped_release>())
      .def("setWidth", &ROPA::SetWidth, py::call_guard<py::gil_scoped_release>())
      .def("setDepth", &ROPA::SetDepth, py::call_guard<py::gil_scoped_release>())
      .def("setThreshold", &ROPA::SetThreshold, py::call_guard<py::gil_scoped_release>())
      .def("setThicknessLogScale", &ROPA::SetThicknessLogScale,
           py::call_guard<py::gil_scoped_release>())
      .def("setLabelType", &ROPA::SetLabelType, py::call_guard<py::gil_scoped_release>())

      .def("reactions", &ROPA::reactions, py::call_guard<py::gil_scoped_release>())
      .def("coefficients", &ROPA::coefficients, py::call_guard<py::gil_scoped_release>())
      .def("indexFirstName", &ROPA::indexFirstName,
           py::call_guard<py::gil_scoped_release>())
      .def("indexSecondName", &ROPA::indexSecondName,
           py::call_guard<py::gil_scoped_release>())
      .def("computedThickness", &ROPA::computedThickness,
           py::call_guard<py::gil_scoped_release>())
      .def("computedLabel", &ROPA::computedLabel,
           py::call_guard<py::gil_scoped_release>())
      .def("formationRates", &ROPA::formationRates,
           py::call_guard<py::gil_scoped_release>())
      .def("reactionRates", &ROPA::reactionRates,
           py::call_guard<py::gil_scoped_release>())
      .def("sumOfRates", &ROPA::sumOfRates, py::call_guard<py::gil_scoped_release>());

  py::enum_<Phase>(m, "Phase").value("Gas", Phase::Gas).value("Surface", Phase::Surface);

  py::class_<SensitivityCalculator>(m, "Sensitivity")
      .def(py::init<>())
      .def("setResults", &SensitivityCalculator::SetResults, py::arg("data"),
           py::call_guard<py::gil_scoped_release>())
      .def("setNormalizationType", &SensitivityCalculator::SetNormalizationType,
           py::arg("normalizationType") = "max-value",
           py::call_guard<py::gil_scoped_release>())
      .def("setSensitivityType", &SensitivityCalculator::SetSensitivityType,
           py::arg("sensitivityType") = "global",
           py::call_guard<py::gil_scoped_release>())
      .def("setOrderingType", &SensitivityCalculator::SetOrderingType,
           py::arg("orderingType") = "peak-values",
           py::call_guard<py::gil_scoped_release>())
      .def("setTarget", &SensitivityCalculator::SetTarget, py::arg("target"),
           py::call_guard<py::gil_scoped_release>())
      .def("setLocalValue", &SensitivityCalculator::SetLocalValue, py::arg("localValue"),
           py::call_guard<py::gil_scoped_release>())
      .def("setLowerBound", &SensitivityCalculator::SetLowerBound, py::arg("lowerBound"),
           py::call_guard<py::gil_scoped_release>())
      .def("setUpperBound", &SensitivityCalculator::SetUpperBound, py::arg("upperBound"),
           py::call_guard<py::gil_scoped_release>())
      .def("prepare", &SensitivityCalculator::Prepare, py::arg("phase") = Phase::Gas,
           py::call_guard<py::gil_scoped_release>())
      .def("sensitivityAnalysis", &SensitivityCalculator::Sensitivity_Analysis,
           py::arg("number_of_reactions") = 10, py::call_guard<py::gil_scoped_release>())
      .def("readSensitivityCoefficients",
           &SensitivityCalculator::ReadSensitivityCoefficients,
           py::call_guard<py::gil_scoped_release>())
      .def("getSensitivityProfile", &SensitivityCalculator::GetSensitivityProfile,
           py::arg("reaction_index"), py::call_guard<py::gil_scoped_release>())
      .def("reactions", &SensitivityCalculator::reactions,
           py::call_guard<py::gil_scoped_release>())
      .def("sensitivityCoefficients", &SensitivityCalculator::sensitivityCoefficients,
           py::call_guard<py::gil_scoped_release>());

  py::class_<SpeciesClass>(m, "SpeciesClass")
      .def(py::init<>())
      .def("setResults", &SpeciesClass::SetResults, py::arg("data"),
           py::call_guard<py::gil_scoped_release>())
      .def("speciesClassesAvailable", &SpeciesClass::speciesClassesAvailable,
           py::call_guard<py::gil_scoped_release>())
      .def("elementalDistribution", &SpeciesClass::ElementalDistribution,
           py::arg("element"), py::arg("normalize") = false,
           py::call_guard<py::gil_scoped_release>())
      .def("elementMolesBySpecies", &SpeciesClass::ElementMolesBySpecies,
           py::arg("element"), py::call_guard<py::gil_scoped_release>())
      .def("speciesNames", &SpeciesClass::speciesNames,
           py::call_guard<py::gil_scoped_release>())
      .def("elementMolesBySpeciesMatrix", &SpeciesClass::elementMolesBySpecies,
           py::call_guard<py::gil_scoped_release>())
      .def("setFluxPerClass", &SpeciesClass::SetFluxPerClass,
           py::arg("flux_per_class") = true, py::call_guard<py::gil_scoped_release>())
      .def("setClassName", &SpeciesClass::SetClassName, py::arg("class_name") = "",
           py::call_guard<py::gil_scoped_release>())
      .def("setSpecies", &SpeciesClass::SetSpecies, py::arg("species") = "",
           py::call_guard<py::gil_scoped_release>())
      .def("setElement", &SpeciesClass::SetElement, py::arg("element"),
           py::call_guard<py::gil_scoped_release>())
      .def("setFluxAnalysisType", &SpeciesClass::SetFluxAnalysisType,
           py::arg("type") = "production", py::call_guard<py::gil_scoped_release>())
      .def("setLocalValue", &SpeciesClass::SetLocalValue, py::arg("localValue") = 0.,
           py::call_guard<py::gil_scoped_release>())
      .def("setThickness", &SpeciesClass::SetThickness, py::arg("thickness") = "relative",
           py::call_guard<py::gil_scoped_release>())
      .def("setLabelType", &SpeciesClass::SetLabelType, py::arg("type") = "relative",
           py::call_guard<py::gil_scoped_release>())
      .def("setThicknessLogScale", &SpeciesClass::SetThicknessLogScale,
           py::arg("thicknesslogscale") = true, py::call_guard<py::gil_scoped_release>())
      .def("setWidth", &SpeciesClass::SetWidth, py::arg("width") = 3,
           py::call_guard<py::gil_scoped_release>())
      .def("setDepth", &SpeciesClass::SetDepth, py::arg("depth") = 2,
           py::call_guard<py::gil_scoped_release>())
      .def("setThreshold", &SpeciesClass::SetThreshold, py::arg("threshold") = 0.1,
           py::call_guard<py::gil_scoped_release>())
      .def("setCarbonWeighted", &SpeciesClass::SetCarbonWeighted,
           py::arg("carbon_weighted") = true, py::call_guard<py::gil_scoped_release>())
      .def("setAutoPruneDiagonal", &SpeciesClass::SetAutoPruneDiagonal,
           py::arg("auto_prune_diagonal") = true,
           py::call_guard<py::gil_scoped_release>())
      .def("fluxByClass", &SpeciesClass::FluxByClass,
           py::call_guard<py::gil_scoped_release>())
      .def("classNames", &SpeciesClass::classNames,
           py::call_guard<py::gil_scoped_release>())
      .def("abscissa", &SpeciesClass::abscissa, py::call_guard<py::gil_scoped_release>())
      .def("elementalFractions", &SpeciesClass::elementalFractions,
           py::call_guard<py::gil_scoped_release>())
      .def("fluxMatrix", &SpeciesClass::fluxMatrix,
           py::call_guard<py::gil_scoped_release>())
      .def("indexFirstClass", &SpeciesClass::indexFirstClass,
           py::call_guard<py::gil_scoped_release>())
      .def("indexSecondClass", &SpeciesClass::indexSecondClass,
           py::call_guard<py::gil_scoped_release>())
      .def("computedThickness", &SpeciesClass::computedThickness,
           py::call_guard<py::gil_scoped_release>())
      .def("computedLabel", &SpeciesClass::computedLabel,
           py::call_guard<py::gil_scoped_release>());

  py::class_<ReactionClass>(m, "ReactionClass")
      .def(py::init<>())
      .def("setHeterogeneous", &ReactionClass::SetHeterogeneous,
           py::arg("heterogeneous") = false, py::call_guard<py::gil_scoped_release>())
      .def("setResults", &ReactionClass::SetResults, py::arg("data"),
           py::call_guard<py::gil_scoped_release>())
      .def("classesAvailable", &ReactionClass::classesAvailable,
           py::call_guard<py::gil_scoped_release>())
      .def("mainClass", &ReactionClass::mainClass,
           py::call_guard<py::gil_scoped_release>())
      .def("subClass", &ReactionClass::subClass, py::call_guard<py::gil_scoped_release>())
      .def("mergeDuplicates",  
          // C++ style parameters-by-reference needs Python translation
          [](const ReactionClass& self,
             const std::vector<std::vector<int>>& species_indices,
             const std::vector<std::vector<double>>& species_coefficients) {
            std::vector<int> representative_indices;
            std::vector<std::vector<double>> merged_coefficients;
            self.MergeDuplicates(species_indices, species_coefficients,
                                 representative_indices, merged_coefficients);
            return std::make_pair(representative_indices, merged_coefficients);
          },
          py::arg("species_indices"), py::arg("species_coefficients"),
          py::call_guard<py::gil_scoped_release>());

  // Raw <SootProperties> (PolimiSoot BIN properties) accessor. The block is parsed
  // by PostProcessorCore during loadKinetics; each getter returns one
  // per-bin vector (all sharing the bin order), e.g. pp.soot.dpp().
  py::class_<Soot>(m, "Soot")
      .def(py::init<>())
      .def("setResults", &Soot::SetResults, py::arg("data"),
           py::call_guard<py::gil_scoped_release>())
      .def("sootAvailable", &Soot::sootAvailable,
           py::call_guard<py::gil_scoped_release>())
      .def("numberOfBins", &Soot::numberOfBins, py::call_guard<py::gil_scoped_release>())
      .def("index", &Soot::index, py::call_guard<py::gil_scoped_release>())
      .def("section", &Soot::section, py::call_guard<py::gil_scoped_release>())
      .def("nc", &Soot::nc, py::call_guard<py::gil_scoped_release>())
      .def("nh", &Soot::nh, py::call_guard<py::gil_scoped_release>())
      .def("no", &Soot::no, py::call_guard<py::gil_scoped_release>())
      .def("htoc", &Soot::htoc, py::call_guard<py::gil_scoped_release>())
      .def("mw", &Soot::mw, py::call_guard<py::gil_scoped_release>())
      .def("density", &Soot::density, py::call_guard<py::gil_scoped_release>())
      .def("volume", &Soot::volume, py::call_guard<py::gil_scoped_release>())
      .def("mass", &Soot::mass, py::call_guard<py::gil_scoped_release>())
      .def("numpp", &Soot::numpp, py::call_guard<py::gil_scoped_release>())
      .def("dsph", &Soot::dsph, py::call_guard<py::gil_scoped_release>())
      .def("dcol", &Soot::dcol, py::call_guard<py::gil_scoped_release>())
      .def("dpp", &Soot::dpp, py::call_guard<py::gil_scoped_release>())
      .def("df", &Soot::df, py::call_guard<py::gil_scoped_release>())
      .def("psd", &Soot::ParticleSizeDistribution, py::arg("local_value"),
           py::arg("particle_type") = "all", py::arg("diameter_type") = "dmob",
           py::arg("min_section") = 1,
           py::arg("correlation_name") = "Kelesidis", py::arg("merge_tol") = 0.20,
           py::call_guard<py::gil_scoped_release>())
      .def("averagedProfile", &Soot::AveragedProfile, py::arg("property"),
           py::arg("min_section") = 5, py::arg("weighting") = "mass",
           py::call_guard<py::gil_scoped_release>())
      .def("ssa", &Soot::SpecificSurfaceArea, py::arg("min_section") = 5,
           py::call_guard<py::gil_scoped_release>());
}
