#include "PostProcessorWrapper_py.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "source/ProfilesDatabase.h"
#include "source/ROPA.h"
#include "source/ROPA_Surface.h"
#include "source/Sensitivities.h"
#include "source/Sensitivities_Surface.h"
#include "source/ReactionClasses.h"
#include "source/SpeciesClasses.h"
#include "source/Soot.h"

namespace py = pybind11;
constexpr auto byref = py::return_value_policy::reference_internal;

PYBIND11_MODULE(pySMOKEPostProcessor, m) {
  m.doc() = "Python Interface to the OpenSMOKEpp Graphical Post Processor";

  py::class_<ProfilesDatabase>(m, "ProfilesDatabase")
      .def(py::init<>())
      .def("readKineticMechanism", &ProfilesDatabase::ReadKineticMechanism,
           py::arg("folder_name") = "kinetics", py::call_guard<py::gil_scoped_release>())
      .def("readHeterogeneousKineticMechanism",
           &ProfilesDatabase::ReadHeterogeneousKineticMechanism,
           py::arg("folder_name") = "kinetics", py::arg("phase_name") = "Surface",
           py::call_guard<py::gil_scoped_release>())
      .def("readFileResults", &ProfilesDatabase::ReadFileResults,
           py::arg("folder_name") = "Output",
           py::arg("isHeterogeneous") = false,
           py::call_guard<py::gil_scoped_release>())
      .def("prepare", &ProfilesDatabase::Prepare,
           py::call_guard<py::gil_scoped_release>())
      .def("prepareHeterogeneous", &ProfilesDatabase::PrepareHeterogeneous,
           py::call_guard<py::gil_scoped_release>());

  py::class_<ROPA>(m, "ROPA")
      .def(py::init<>())
      .def("setDataBase", &ROPA::SetDatabase, py::arg("data"),
           py::call_guard<py::gil_scoped_release>())
      .def("rateOfProductionAnalysis", &ROPA::RateOfProductionAnalysis,
           py::arg("number_of_reactions") = 10, py::call_guard<py::gil_scoped_release>())
      .def("ropa", &ROPA::RateOfProductionAnalysis2D,
           py::call_guard<py::gil_scoped_release>())  // TODO keyword arguments
      .def("fluxAnalysis", &ROPA::FluxAnalysis,
           py::call_guard<py::gil_scoped_release>())  // No arguments
      .def("getReactionRates", &ROPA::GetReactionRates, py::arg("reaction_indices"),
           py::arg("sum_rates") = false, py::call_guard<py::gil_scoped_release>())
      .def("getFormationRates", &ROPA::GetFormationRates, py::arg("specie"),
           py::arg("units"), py::arg("type"), py::call_guard<py::gil_scoped_release>())
      .def("setKineticFolder", &ROPA::SetKineticFolder,
           py::arg("kineticFolder") = "kinetics",
           py::call_guard<py::gil_scoped_release>())
      .def("setOutputFolder", &ROPA::SetOutputFolder, py::arg("outputFolder") = "Output",
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

  py::class_<ROPA_Surface, ROPA>(m, "ROPA_Surface")
      .def(py::init<>())
      .def("rateOfProductionAnalysis", &ROPA_Surface::RateOfProductionAnalysis,
           py::arg("number_of_reactions") = 10,
           py::arg("heterogeneous_reactions") = false,
           py::call_guard<py::gil_scoped_release>())
      .def("getReactionRates", &ROPA_Surface::GetReactionRates,
           py::arg("reaction_indices"), py::arg("sum_rates") = false,
           py::arg("heterogeneous_reactions") = false,
           py::call_guard<py::gil_scoped_release>())
      .def("getFormationRates", &ROPA_Surface::GetFormationRates, py::arg("specie"),
           py::arg("units") = "mole", py::arg("type") = "net",
           py::arg("heterogeneous_reactions") = false,
           py::call_guard<py::gil_scoped_release>())
      .def("setROPAPhase", &ROPA_Surface::SetROPAPhase,
           py::arg("heterogeneous_reactions") = false,
           py::call_guard<py::gil_scoped_release>());

  py::class_<Sensitivities>(m, "Sensitivity")
      .def(py::init<>())
      .def("setDataBase", &Sensitivities::SetDatabase, py::arg("data"),
           py::call_guard<py::gil_scoped_release>())
      .def("setNormalizationType", &Sensitivities::SetNormalizationType,
           py::arg("normalizationType") = "max-value",
           py::call_guard<py::gil_scoped_release>())
      .def("setSensitivityType", &Sensitivities::SetSensitivityType,
           py::arg("sensitivityType") = "global",
           py::call_guard<py::gil_scoped_release>())
      .def("setOrderingType", &Sensitivities::SetOrderingType,
           py::arg("orderingType") = "peak-values",
           py::call_guard<py::gil_scoped_release>())
      .def("setTarget", &Sensitivities::SetTarget, py::arg("target"),
           py::call_guard<py::gil_scoped_release>())
      .def("setLocalValue", &Sensitivities::SetLocalValue, py::arg("localValue"),
           py::call_guard<py::gil_scoped_release>())
      .def("setLowerBound", &Sensitivities::SetLowerBound, py::arg("lowerBound"),
           py::call_guard<py::gil_scoped_release>())
      .def("setUpperBound", &Sensitivities::SetUpperBound, py::arg("upperBound"),
           py::call_guard<py::gil_scoped_release>())
      .def("prepare", &Sensitivities::Prepare, py::call_guard<py::gil_scoped_release>())
      .def("sensitivityAnalysis", &Sensitivities::Sensitivity_Analysis,
           py::arg("number_of_reactions") = 10, py::call_guard<py::gil_scoped_release>())
      .def("readSensitivityCoefficients", &Sensitivities::ReadSensitivityCoefficients,
           py::call_guard<py::gil_scoped_release>())
      .def("getSensitivityProfile", &Sensitivities::GetSensitivityProfile,
           py::arg("reaction_index"), py::call_guard<py::gil_scoped_release>())
      .def("reactions", &Sensitivities::reactions,
           py::call_guard<py::gil_scoped_release>())
      .def("sensitivityCoefficients", &Sensitivities::sensitivityCoefficients,
           py::call_guard<py::gil_scoped_release>());

  // sensitivityAnalysis/readSensitivityCoefficients are inherited unchanged from
  // Sensitivities (see Sensitivities_Surface.h) - only Prepare() differs.
  py::class_<Sensitivities_Surface, Sensitivities>(m, "Sensitivity_Surface")
      .def(py::init<>())
      .def("prepare", &Sensitivities_Surface::Prepare,
           py::arg("heterogeneousSensitivity") = false,
           py::call_guard<py::gil_scoped_release>())
      .def("sensitivityAnalysis", &Sensitivities::Sensitivity_Analysis,
           py::arg("number_of_reactions") = 10, py::call_guard<py::gil_scoped_release>())
      .def("readSensitivityCoefficients", &Sensitivities::ReadSensitivityCoefficients,
           py::call_guard<py::gil_scoped_release>());

  py::class_<SpeciesClass>(m, "SpeciesClass")
      .def(py::init<>())
      .def("setDataBase", &SpeciesClass::SetDatabase, py::arg("data"),
           py::call_guard<py::gil_scoped_release>())
      .def("speciesClassesAvailable", &SpeciesClass::speciesClassesAvailable,
           py::call_guard<py::gil_scoped_release>())
      .def("elementalDistribution", &SpeciesClass::ElementalDistribution,
           py::arg("element"), py::arg("normalize") = false,
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
      .def("speciesToClass", &SpeciesClass::speciesToClass,
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
      .def("setDataBase", &ReactionClass::SetDatabase, py::arg("data"),
           py::call_guard<py::gil_scoped_release>())
      .def("classesAvailable", &ReactionClass::classesAvailable,
           py::call_guard<py::gil_scoped_release>())
      .def("mainClass", &ReactionClass::mainClass,
           py::call_guard<py::gil_scoped_release>())
      .def("subClass", &ReactionClass::subClass, py::call_guard<py::gil_scoped_release>())
      .def("mergeDuplicates", // C++ style parameters-by-reference needs Python translation
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
  // by ProfilesDatabase during readKineticMechanism; each getter returns one
  // per-bin vector (all sharing the bin order), e.g. pp.soot.dpp().
  py::class_<Soot>(m, "Soot")
      .def(py::init<>())
      .def("setDataBase", &Soot::SetDatabase, py::arg("data"),
           py::call_guard<py::gil_scoped_release>())
      .def("sootAvailable", &Soot::sootAvailable,
           py::call_guard<py::gil_scoped_release>())
      .def("numberOfBins", &Soot::numberOfBins,
           py::call_guard<py::gil_scoped_release>())
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
      .def("psd", &Soot::ParticleSizeDistribution, py::arg("local_value") = 0.0,
           py::arg("particle_type") = "all", py::arg("diameter_type") = "dmob",
           py::arg("min_section") = 5, py::arg("mobility_exponent") = 0.45,
           py::arg("merge_tol") = 0.20, py::call_guard<py::gil_scoped_release>());
}
