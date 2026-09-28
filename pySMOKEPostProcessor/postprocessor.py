import numpy as np
import pandas as pd

from .graph_writer import GraphWriter
from .soot_utilities.psd import compute_psd
from .pySMOKEPostProcessor import (
    ROPA,
    PostProcessorCore,
    Phase,
    Phase2Kind,
    Sensitivity,
    SpeciesClass,
    Soot,
)

# np.trapz was removed in numpy 2.0 (renamed np.trapezoid); np.trapezoid does not exist
# before numpy 2.0. environment.yml requires numpy>=1.20, so support both.
# An unrelated `conda install <package>` (es seaborn) can shadow the error by installing
#  a newer numpy
_trapz = getattr(np, "trapezoid", None) or np.trapz


class PostProcessor:
    """
    Main Class of the package, needed to call the C++ backend.

    Buildable from kinetics only, output only, or both - a method that needs a
    missing input raises a clear RuntimeError (see _require_kinetics/_require_output).

    Attributes:
        db: PostProcessorCore, the C++ orchestrator owning both the kinetics and
            the output state.
        kineticFolder: last kinetics folder passed to __init__ (None if never loaded).
        outputFolder: last output folder passed to __init__/updateOutput (None if
            never loaded).
        km: KineticMapReader_Gas (pybind object) once kinetics is loaded, else None.
        kmhet: whichever phase-2 kinetics reader is loaded - KineticMapReader_Surface/
            _Liquid/_Solid - else None.
        isHeterogeneous: True when the loaded mechanism has a *surface* phase
            specifically (existing call sites, e.g. RateOfProductionAnalysis_Surface,
            branch on this to mean "use the *_Surface widgets"); a Liquid/Solid
            phase-2 leaves this False. Use db.phase2Kind() for the general case.
            # TODO this should be extended/modified for Liquid/Solid
    """

    def __init__(
        self,
        kineticFolder: str = None,
        outputFolder: str = None,
        isHeterogeneous: bool = None,
    ) -> None:
        self.db = PostProcessorCore()
        self.kineticFolder = kineticFolder
        self.outputFolder = outputFolder

        self.km = None
        self.kmhet = None
        self.isHeterogeneous = bool(isHeterogeneous)
        self.soot = None
        self.dfSootProperties = None

        if kineticFolder is not None:
            # Reads kinetics (gas + auto-detected phase-2, e.g. surface) and creates kineticMap
            self.db.loadKinetics(kineticFolder)
            self.km = self.db.gasKinetics()
            if self.db.hasPhase2Kinetics():
                self.kmhet = self.db.phase2Kinetics()
                if isHeterogeneous is None:
                    self.isHeterogeneous = self.db.phase2Kind() == Phase2Kind.Surface

        if outputFolder is not None:
            # isHeterogeneous can only be auto-detected from the kinetics folder
            # (kinetics.surface.xml presence); for output-only construction (no
            # kineticFolder) pass it explicitly if the run is heterogeneous.
            self.db.readFileResults(outputFolder, self.isHeterogeneous)

        if self.db.hasKinetics() and self.db.hasOutput():
            # Raw <SootProperties> (PolimiSoot BIN properties). 
            # The optional block is parsed while reading the kinetics.xml; 
            # self.soot is always present once both inputs are loaded, 
            # self.soot.sootAvailable() is False when the mechanism has no soot bins. 
            # Access the per-bin arrays directly, e.g. self.soot.dpp().
            self.soot = Soot()
            self.soot.setResults(self.db)

            # For convenience, the BinProperties is reconstructed from the kinetics.xml
            # as a pandas.DataFrame (if available).
            # None when the mechanism has no soot bins.
            self.dfSootProperties = self._build_soot_properties_dataframe()

    def _require_kinetics(self, method_name: str) -> None:
        if self.km is None:
            raise RuntimeError(
                "{}() needs a kinetic mechanism, but this PostProcessor was built "
                "without one (kineticFolder=None). Build with a kineticFolder to use "
                "it.".format(method_name)
            )

    def _require_output(self, method_name: str) -> None:
        if not self.db.hasOutput():
            raise RuntimeError(
                "{}() needs simulation output, but this PostProcessor was built "
                "without one (outputFolder=None). Build with an outputFolder, or "
                "call updateOutput(...), to use it.".format(method_name)
            )

    def updateOutput(self, outputFolder: str) -> None:
        """
        Re-points this PostProcessor at a different Output.xml under the same
        mechanism, without re-parsing kinetics.
        To be used for parametric analyses or comparison case-to-case
        with a fixed mechanism.

        Requires an output to already be loaded (outputFolder was given to
        __init__, or a previous updateOutput call succeeded); use
        PostProcessor(kineticFolder, outputFolder) for the first load.
        """
        self._require_output("updateOutput")
        self.db.updateOutput(outputFolder)
        self.outputFolder = outputFolder

    def getSpeciesProfile(self, name: str, basis: str = "mass", as_dataframe: bool = False):
        """
        (independent_variable, profile) tuple for one species, 
        read from Output.xml - no kinetics.xml dependency.
        These methods are the only callables in case the pp was built without kinetics.
        The independent_variable is to be interpreted as time for batch reactors, space for
        flames, and it depends on input setup for plug flow reactors (@Length vs @ResidenceTime)

        Args:
            name: species name.
            basis: "mass" or "mole".
            as_dataframe: return a single-column DataFrame indexed by the
                independent variable instead of the (x, y) tuple.
        """
        self._require_output("getSpeciesProfile")
        x, y = self.db.getSpeciesProfile(name, basis)
        if as_dataframe:
            return pd.DataFrame({name: y}, index=np.array(x))
        return np.array(x), np.array(y)

    def getIndependentVariableProfile(self, as_dataframe: bool = False):
        """
        The independent_variable is to be interpreted as time for batch reactors, space for
        flames, and it depends on input setup for plug flow reactors (@Length vs @ResidenceTime)
        """
        self._require_output("getIndependentVariableProfile")
        x = np.array(self.db.getIndependentVariableProfile())
        if as_dataframe:
            return pd.DataFrame({"x": x})
        return x

    def _getAdditionalProfile(self, key: str, as_dataframe: bool = False):
        """Shared implementation behind the named profile getters below.
        This is a common interface to get the profiles of the quantities in the 
        <additional> XML leaf.
        """
        self._require_output(key)
        x = self.getIndependentVariableProfile()
        y = np.array(self.db.additionalProfile(key))
        if as_dataframe:
            return pd.DataFrame({key: y}, index=x)
        return x, y

    def getTemperatureProfile(self, as_dataframe: bool = False):
        return self._getAdditionalProfile("temperature", as_dataframe)

    def getPressureProfile(self, as_dataframe: bool = False):
        return self._getAdditionalProfile("pressure", as_dataframe)

    def getDensityProfile(self, as_dataframe: bool = False):
        return self._getAdditionalProfile("density", as_dataframe)

    def getViscosityProfile(self, as_dataframe: bool = False):
        return self._getAdditionalProfile("viscosity", as_dataframe)

    # LG this is likely from previous OpenSMOKE versions, at the moment I think
    #    no solver prints the YSoot in the additional columns.
    def getYSootProfile(self, as_dataframe: bool = False): 
        return self._getAdditionalProfile("YSoot", as_dataframe)


    def _getAveragedSootProfile(
        self,
        key: str,
        prop,
        as_dataframe: bool = False,
        min_section: int = 5,
        weighting: str = "mass",
    ):
        """Shared implementation behind the averaged soot getters: `prop` is a
        per-BIN property (one value per BIN, e.g. self.soot.htoc()), averaged along
        the abscissa over the BINs with Bin_section >= min_section. weighting="mass"
        weighs each BIN by its soot mass, "carbon" by its moles of carbon. NaN where
        there is no soot.
        Note: default weighting should be 'carbon' only for H/C ratio
        Args:
            key:

        """
        self._require_output(key)
        cols = self.soot.averagedProfile(prop, min_section, weighting)
        x = np.array(cols["abscissa"])
        y = np.array(cols["mean"])
        if as_dataframe:
            return pd.DataFrame(
                {key: y, "soot_mass[kg/m3]": np.array(cols["mass_kg_per_m3"])}, index=x
            )
        return x, y
    
    def getFvSootProfile(self, as_dataframe: bool = False, min_section: int = 5):
        """Soot volume fraction [-] along the abscissa.
        Needs a <SootProperties> block.

        Args:
            as_dataframe: return a DataFrame ("fv[-]", "soot_mass[kg/m3]",
                indexed by the independent variable) instead of the (x, fv)
                tuple.
            min_section: only BINs with Bin_section >= min_section count as soot.

        Returns:
            (x, fv) tuple of np.ndarray (independent variable, volume fraction),
            or a DataFrame (columns "fv[-]", "soot_mass[kg/m3]") when
            as_dataframe=True.
        """
        with np.errstate(divide="ignore"):
            specific_volume = 1.0 / np.array(self.soot.density())  # m3/kg
        df = self._getAveragedSootProfile(
            "1/rho_soot[m3/kg]", specific_volume, True, min_section, "mass"
        )
        mass = df["soot_mass[kg/m3]"].to_numpy()
        has_soot = mass > 0.0
        fv = np.zeros_like(mass)
        fv[has_soot] = df["1/rho_soot[m3/kg]"].to_numpy()[has_soot] * mass[has_soot]
        x = df.index.to_numpy()
        if as_dataframe:
            return pd.DataFrame({"fv[-]": fv, "soot_mass[kg/m3]": mass}, index=x)
        return x, fv



    def getSSASootProfile(self, as_dataframe: bool = False, min_section: int = 5):
        """Mean soot specific surface area [m2/kg] along the abscissa, mass-weighted
        over the BINs with Bin_section >= min_section (total soot surface / total
        soot mass), each BIN having SSA = pi*Dpp**2*numPP/Bin_mass (numPP <= 0
        counted as 1). NaN where there is no soot. Needs a <SootProperties> block.

        Args:
            as_dataframe: return a DataFrame ("SSA[m2/kg]", "soot_mass[kg/m3]",
                indexed by the independent variable) instead of the (x, y) tuple.
            min_section: only BINs with Bin_section >= min_section count as soot.

        Returns:
            (x, y) tuple of np.ndarray (independent variable, SSA), or a
            DataFrame (columns "SSA[m2/kg]", "soot_mass[kg/m3]") when
            as_dataframe=True.
        """
        return self._getAveragedSootProfile(
            "SSA[m2/kg]", self.soot.ssa(min_section), as_dataframe, min_section, "mass"
        )



    def getHtoCSootProfile(
        self, as_dataframe: bool = False, min_section: int = 5, weighting: str = "mass"
    ):
        """Mean soot H/C [-] along the abscissa over the BINs with
        Bin_section >= min_section. weighting="mass" is the mass-weighted mean of the
        per-BIN H/C; "carbon" gives exactly total H atoms / total C atoms. Needs a
        <SootProperties> block.

        Args:
            as_dataframe: return a DataFrame ("H/C[-]", "soot_mass[kg/m3]",
                indexed by the independent variable) instead of the (x, y) tuple.
            min_section: only BINs with Bin_section >= min_section count as soot.
            weighting: "mass" (mass-weighted mean of per-BIN H/C) or "carbon"
                (moles of H / moles of C).

        Returns:
            (x, y) tuple of np.ndarray (independent variable, H/C), or a
            DataFrame (columns "H/C[-]", "soot_mass[kg/m3]") when
            as_dataframe=True.
        """
        return self._getAveragedSootProfile(
            "H/C[-]", self.soot.htoc(), as_dataframe, min_section, weighting
        )

    def _build_soot_properties_dataframe(self):
        s = self.soot
        if not s.sootAvailable():
            return None
        species = self.km.speciesNames()
        names = [species[i] if 0 <= i < len(species) else None for i in s.index()]
        return pd.DataFrame(
            {
                "Bin_index": list(s.index()),
                "Bin_name": names,
                "nC[-]": list(s.nc()),
                "nH[-]": list(s.nh()),
                "nO[-]": list(s.no()),
                "H/C[-]": list(s.htoc()),
                "MW[kg/kmol]": list(s.mw()),
                "density[kg/m3]": list(s.density()),
                "mass[kg]": list(s.mass()),
                "volume[m3]": list(s.volume()),
                "Dsph[m]": list(s.dsph()),
                "Dcol[m]": list(s.dcol()),
                "Dpp[m]": list(s.dpp()),
                "Df[-]": list(s.df()),
                "numPP[#]": list(s.numpp()),
                "SSA[m2/kg]": list(s.ssa()),
                "Bin_section": list(s.section()),
            }
        )

    def SootPSD(
        self,
        local_value: float,
        particle_type: str = "all",
        diameter_type: str = "dmob",
        min_section: int = 1,
        correlation_name: str = "Kelesidis",
        merge_tol: float = 0.20,
    ):
        """Soot particle size distribution.
        Pre-requisite: the mechanism has to be compiled with SootProperties enabled.

        Args:
            local_value: coordinate at which the analysis is performed.
            particle_type: 
                "all" keeps every BIN with Bin_section >= min_section, each 
                    counted as one particle (BINliq are counted as numPP=1) 
                    -> PSD, dN/dlog10(d).
                "primary" -> the primary-particle size distribution (PPSD):
                    every BIN (single particles AND aggregates) contributes 
                    the primary particles it is made of, classified by that 
                    BIN's OWN dpp
                    Only diameter_type="dpp" is valid with this.
            diameter_type: "dmob" mobility diameter dm = f(Dpp, numPP) - see
                correlation_name, "dpp" primary-particle diameter, "dcol"
                collision diameter, "dva" volume-equivalent sphere diameter.
            min_section: only BINs with Bin_section >= min_section are kept.
                Defaults to 1 here (unlike getSSASootProfile/getHtoCSootProfile/
                getFvSootProfile's 5): a size distribution is meant to show the
                full population down to the smallest sections, not just the
                bulk soot-only ones.
            correlation_name: only used when diameter_type="dmob" - 
                Literature correlation used to describe the mobility diameter:
                "Kelesidis" -> dm = Dpp*numPP**0.45 (Kelesidis et al. 2017,
                    Carbon, 10.1016/j.carbon.2017.06.004).
                "Sorensen"  -> dm = Dpp*numPP**0.465 (Sorensen 2011, Aerosol
                    Sci. Technol., 10.1080/02786826.2011.560909).
                "Rissler"   -> dm = 0.794*Dpp*numPP**0.51 (Rissler et al. 2013,
                    Aerosol Sci. Technol., 10.1080/02786826.2013.791381).
            merge_tol: relative tolerance used to merge near-equal diameters
                into fixed sections (representative = geometric mean) before
                dN/dlog10(d) is formed.

        Returns:
            A DataFrame (<d>[nm], <d>_min[nm], <d>_max[nm], N[#/m3], n_bins,
            Dlog10(<d>[nm]), dN/dlog10(<d>[nm])[#/m3]); the gas state actually
            used (abscissa, T, P, rho, MW) plus particle_type, diameter_type,
            min_section, correlation_name, local_value are in df.attrs. For
            particle_type="primary", N[#/m3] is a primary-particle number
            density and is what should be plotted directly (as
            ``plot_distribution`` does)
        """
        return compute_psd(
            self,
            local_value,
            particle_type=particle_type,
            diameter_type=diameter_type,
            min_section=min_section,
            correlation_name=correlation_name,
            merge_tol=merge_tol,
        )

    def RateOfProductionAnalysis(
        self,
        species: str,
        ropa_type: str,
        local_value: float = 0,
        lower_value: float = 0,
        upper_value: float = 0,
        number_of_reactions: int = 10,
        two_dimensions: bool = False,
        region_location: dict = None,
        mass_ropa: bool = False,
        include_names: bool = True,
        heterogeneous_reactions: bool = False,
    ) -> dict:
        """
        Function that performs the [R]ate [O]f [P]roduction [A]nalysis. Covers both gas
        and heterogeneous reactions (currently tested for surface-only).
        Args:
            species: Name of the target species for the ROPA
            ropa_type: Type of ROPA to be performed available are: global | local | region
            local_value: Local value of the domain in where perform the ROPA
            lower_value: Lower value of the domain for the region ROPA
            upper_value: Upper value of the domain for the region ROPA
            number_of_reactions: Number Of Reactions to return after the ROPA
            two_dimensions: Activate the support for the ROPA for 2D/3D simulations (TODO: find a better name)
            region_location: 2D option
            mass_ropa: Return the ROPA coefficients in mass unit. Not supported together
                with heterogeneous_reactions=True.
            include_names: Resolve reaction_names (one ReactionNameFromIndex lookup per
                returned reaction). Callers that only need reaction_indices/coefficients
                (e.g. reaction-class flux aggregation, which re-derives names from its own
                per-mechanism table) can set this False to skip it - cheap per call, but it
                adds up across the many species x many-timesteps loops in
                script_utils.process_classes / surfacereactions_utilities.deposition_plot.
            heterogeneous_reactions: Perform the ROPA on the surface (deposition)
                reactions instead of gas reactions. Requires a heterogeneous mechanism
                (self.isHeterogeneous); not supported together with two_dimensions=True.

        Returns:
            A dictionary as the following one:
                ropa_results = {'coefficients': [...],
                                'reaction_names': [...],
                                'reaction_indices': [...]}
                Containing the ROPA coefficients, the reaction names (None if include_names
                is False) and the indices of the reactions (in the selected reaction phase,
                when heterogeneous_reactions=True).
        """

        if two_dimensions:
            if heterogeneous_reactions:
                raise ValueError("two_dimensions and heterogeneous_reactions cannot both be True")
            return self.RateOfProductionAnalysis2D(species, ropa_type, number_of_reactions, region_location, mass_ropa)

        widget = ROPA()
        widget.setResults(self.db)
        widget.setROPAType(ropa_type)
        widget.setSpecies(species)
        widget.setLocalValue(local_value)
        widget.setLowerBound(lower_value)
        widget.setUpperBound(upper_value)

        widget.rateOfProductionAnalysis(number_of_reactions, heterogeneous_reactions)

        reaction_indices = widget.reactions()
        ropa_coefficients = widget.coefficients()

        reaction_names = None
        if include_names:
            kinmap = self.kmhet if heterogeneous_reactions else self.km
            reaction_names = [kinmap.formattedReactionNameFromIndex(i) for i in reaction_indices]

        if mass_ropa:
            if heterogeneous_reactions:
                raise ValueError("mass_ropa is not supported with heterogeneous_reactions=True")
            ropa_coefficients = self.convert_tomass(ropa_coefficients, species)

        ropa_result = {
            "coefficients": ropa_coefficients,
            "reaction_names": reaction_names,
            "reaction_indices": reaction_indices,
        }

        return ropa_result

    def RateOfProductionAnalysis2D(
        self,
        species: str,
        ropa_type: str,
        number_of_reactions: int = 10,
        region_location: dict = None,
        mass_ropa: bool = False,
    ) -> dict:
        widget = ROPA()
        widget.setResults(self.db)
        widget.setROPAType(ropa_type)
        widget.setSpecies(species)
        widget.setLocalValue(0)
        widget.setLowerBound(0)
        widget.setUpperBound(0)

        local_value_x = region_location["local_value_x"]
        local_value_z = region_location["local_value_z"]
        lower_value_x = region_location["lower_value_x"]
        lower_value_z = region_location["lower_value_z"]
        upper_value_x = region_location["upper_value_x"]
        upper_value_z = region_location["upper_value_z"]

        widget.RateOfProductionAnalysis2D(
            number_of_reactions,
            local_value_x,
            local_value_z,
            lower_value_x,
            upper_value_x,
            lower_value_z,
            upper_value_z,
        )

        reaction_indices = widget.reactions()
        ropa_coefficients = widget.coefficients()

        reaction_names = []
        for i in reaction_indices:
            reaction_names.append(self.km.formattedReactionNameFromIndex(i))

        if mass_ropa:
            ropa_coefficients = self.convert_tomass(ropa_coefficients, species)

        ropa_result = {
            "coefficients": ropa_coefficients,
            "reaction_names": reaction_names,
            "reaction_indices": reaction_indices,
        }

        return ropa_result

    def SensitivityAnalysis(
        self,
        target: str,
        sensitivity_type: str,
        ordering_type: str = 'peak-values',
        normalization_type: str = 'max-value',
        local_value: float = 0,
        lower_value: float = 0,
        upper_value: float = 0,
        number_of_reactions: int = 10,
    ) -> dict:
        # SENSITIVITY HERE
        widget = Sensitivity()

        widget.setResults(self.db)
        widget.setSensitivityType(sensitivity_type)
        widget.setOrderingType(ordering_type)
        widget.setNormalizationType(normalization_type)
        widget.setTarget(target)
        widget.setLocalValue(local_value)
        widget.setLowerBound(lower_value)
        widget.setUpperBound(upper_value)
        widget.prepare()
        widget.readSensitivityCoefficients()
        widget.sensitivityAnalysis(number_of_reactions)

        reaction_indices = widget.reactions()
        sensitivity_coefficients = widget.sensitivityCoefficients()

        reaction_names = []
        for i in reaction_indices:
            reaction_names.append(self.km.formattedReactionNameFromIndex(i))

        sensitivity_result = {
            "coefficients": sensitivity_coefficients,
            "reaction_names": reaction_names,
            "reaction_indices": reaction_indices,
        }

        return sensitivity_result

    def SensitivityAnalysis_Surface(
        self,
        target: str,
        sensitivity_type: str = 'global',
        ordering_type: str = 'peak-values',
        normalization_type: str = 'max-value',
        local_value: float = 0,
        lower_value: float = 0,
        upper_value: float = 0,
        number_of_reactions: int = 10,
        heterogeneous_sensitivity: bool = False
    ) -> dict:
        widget = Sensitivity()

        widget.setResults(self.db)
        widget.setSensitivityType(sensitivity_type)
        widget.setOrderingType(ordering_type)
        widget.setNormalizationType(normalization_type)
        widget.setTarget(target)
        widget.setLocalValue(local_value)
        widget.setLowerBound(lower_value)
        widget.setUpperBound(upper_value)
        widget.prepare(Phase.Surface if heterogeneous_sensitivity else Phase.Gas)
        widget.readSensitivityCoefficients()
        widget.sensitivityAnalysis(number_of_reactions)

        reaction_indices = widget.reactions()
        sensitivity_coefficients = widget.sensitivityCoefficients()

        reaction_names = []
        for i in reaction_indices:
            kinmap = self.kmhet if heterogeneous_sensitivity else self.km
            reaction_names.append(kinmap.formattedReactionNameFromIndex(i))

        sensitivity_result = {
            "coefficients": sensitivity_coefficients,
            "reaction_names": reaction_names,
            "reaction_indices": reaction_indices,
        }

        return sensitivity_result

    def FluxAnalysis(
        self,
        species: str,
        element: str,
        flux_analysis_type: str,
        thickness: str,
        thickness_log_scale: bool,
        label_type: str,
        depth: int = 2,
        width: int = 5,
        threshold: float = 0,
        local_value: float = 0.01,
    ):
        widget = ROPA()

        widget.setResults(self.db)
        widget.setSpecies(species)
        widget.setElement(element)
        widget.setFluxAnalysisType(flux_analysis_type)
        widget.setLocalValue(local_value)
        widget.setThickness(thickness)
        widget.setThicknessLogScale(thickness_log_scale)
        widget.setLabelType(label_type)
        widget.setDepth(depth)
        widget.setWidth(width)
        widget.setThreshold(threshold)

        widget.fluxAnalysis()

        indexFirstName = widget.indexFirstName()
        indexSecondName = widget.indexSecondName()
        computedThickness = widget.computedThickness()
        computedLabel = widget.computedLabel()

        firstNames = []
        secondNames = []
        for i, j in enumerate(indexFirstName):
            firstNames.append(self.km.speciesNameFromIndex(j))
            secondNames.append(self.km.speciesNameFromIndex(indexSecondName[i]))
        Graph = GraphWriter(flux_analysis_type)  # , species, element)
        Graph = Graph.CreateGraph(firstNames, secondNames, computedThickness, computedLabel)

        return Graph

    def _species_class_widget(self):
        """
        Builds the C++ SpeciesClass widget.
        The <SpeciesClasses> block inside kinetics.xml is optional,
        This function is an intermediate to raise a clear python error
        instead of core-dumping as soon as the block is not found.
        """
        widget = SpeciesClass()
        widget.setResults(self.db)
        if not widget.speciesClassesAvailable():
            raise Exception(
                "The kinetic mechanism does not contain a <SpeciesClasses> block."
                "Species-class post-processing is not available for this mechanism."
            )
        return widget

    def ElementalDistributionByClass(self, element: str, normalize: bool = True) -> pd.DataFrame:
        """
        Distribution of an atomic element across the species classes along the
        independent variable (e.g. time for a batch reactor).

        Args:
            element: element symbol as written in the mechanism (e.g. "C", "H", "O").
            normalize: if True every abscissa column sums to 1 (fraction of the
                element carried by each class); if False the values are moles of
                the element per unit mass of mixture. Moles conserve over all phases,
                but it does not necessarily conserve for single phases (e.g. surface depo)

        Returns:
            DataFrame indexed by the independent variable, one column per class.
        """
        widget = self._species_class_widget()
        widget.elementalDistribution(element, normalize)

        class_names = widget.classNames()
        fractions = np.array(widget.elementalFractions(), dtype=np.float64)  # [class][point]
        return pd.DataFrame(fractions.T, index=np.array(widget.abscissa()), columns=class_names)

    def ElementMolesBySpecies(self, element: str) -> pd.DataFrame:
        """
        Per-species (not per-class) moles-of-element-per-unit-mass-of-mixture along
        the independent variable - every species in the mechanism, classified or
        not. Shares its C++ numeric core with ElementalDistributionByClass (which
        aggregates this same computation by class); this is the per-species view
        elements_balance needs for its lumping/threshold logic. Does not require a
        <SpeciesClasses> block.

        Returns:
            DataFrame indexed by the independent variable, one column per species.
        """
        widget = SpeciesClass()
        widget.setResults(self.db)
        widget.elementMolesBySpecies(element)

        species_names = widget.speciesNames()
        moles = np.array(widget.elementMolesBySpeciesMatrix(), dtype=np.float64)  # [species][point]
        return pd.DataFrame(moles.T, index=np.array(widget.abscissa()), columns=species_names)

    def FluxAnalysisByClass(
        self,
        element: str,
        flux_analysis_type: str,
        class_name: str = None,
        species_name: str = None,
        flux_per_class: bool = True,
        thickness: str = "relative",
        thickness_log_scale: bool = True,
        label_type: str = "relative",
        depth: int = 2,
        width: int = 3,
        threshold: float = 0.02,
        local_value: float = 0.01,
        carbon_weighted: bool = True,
        auto_prune_diagonal: bool = True,
    ) -> dict:
        """
        Element flux between species classes at a point of the independent
        variable.

        The full directed class -> class element-flux matrix is computed from
        every reaction (no species-level pruning). ``depth``, ``width`` and
        ``threshold`` only prune the *graph*, which is walked breadth-first from
        a seed class: ``depth`` class generations, at most ``width`` links per
        class node, each link at least ``threshold`` % of that node's flux.
        ``flux_analysis_type`` picks the walk direction ("destruction" = where
        the element goes, "production" = where it comes from).

        auto_prune_diagonal (default True):
            Intra-class flux (a class's flux to itself) is typically far larger
            than any cross-class one and makes a heatmap of ``matrix``
            unreadable; the graph never draws it either way (self-loops are
            excluded from the walk). With the default, the diagonal of
            ``matrix`` is zeroed too, so the returned matrix is heatmap-ready.
            Pass False to get the raw diagonal back (e.g. to inspect intra-class
            activity numerically).

        Seed selection -- pass exactly one of:
            flux_per_class=True  (default): ``class_name`` -- walk from that class.
            flux_per_class=False:           ``species_name`` -- walk from the
                                            class that contains that species.

        carbon_weighted (default True):
            True  -- OpenSMOKE's element-flux weighting n_i*n_j*|R_r| / N_C,r
                (maps/FluxAnalysisMap.hpp::AnalyzeNetFluxes): the actual
                carbon-atom throughput. Faithful; note (bin carbon count)^2 makes
                soot-soot terms dominate a mechanism with lumped BINs.
            False -- each reaction contributes its own molar rate |R_r|, split
                between (source class, target class) pairs by the carbon share of
                each species. Bounded per reaction, so lumped soot/BIN classes
                stay comparable to small gas species.

        Returns:
            {"matrix": DataFrame, rows = source class, columns = target class,
                       full directed element flux (diagonal zeroed unless
                       auto_prune_diagonal=False);
             "graph":  graphviz.Digraph of the pruned class walk}
        """
        widget = self._species_class_widget()
        class_names = widget.classNames()

        if flux_per_class:
            if class_name is None or species_name is not None:
                raise ValueError(
                    "flux_per_class=True requires 'class_name' and no 'species_name'"
                )
            if class_name not in class_names:
                raise ValueError(
                    "class_name '{}' is not a species class in this mechanism; available: {}".format(
                        class_name, ", ".join(class_names)
                    )
                )
        else:
            if species_name is None or class_name is not None:
                raise ValueError(
                    "flux_per_class=False requires 'species_name' and no 'class_name'"
                )
            if species_name not in self.km.speciesNames():
                raise ValueError(
                    "species_name '{}' is not in the mechanism".format(species_name)
                )

        widget.setFluxPerClass(flux_per_class)
        widget.setClassName(class_name or "")
        widget.setSpecies(species_name or "")
        widget.setElement(element)
        widget.setFluxAnalysisType(flux_analysis_type)
        widget.setThickness(thickness)
        widget.setThicknessLogScale(thickness_log_scale)
        widget.setLabelType(label_type)
        widget.setDepth(depth)
        widget.setWidth(width)
        widget.setThreshold(threshold)
        widget.setLocalValue(local_value)
        widget.setCarbonWeighted(carbon_weighted)
        widget.setAutoPruneDiagonal(auto_prune_diagonal)
        widget.fluxByClass()

        matrix = pd.DataFrame(
            np.array(widget.fluxMatrix(), dtype=np.float64), index=class_names, columns=class_names
        )

        first_names = [class_names[a] for a in widget.indexFirstClass()]
        second_names = [class_names[b] for b in widget.indexSecondClass()]
        graph = GraphWriter(flux_analysis_type).CreateGraph(
            first_names, second_names, widget.computedThickness(), widget.computedLabel()
        )
        return {"graph": graph, "matrix": matrix}

    def GetReactionRates(self, reaction_name: list = None, reaction_index: list = None, sum_rates: bool = False, heterogeneous_reactions = False):
        if reaction_name is not None:
            if not heterogeneous_reactions: # If homogeneous, it will be false anyway
                reaction_index = [self.km.reactionIndexFromName(name=i) for i in reaction_name]
            else:
                reaction_index = [self.kmhet.reactionIndexFromName(name=i) for i in reaction_name]

        widget = ROPA()
        widget.setResults(self.db)
        widget.getReactionRates(reaction_index, sum_rates, heterogeneous_reactions)

        if sum_rates:
            reaction_rates = [widget.sumOfRates()]
        else:
            reaction_rates = widget.reactionRates()

        return reaction_rates

    def GetFormationRates(
        self,
        formation_rate_type: str,
        species: str,
        units: str = "mole",
        heterogeneous_reactions: bool = False,
    ):
        widget = ROPA()
        widget.setResults(self.db)
        widget.getFormationRates(species, units, formation_rate_type, heterogeneous_reactions)
        formationRates = widget.formationRates()

        return formationRates

    def SensitivityCoefficients(
        self,
        target: str,
        normalization_type: str,
        reaction_name: str = None,
        reaction_index: int = None,
    ):
        if reaction_name is not None:
            reaction_index = self.km.reactionIndexFromName(name=reaction_name)

        widget = Sensitivity()
        widget.setResults(self.db)
        widget.setSensitivityType("global")
        widget.setOrderingType("peak-values")
        widget.setNormalizationType(normalization_type)
        widget.setTarget(target)
        widget.setLocalValue(0.0)
        widget.setLowerBound(0.0)
        widget.setUpperBound(0.0)
        widget.prepare()
        widget.readSensitivityCoefficients()
        widget.getSensitivityProfile(reaction_index)
        sensitivity_coefficients = widget.sensitivityCoefficients()

        return sensitivity_coefficients

    def convert_tomass(self, ropa_coefficients: list, species: str) -> list:
        """
        Function that converts the rate of production coefficients from mole
        units to mass units.
        Args:
            ropa_coefficients: list containing the rate of production
            coefficients
            species: species name
        Returns:
            list containing the rate of production coefficients in mass unit.
        """
        sp_idx = self.km.speciesIndexFromName(species)
        mwi = self.km.mw(sp_idx)
        ropa_coefficients = [c * mwi for c in ropa_coefficients]

        return ropa_coefficients

    def reactionrategroups(self, rxnnames_sr, xaxis: list, threshold: float = 0.01, heterogeneous_reactions = False):
        # reaction rates by groups (example: by class)
        # rxnnames_sr: series with labels and reaction names
        # xaxis
        # threshold: plot only if contributes above threshold% to the total rate
        rr = dict.fromkeys(rxnnames_sr.index)
        rrsum = pd.Series(index=rxnnames_sr.index, dtype=np.float64)
        for label, rxnnames in rxnnames_sr.items():
            rr[label] = np.array(self.GetReactionRates(
                reaction_name=rxnnames, sum_rates=True,
                heterogeneous_reactions=heterogeneous_reactions)[0])
            rrsum[label] = _trapz(y=rr[label], x=xaxis)

        # check cumulative contribution and filter based on threshold
        rrsum /= np.sum(abs(rrsum))  # abs?
        filteredlabels = list(rrsum[abs(rrsum) > threshold].index)
        filtered_rr = [rr[key] for key in filteredlabels]
        rates_df = pd.DataFrame(np.array(filtered_rr).T, columns=filteredlabels, index=xaxis)

        return rates_df

    def cumulativerates(
        self,
        xaxis: list,
        ropa_dct: dict,
        rate_type: str = "PC",
        threshold: float = 0.01,
        heterogeneous_reactions: bool = False
    ):
        """""
        cumulative reaction rate matrix extract
        rate_type: 'PC' = net, 'P' = production, 'C' = consumption
        threshold: delete rates based on % contribution (default: keep only those contributing > 1%)
        xaxis: derived from output
        """""
        # 0. ropa DCT: sum coefficients for duplicates
        coefficients, indices, names, names_split = [], [], [], []
        allindices_array = np.array(ropa_dct["reaction_indices"])
        allcoeffs_array = np.array(ropa_dct["coefficients"])
        for i, idx in enumerate(allindices_array):
            if idx not in indices:
                indices.append(idx)
                names.append(ropa_dct["reaction_names"][i])
                names_split.append(ropa_dct["reaction_names"][i].split(": ")[1])
                if ropa_dct["reaction_indices"].count(idx) > 1:
                    positions = np.where((allindices_array == idx))[0]
                    coefficients.append(np.sum(allcoeffs_array[positions]))
                else:
                    coefficients.append(allcoeffs_array[i])

        # 1. ropa DF: indexes (positive or negative) and reaction names
        factor = [1.0 - 2.0 * (coeff < 0) for coeff in coefficients]
        ropa_df = pd.DataFrame(
            np.array([factor, names, names_split], dtype=object).T,
            index=indices,
            columns=["factor", "reaction_names", "reaction_names_split"],
        )

        if rate_type == "P":
            ropa_df = ropa_df[ropa_df["factor"].values > 0]
        elif rate_type == "C":
            ropa_df = ropa_df[ropa_df["factor"] < 0]

        # 2. get reaction rates
        # if cumulative rate and ropa sign agree: don't change sign; otherwise, do
        rr = dict.fromkeys(ropa_df.index)
        rrsum = pd.Series(index=ropa_df.index, dtype=np.float64)

        rrdel = []
        scannedidxs = []
        for idx in ropa_df.index:
            check = False

            rr_idx = np.array(self.GetReactionRates(
                reaction_index=[idx], heterogeneous_reactions=heterogeneous_reactions)[0])

            rrsum_idx = _trapz(y=rr_idx, x=xaxis)
            if (rrsum_idx * float(ropa_df["factor"][idx])) < 0:
                # integral and ropa have opposite signs: change sign
                rr_idx *= -1
                rrsum_idx *= -1

            # if the reaction corresponds to the bw rxn of an irrev rxn already treated: sum
            # look for bw rxn
            if "=>" in ropa_df["reaction_names_split"][idx]:
                prod, reac = ropa_df["reaction_names_split"][idx].split("=>")
                bwname = "=>".join([reac, prod])
                idxrs = ropa_df.index[ropa_df["reaction_names_split"] == bwname].tolist()
                for idxr in idxrs:
                    if idxr in scannedidxs:
                        rr[idxr] += rr_idx
                        rrsum[idxr] += rrsum_idx
                        ropa_df.loc[idxr, "reaction_names"] = ropa_df["reaction_names"][idx].replace("=>", "=")
                        rrdel.append(idx)
                        check = True

            # if a reaction with the same name was already analyzed: add flux to that
            rxnname = ropa_df["reaction_names_split"][idx]
            idxrs = ropa_df.index[ropa_df["reaction_names_split"] == rxnname].tolist()
            idxr_scanned = [idxr for idxr in idxrs if idxr in scannedidxs]
            if len(idxr_scanned) > 0:
                idxr = idxr_scanned[0]
                rr[idxr] += rr_idx
                rrsum[idxr] += rrsum_idx
                ropa_df.loc[idxr, "reaction_names"] = ropa_df["reaction_names"][idx].replace("=>", "=")
                rrdel.append(idx)
                check = True

            if check == False:
                # assign index
                rr[idx] = rr_idx
                rrsum[idx] = rrsum_idx
                scannedidxs.append(idx)
            elif check == True:
                del rr[idx]

        # drop indexes of bw rxn
        rrsum = rrsum.drop(labels=rrdel)
        # check cumulative contribution and delete based on threshold
        rrsum /= np.sum(abs(rrsum))  # abs?
        filteredidxs = list(rrsum[abs(rrsum) > threshold].index)
        filtered_rr = [rr[key] for key in filteredidxs]

        # add to the names the cumulative % contribution
        names_wpct = []
        for idx in filteredidxs:
            name = ropa_df["reaction_names_split"][idx]
            pct = rrsum[idx] * 100
            names_wpct.append("{} {:.5f}%".format(name, pct))
        cumulativerates_df = pd.DataFrame(np.array(filtered_rr).T, columns=names_wpct, index=xaxis)

        return cumulativerates_df


# alternative : do a local ropa to set values; if you don't find the reaction,
# set the value to zero
