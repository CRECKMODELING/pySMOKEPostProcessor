import numpy as np
import pandas as pd
import os

from .graph_writer import GraphWriter
from .maps.KineticMap import KineticMap, KineticMapSurface
from .soot_utilities.psd import compute_psd
from .pySMOKEPostProcessor import ROPA, ProfilesDatabase, Sensitivity, ROPA_Surface, Sensitivity_Surface, SpeciesClass, Soot

# np.trapz was removed in numpy 2.0 (renamed np.trapezoid); np.trapezoid does not exist
# before numpy 2.0. environment.yml only pins numpy>=1.20, so support both - and don't
# assume the env stays on whatever numpy pip/conda happened to resolve today: an
# unrelated `conda install <package>` can pull in a newer numpy as a dependency and
# silently drop trapz (this is what actually happened here - see conda-meta/history).
_trapz = getattr(np, "trapezoid", None) or np.trapz


class PostProcessor:
    """
    Main Class of the package, needed to call the C++ backend
    Attributes:
        db: Object representing the Output.xml file inside the C++ code
        kineticFolder: Path pointing to the folder containing the kinetic mechanism.
        outputFolder: Path pointing to the folder containing the Simulations Output.
        km: Internal object representing the kinetic map. (TODO: replace in favour of OpenSMOKE_Interfaces)
    """

    def __init__(self, kineticFolder: str, outputFolder: str) -> None:
        self.db = ProfilesDatabase()
        self.kineticFolder = kineticFolder
        self.outputFolder = outputFolder

        #   Check if surface kinetics is available
        surface_file = os.path.join(self.kineticFolder, "kinetics.surface.xml")
        self.isHeterogeneous = os.path.exists(surface_file)

        # Reads Output
        self.db.readFileResults(self.outputFolder, self.isHeterogeneous)
        
        # Reads kinetics and creates kineticMap
        self.db.readKineticMechanism(self.kineticFolder)
        self.km = KineticMap(self.kineticFolder)

        # Reads and creates kineticMapSurface (only if available)
        if self.isHeterogeneous:
            self.db.readHeterogeneousKineticMechanism(self.kineticFolder,"Surface")
            self.kms = KineticMapSurface(self.kineticFolder)

        # Raw <SootProperties> (PolimiSoot BIN properties). The optional block is
        # parsed by while reading the kinetics.xml; self.soot is always
        # present, self.soot.sootAvailable() is False when the mechanism has no soot
        # bins. Access the per-bin arrays directly, e.g. self.soot.dpp().
        self.soot = Soot()
        self.soot.setDataBase(self.db)

        # One-row-per-bin view of the whole <SootProperties> block. Column names
        # match OpenSMOKE's BinProperties.txt so the two are directly comparable;
        # Bin_index is the bin's species index in the gas-phase scheme, Bin_name
        # its resolved species name. None when the mechanism has no soot bins.
        self.dfSootProperties = self._build_soot_properties_dataframe()

    def _build_soot_properties_dataframe(self):
        s = self.soot
        if not s.sootAvailable():
            return None
        species = self.km.species
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
                "Bin_section": list(s.section()),
            }
        )

    def SootPSD(
        self,
        local_value: float = 0.0,
        particle_type: str = "all",
        diameter_type: str = "dmob",
        min_section: int = 5,
        mobility_exponent: float = 0.45,
        merge_tol: float = 0.20,
    ):
        """Soot particle size distribution.

        particle_type: "all" keeps every BIN with Bin_section >= min_section
            (numPP <= 0 nascent bins counted as one spherule); "primary" keeps only
            the free primary particles (numPP == 1) - the classic PPSD - and
            ignores min_section.
        diameter_type: "dmob" mobility diameter dm = Dpp*numPP**mobility_exponent,
            "dpp" primary-particle diameter, "dcol" collision diameter.
        local_value picks the profile point like local ROPA (first point whose
            abscissa >= local_value; abscissa is time for a reactor, a coordinate
            for a flame).

        Returns a DataFrame (<d>[nm], N[#/m3], dN/dlog10(<d>[nm])[#/m3], ...); the
        gas state used is in df.attrs. Needs a <SootProperties> block.
        """
        return compute_psd(
            self,
            local_value,
            particle_type=particle_type,
            diameter_type=diameter_type,
            min_section=min_section,
            mobility_exponent=mobility_exponent,
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
    ) -> dict:
        """
        Function that performs the [R]ate [O]f [P]roduction [A]nalysis
        Args:
            species: Name of the target species for the ROPA
            ropa_type: Type of ROPA to be performed available are: global | local | region
            local_value: Local value of the domain in where perform the ROPA
            lower_value: Lower value of the domain for the region ROPA
            upper_value: Upper value of the domain for the region ROPA
            number_of_reactions: Number Of Reactions to return after the ROPA
            two_dimensions: Activate the support for the ROPA for 2D/3D simulations (TODO: find a better name)
            region_location: 2D option
            mass_ropa: Return the ROPA coefficients in mass unit
            include_names: Resolve reaction_names (one ReactionNameFromIndex lookup per
                returned reaction). Callers that only need reaction_indices/coefficients
                (e.g. reaction-class flux aggregation, which re-derives names from its own
                per-mechanism table) can set this False to skip it - cheap per call, but it
                adds up across the many species x many-timesteps loops in
                script_utils.process_classes / surfacereactions_utilities.deposition_plot.

        Returns:
            A dictionary as the following one:
                ropa_results = {'coefficients': [...],
                                'reaction_names': [...],
                                'reaction_indices': [...]}
                Containing the ROPA coefficients, the reaction names (None if include_names
                is False) and the indices of the reactions.
        """

        if two_dimensions is False:
            widget = ROPA()
            widget.setDataBase(self.db)
            widget.setROPAType(ropa_type)
            widget.setSpecies(species)
            widget.setLocalValue(local_value)
            widget.setLowerBound(lower_value)
            widget.setUpperBound(upper_value)

            widget.rateOfProductionAnalysis(number_of_reactions)

            reaction_indices = widget.reactions()
            ropa_coefficients = widget.coefficients()

            reaction_names = None
            if include_names:
                reaction_names = [self.km.ReactionNameFromIndex(i) for i in reaction_indices]

            if mass_ropa:
                ropa_coefficients = self.convert_tomass(ropa_coefficients, species)

            ropa_result = {
                "coefficients": ropa_coefficients,
                "reaction_names": reaction_names,
                "reaction_indices": reaction_indices,
            }

            return ropa_result
        else:
            return self.RateOfProductionAnalysis2D(species, ropa_type, number_of_reactions, region_location, mass_ropa)

    def RateOfProductionAnalysis_Surface(
        self,
        species: str,
        ropa_type: str,
        local_value: float = 0,
        lower_value: float = 0,
        upper_value: float = 0,
        number_of_reactions: int = 10,
        heterogeneous_reactions: bool = False,
        #   mass_ropa: bool = False     # Not sure we need the conversion to mass. I see the point, but for now not important
        include_names: bool = True,
    ) -> dict:
        """
        Function that performs the [R]ate [O]f [P]roduction [A]nalysis for the coupled gas mechanism with deposition (Surface)
        Args:
            species: Name of the target species for the ROPA
            ropa_type: Type of ROPA to be performed available are: gloabal | local | region
            local_value: Local value of the domain in where perform the ROPA
            lower_value: Lower value of the domain for the region ROPA
            upper_value: Upper value of the domain for the region ROPA
            number_of_reactions: Number Of Reactions to return after the ROPA
            heterogeneous_reactions: Decide if ROPA is performed on Gas (False) or Surface (True) reactions
            include_names: Resolve reaction_names (one ReactionNameFromIndex lookup per
                returned reaction). See RateOfProductionAnalysis for when to skip it.

        Returns:
            A dictionary as the following one:
                ropa_results = {'coefficients': [...],
                                'reaction_names': [...],
                                'reaction_indices': [...]}
                Containing the ROPA coefficients, the reaction names (None if include_names
                is False) and the indices of the reactions in the selected reaction phase
        """

        widget = ROPA_Surface()
        widget.setDataBase(self.db)
        widget.setROPAType(ropa_type)
        widget.setROPAPhase(heterogeneous_reactions)
        widget.setSpecies(species)
        widget.setLocalValue(local_value)
        widget.setLowerBound(lower_value)
        widget.setUpperBound(upper_value)

        widget.rateOfProductionAnalysis(number_of_reactions,heterogeneous_reactions)

        reaction_indices = widget.reactions()
        ropa_coefficients = widget.coefficients()
        reaction_names = None
        if include_names:
            name_lookup = self.kms.ReactionNameFromIndex if heterogeneous_reactions else self.km.ReactionNameFromIndex
            reaction_names = [name_lookup(i) for i in reaction_indices]

        # if mass_ropa:
        #     ropa_coefficients = self.convert_tomass(ropa_coefficients, species)

        ropa_result = {
            "coefficients": ropa_coefficients,
            "reaction_names": reaction_names,
            "reaction_indices": reaction_indices
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
        widget.setDataBase(self.db)
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
            reaction_names.append(self.km.ReactionNameFromIndex(i))

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

        widget.setDataBase(self.db)
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
            reaction_names.append(self.km.ReactionNameFromIndex(i))

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
        widget = Sensitivity_Surface()

        widget.setDataBase(self.db)
        widget.setSensitivityType(sensitivity_type)
        widget.setOrderingType(ordering_type)
        widget.setNormalizationType(normalization_type)
        widget.setTarget(target)
        widget.setLocalValue(local_value)
        widget.setLowerBound(lower_value)
        widget.setUpperBound(upper_value)
        widget.prepare(heterogeneous_sensitivity)
        widget.readSensitivityCoefficients()
        widget.sensitivityAnalysis(number_of_reactions)

        reaction_indices = widget.reactions()
        sensitivity_coefficients = widget.sensitivityCoefficients()

        reaction_names = []
        for i in reaction_indices:
            if heterogeneous_sensitivity:
                reaction_names.append(self.kms.ReactionNameFromIndex(i))
            else:
                reaction_names.append(self.km.ReactionNameFromIndex(i))

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

        widget.setDataBase(self.db)
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
            firstNames.append(self.km.SpeciesNameFromIndex(j))
            secondNames.append(self.km.SpeciesNameFromIndex(indexSecondName[i]))
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
        widget.setDataBase(self.db)
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
            if species_name not in self.km.species:
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
                reaction_index = [self.km.ReactionIndexFromName(name=i) for i in reaction_name]
            else:
                reaction_index = [self.kms.ReactionIndexFromName(name=i) for i in reaction_name]

        if not self.isHeterogeneous:
            widget = ROPA()
            widget.setDataBase(self.db)
            widget.getReactionRates(reaction_index, sum_rates)
        else:
            widget = ROPA_Surface()
            widget.setDataBase(self.db)
            widget.getReactionRates(reaction_index, sum_rates,heterogeneous_reactions)

        if sum_rates:
            reaction_rates = [widget.sumOfRates()]
        else:
            reaction_rates = widget.reactionRates()

        return reaction_rates

    def GetFormationRates(self, formation_rate_type: str, species: str, units: str = "mole"):
        # TODO @lgiardini implement this in the hpp and update for surface reactions
        widget = ROPA()
        widget.setDataBase(self.db)
        widget.getFormationRates(species, units, formation_rate_type)
        formationRates = widget.formationRates()

        return formationRates

    # [LG] This function seems unused and redundant with respect to function within the sensitivity class in C++
    #       plus it is not updated for heterogeneous stuff.
    def SensitivityCoefficients(
        self,
        target: str,
        normalization_type: str,
        reaction_name: str = None,
        reaction_index: int = None,
    ):
        if reaction_name is not None:
            reaction_index = self.km.ReactionIndexFromName(name=reaction_name)

        widget = Sensitivity()
        widget.setDataBase(self.db)
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
        sp_idx = self.km.IndexFromSpeciesName(species)
        mwi = self.km.mws[sp_idx]
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
            if not self.isHeterogeneous:
                rr[label] = np.array(self.GetReactionRates(reaction_name=rxnnames, sum_rates=True)[0])
            else:
                rr[label] = np.array(self.GetReactionRates(reaction_name=rxnnames, sum_rates=True,heterogeneous_reactions=heterogeneous_reactions)[0])
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

            if not self.isHeterogeneous:
                rr_idx = np.array(self.GetReactionRates(reaction_index=[idx])[0])
            else:
                rr_idx = np.array(self.GetReactionRates(reaction_index=[idx],heterogeneous_reactions=heterogeneous_reactions)[0])

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
