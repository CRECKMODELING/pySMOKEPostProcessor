from pySMOKEPostProcessor.surfacereactions_utilities.time_history_utilities import getTimeProfile, getROPAIntegralTimeHistory
from pySMOKEPostProcessor import script_utils
from pySMOKEPostProcessor.postprocessor import PostProcessor
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


class CumulativeROPAByClass:
    """
    Shared parent class above CumulativeDeposition and CumulativeSootProduction: 
    both run a local (per-timestep) ROPA-by-class over the whole simulation 
    for a fixed set of target species, collapse those species onto one combined 
    row per timestep, time-integrate the result, and stack-plot it by reaction class. 
    What actually differs between the two is domain-specific:
     - which species/phase to run the ROPA on
     - how to convert the integrated moles/mass into plotted units
    """

    def __init__(self,
                kineticFolder: str,
                outputFolder: str,
                class_group_file: str,
                targets: list,
                heterogeneous_reactions: bool,
                sort_type = [['reactiontype']],
                mass_ropa: bool = False,
                pp: PostProcessor = None):

        self.kineticFolder = kineticFolder
        self.outputFolder = outputFolder
        self.targets = targets
        self.heterogeneous_reactions = heterogeneous_reactions
        self.mass_ropa = mass_ropa
        self.sort_type = sort_type    # Default classification option
        # Built once and reused for every timestep. A subclass that already needs a
        # PostProcessor before calling super().__init__ (e.g. to resolve its own
        # target species off pp.km) can hand it over here instead of paying for a
        # second Output.xml/kinetics.xml parse.
        self.pp = pp if pp is not None else PostProcessor(kineticFolder, outputFolder)
        self.rxns_sorted = script_utils.get_sortedrxns(self.pp, class_group_file, heterogeneous_reactions=heterogeneous_reactions)

    def setClassification(self, sort_type):
        self.sort_type = sort_type

    def _cumulativeROPAIntegral(self, profilesTag: str, lump_steps: int) -> pd.DataFrame:
        """
        Runs a local ROPA-by-class at every (lumped) timestep and collapses every
        target species onto one combined row per timestep

        Returns the time-integrated result (see getROPAIntegralTimeHistory);
        stashes the un-integrated local values on self.ROPAbyClass_t, as before.
        """
        # The tag argument is only useful for phases where it has a different name (Riccardo...)
        timesteps = getTimeProfile(outputFolder=self.outputFolder, tag=profilesTag)

        collapsed = []
        for i in range(int(len(timesteps)/lump_steps)): # process_classes(pp=self.pp) reuses the
                            # PostProcessor/classification built once in __init__ instead of
                            # re-parsing Output.xml/kinetics.xml on every timestep
            i_lump = i*lump_steps
            dfs = script_utils.process_classes(
                    simul_fld=self.outputFolder, kin_xml_fld=self.kineticFolder, rxns_sorted_obj=self.rxns_sorted,
                    species_list=self.targets,   sortlists=self.sort_type,       heterogeneous_reactions=self.heterogeneous_reactions,
                    ropa_type='local',           local_value=timesteps[i_lump],  weigh='false', mass_ropa=self.mass_ropa,
                    pp=self.pp)
            for df in dfs:
                # Adding "time" (float64) onto df before summing promotes the summed
                # Series to float64 too - df's other columns are float32 (see
                # reaction_classes.py), so summing them alone would keep the Series
                # float32 and silently truncate the timestep assigned into it below.
                df["time"] = timesteps[i_lump]
                s = df.fillna(0).sum(axis=0)
                s["time"] = timesteps[i_lump]
                collapsed.append(s)

        df_ROPAt = pd.DataFrame(collapsed).fillna(0)
        self.ROPAbyClass_t = df_ROPAt
        return getROPAIntegralTimeHistory(df_ROPAt)

    def _plotStacked(self,
                     df_ROPAintegral: pd.DataFrame,
                     ylab: str,
                     threshold: float,
                     fig: plt.Figure,
                     ax: plt.Axes,
                     colors) -> tuple:
        cols_unfiltered = [c for c in df_ROPAintegral.columns if c != "time"]
        total = df_ROPAintegral[cols_unfiltered].sum(axis=1)
        cols = [col for col in cols_unfiltered if (abs(df_ROPAintegral[col].max()) > threshold*total.max())]

        if fig is None or ax is None:
            fig,ax = plt.subplots(figsize=(7.,6.))
        ax.stackplot( df_ROPAintegral["time"], *[df_ROPAintegral[col] for col in cols], labels=[col for col in cols],
                    alpha=0.87, edgecolor = 'k', linewidth = 0.6, colors = colors[:len(cols)] )
        fontsz = 14
        ax.set_xlim(left=0,right=df_ROPAintegral['time'].values[-1])
        ax.set_xlabel("Time [s]",fontsize=fontsz)
        ax.set_ylabel(ylabel=ylab,fontsize=fontsz)
        ax.grid(True,alpha=0.2)
        ax.legend(fontsize=fontsz,loc='center left', bbox_to_anchor=(1, 0.5))
        fig.tight_layout()
        return fig,ax


class CumulativeDeposition(CumulativeROPAByClass):
    """
    Parameters
    ----------
    kineticFolder : string
                Pointer to the folder containing the kinetics.surface.xml and surface_reaction_names.xml files
    outputFolder:   string
                Pointer to the folder containing the Output.xml file

    allCarbon:      bool
                Decide if deposition is evaluated on all bulk species or on C(B) only.

    """
    def __init__(self,
                kineticFolder: str,
                outputFolder: str,
                class_group_file: str,
                allCarbon: bool = False,
                sort_type = [['reactiontype']],
                carbon_density: float = 2000.,
                area: float = None):

        self.bulk_density = carbon_density  # Carbon density [kg/m3]
        self.MW = 12.01     # Carbon MW [kg/kmol]
        self.area = area    # Hp constant area, will be eventually read from Output.xml when variable
        self.iExtensiveEnabled = self.area is not None

        targets = ['C(B)']
        self.allCarbon = allCarbon
        if allCarbon:
            targets.append('c(B)')  # Non mi piace nemmeno hardcoded, va fatto meglio

        super().__init__(kineticFolder, outputFolder, class_group_file, targets,
                         heterogeneous_reactions=True, sort_type=sort_type)

    def plotCumulativeDeposition(self,
                    profilesTag: str = 'profiles',
                    lump_steps: int = 1,
                    units: str = "mass-specific",
                    threshold: float = 0.01,
                    fig: plt.Figure = None,
                    ax: plt.Axes = None,
                    colors = sns.color_palette("muted") ):

        """
        Parameters
        ----------
            profileTag: string
                Name of the tag to look for in the Output.xml file. Standard is 'profiles'.
            lump_steps: int
                If different from 1, reduces the number of steps used for integral ROPA evaluation, Nsteps = N/lump_steps
            units: str
                Units used as y-axis in the plot. Alternatives: 'mass-specific', 'moles-specific', 'mass', 'moles', 'thickness'
                For extensive units (e.g. mass), it is required to also set the area in the class constructor.
            threshold: float
                Threshold for filtering reaction classes to plot. Only classes whose maximum contribution is above the (relative) threshold are plotted.
            fig,ax: matplotlib.pyplot objects
                If not provided, they are created by the function. Intended for subplots usage.

        Returns
        -------
        figure,axis (matplotlib.pyplot objects) for the figure generated by the function
        """
        df_ROPAintegral = self._cumulativeROPAIntegral(profilesTag, lump_steps)
        df_ROPAintegral = self.convertToDesiredUnits(data = df_ROPAintegral, units = units)
        ylab = self.getLabelFromUnits(units)

        return self._plotStacked(df_ROPAintegral, ylab, threshold, fig, ax, colors)

    def convertToDesiredUnits(self, data:pd.DataFrame, units:str):
        # TODO: Add specific-volumetric quantities (for soot comparison)
        if not self.iExtensiveEnabled:
            if units == 'mass' or units == 'moles':
                raise ValueError("Error! Extensive units require area as input!")
        cols = [c for c in data.columns if c != "time"]
        match units:
            case 'mass-specific':
                data[cols] *= self.MW
                return data
            case 'moles':
                data[cols] *= self.area
            case 'mass':
                data[cols] *= self.area * self.MW
            case 'mass-specific-volumetric':
                raise ValueError("'specific-volumetric' units not implemented yet!")
            case 'moles-specific-volumetric':
                raise ValueError("'specific-volumetric' units not implemented yet!")
            case 'thickness':
                data[cols] *= self.MW/self.bulk_density * 1.E6
            case 'moles-specific':
                dummy = True
            case _:
                raise ValueError("Units not recognised. Options available: 'mass-specific', 'moles-specific', 'mass', 'moles', 'thickness' ")

        return data

    @staticmethod
    def getLabelFromUnits(units):
        ylab = "Deposited carbon "
        match units:
            case 'mass-specific':
                ylab += '[kg/m2]'
            case 'moles-specific':
                dummy = True
            case 'moles':
                ylab += '[kmol]'
            case 'mass':
                ylab += '[kg]'
            case 'mass-specific-volumetric':
                raise ValueError("'specific-volumetric' units not implemented yet!")
            case 'moles-specific-volumetric':
                raise ValueError("'specific-volumetric' units not implemented yet!")
            case 'thickness':
                ylab += r'[$\mu$m]'
            case _:
                raise ValueError("Units not recognised")
        return ylab


# Before Nobili, there was no <SpeciesClasses> block at all, so legacy BIN species
# (e.g. Kik monodisperse model) are hardcoded.
#  Kept only as a fallback for those older mechanisms.
# New soot mechanisms compiled from modules already use the SpeciesClass to classify soot
# for postProcessing as well (e.g. @SootProperties requires the SpeciesClass keyword to work)
# Mechanisms compiled not-from modules don't have this, so the soot species names have
# to be passed as parameters.
_LEGACY_SOOT_BIN_SPECIES = ["BIN5AJ", "BIN5BJ", "BIN5CJ", "BIN6AJ", "BIN6BJ", "BIN6CJ", "BIN7AJ", "BIN7BJ", "BIN7CJ", "BIN8AJ", "BIN8BJ", "BIN8CJ", "BIN9AJ", "BIN9BJ", "BIN9CJ", "BIN10AJ", "BIN10BJ", "BIN10CJ", "BIN11AJ", "BIN11BJ", "BIN11CJ", "BIN12AJ", "BIN12BJ", "BIN12CJ", "BIN13AJ", "BIN13BJ", "BIN13CJ", "BIN14AJ", "BIN14BJ", "BIN14CJ", "BIN15AJ", "BIN15BJ", "BIN15CJ", "BIN16AJ", "BIN16BJ", "BIN16CJ", "BIN17AJ", "BIN17BJ", "BIN17CJ", "BIN18AJ", "BIN18BJ", "BIN18CJ", "BIN19AJ", "BIN19BJ", "BIN19CJ", "BIN20AJ", "BIN20BJ", "BIN20CJ", "BIN21AJ", "BIN21BJ", "BIN21CJ", "BIN22AJ", "BIN22BJ", "BIN22CJ", "BIN23AJ", "BIN23BJ", "BIN23CJ", "BIN24AJ", "BIN24BJ", "BIN24CJ", "BIN25AJ", "BIN25BJ", "BIN25CJ"]

class CumulativeSootProduction(CumulativeROPAByClass):
    def __init__(self,
            kineticFolder: str,
            outputFolder: str,
            class_group_file: str,
            soot_species: list[str] = None,
            sort_type = [['reactiontype']],
            volume: float = None,
            ):

        self.volume = volume   # Hopefully constant.
        self.iExtensiveEnabled = self.volume is not None

        # Needed to resolve soot_species (off pp.km) before the rest of __init__ can
        # run - handed to super().__init__ below so it doesn't get built twice.
        pp = PostProcessor(kineticFolder, outputFolder)
        if soot_species is None:
            soot_species = self._resolveSootSpecies(pp)

        super().__init__(kineticFolder, outputFolder, class_group_file, soot_species,
                         heterogeneous_reactions=False, sort_type=sort_type, mass_ropa=True, pp=pp)

    @staticmethod
    def _resolveSootSpecies(pp: PostProcessor) -> list:
        """
        Figures out which species make up soot when soot_species isn't given explicitly.

        Preferred (mechanisms from 2024 on, module-compiled): 
        the mechanism's own <SpeciesClasses> block groups the different kinds of 
        soot particles (primary particles, liquid-like, aggregates...) 
        into one or more classes named SOOT* - use whichever of those are
        actually present, so this tracks the mechanism instead of a hardcoded list.

        Fallback (pre-2024 mechanisms, no <SpeciesClasses> block at all): match
        _LEGACY_SOOT_BIN_SPECIES against this mechanism's actual species and keep
        whichever are found.

        Raises:
            ValueError: neither approach found anything - pass soot_species explicitly.
        """
        try:
            pp.km.SpeciesClasses()
            soot_classes = [c for c in pp.km.speciesClasses if c.startswith("SOOT")]
            if not soot_classes:
                raise ValueError("no SOOT* entry in <SpeciesClasses>")
            species = [sp for c in soot_classes for sp in pp.km.speciesClasses[c]]
            return species
        except Exception:
            pass

        matched = [sp for sp in _LEGACY_SOOT_BIN_SPECIES if sp in pp.km.species]
        if matched:
            if len(matched) < len(_LEGACY_SOOT_BIN_SPECIES):
                print(' * Warning: only {}/{} legacy soot BIN species names found in this '
                      'mechanism'.format(len(matched), len(_LEGACY_SOOT_BIN_SPECIES)))
            return matched

        raise ValueError(
            "Could not determine soot species automatically: this mechanism has no "
            "SOOT* entry in <SpeciesClasses>, and none of the legacy hardcoded BIN "
            "species names (pre-2024 mechanisms) were found either. "
            "Pass soot_species explicitly."
        )

    def plotSootProduction(self,
                    profilesTag: str = 'profiles',
                    lump_steps: int = 1,
                    units: str = "mass-specific",
                    threshold: float = 0.01,
                    fig: plt.Figure = None,
                    ax: plt.Axes = None,
                    colors = sns.color_palette("muted") ):
        """
        Parameters
        ----------
            profileTag: string
                Name of the tag to look for in the Output.xml file. Standard is 'profiles'.
            lump_steps: int
                If different from 1, reduces the number of steps used for integral ROPA evaluation, Nsteps = N/lump_steps
            units: str
                Units used as y-axis in the plot. Alternatives: 'mass-specific', 'moles-specific', 'mass', 'moles', 'thickness'
                For extensive units (e.g. mass), it is required to also set the area in the class constructor.
            threshold: float
                Threshold for filtering reaction classes to plot. Only classes whose maximum contribution is above the (relative) threshold are plotted.
            fig,ax: matplotlib.pyplot objects
                If not provided, they are created by the function. Intended for subplots usage.

        Returns
        -------
        figure,axis (matplotlib.pyplot objects) for the figure generated by the function
        """
        df_ROPAintegral = self._cumulativeROPAIntegral(profilesTag, lump_steps)
        df_ROPAintegral = self.convertToDesiredUnits(data = df_ROPAintegral, units = units)
        ylab = self.getLabelFromUnits(units)

        return self._plotStacked(df_ROPAintegral, ylab, threshold, fig, ax, colors)

    def convertToDesiredUnits(self, data:pd.DataFrame, units:str):
        if not self.iExtensiveEnabled:
            if units == 'mass':
                raise ValueError("Error! Extensive units require volume as input!")
        cols = [c for c in data.columns if c != "time"]
        match units:
            case 'mass-specific':
                dummy = True
            case 'mass':
                data[cols] *= self.volume
            case _:
                raise ValueError("Units not recognised. Options available: 'mass-specific', 'mass'")
        return data

    @staticmethod
    def getLabelFromUnits(units):
        ylab = "Soot "
        match units:
            case 'mass-specific':
                ylab += '[kg/m3]'
            case 'mass':
                ylab += '[kg]'
            case _:
                raise ValueError("Units not recognised")
        return ylab
