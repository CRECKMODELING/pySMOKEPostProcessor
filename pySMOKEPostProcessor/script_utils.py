"""
wrapper functions calling multiple functionalities
"""

from .maps.OpenSMOKEppXMLFile import OpenSMOKEppXMLFile
from .postprocessor import PostProcessor
from .reaction_classes import FluxByClass, assignclass
from .reaction_classes_utilities.reaction_classes_calc import filter_class0, sortby0


def get_sortedrxns(pp: PostProcessor, class_groups_file, heterogeneous_reactions: bool = False):
    """
    return dataframe of sorted reaction classes based on the specified class groups file

    Args:
        pp: a PostProcessor for the mechanism/simulation pair to classify. Classification
            itself only depends on the mechanism, but reading it requires the same
            ProfilesDatabase every other analysis on `pp` already uses (no output-folder
            re-parse per call: build `pp` once, reuse it here and for the ROPA calls that
            follow, e.g. via process_classes(..., pp=pp)).
        class_groups_file: path to the plain-text class-groups file.
        heterogeneous_reactions: classify the surface mechanism instead of the gas one.
    """
    rxns_sorted_obj, _ = assignclass(pp, class_groups_file, heterogeneous_reactions)

    return rxns_sorted_obj


# generic function for processing


def process_classes(
    simul_fld,
    kin_xml_fld,
    rxns_sorted_obj,
    species_list,
    sortlists,
    ropa_type,
    n_of_rxns: int = 100,
    filter_dcts=None,
    threshs=None,
    weigh="normbyspecies",
    local_value: float =0.0,
    upper_value: float =0.0,
    lower_value: float =0.0,
    mass_ropa: bool = False,
    heterogeneous_reactions: bool = False,
    pp: PostProcessor = None,
):
    sortdfs = []

    if filter_dcts is None:
        filter_dcts = [{}] * len(sortlists)
    if threshs is None:
        threshs = [1e-3] * len(sortlists)

    # If a pre-existing pp is passed, no reason to re-build it.
    if pp is None:
        pp = PostProcessor(kin_xml_fld, simul_fld)
    # initialize
    fluxbyclass = FluxByClass(rxns_sorted_obj, verbose=False)

    # ROPA for each species - if species_list contains dictionary, extract flux for each
    if isinstance(species_list[0], dict):
        flat_species_list = sum([list(d.values())[0] for d in species_list], [])
    else:
        flat_species_list = species_list

    tot_rop_dct = dict.fromkeys(flat_species_list)
    if not pp.isHeterogeneous:
        for sp in flat_species_list:
            tot_rop_dct[sp] = pp.RateOfProductionAnalysis(
                sp,
                ropa_type,
                local_value=local_value,
                lower_value=lower_value,
                upper_value=upper_value,
                number_of_reactions=n_of_rxns,
                mass_ropa=mass_ropa,
                include_names=False )
    else:
        for sp in flat_species_list:
            tot_rop_dct[sp] = pp.RateOfProductionAnalysis_Surface(
                sp,
                ropa_type,
                local_value=local_value,
                lower_value=lower_value,
                upper_value=upper_value,
                number_of_reactions=n_of_rxns,
                #mass_ropa=mass_ropa,
                heterogeneous_reactions=heterogeneous_reactions,
                include_names=False )

    # assign flux and process according to selected criteria
    fluxbyclass.process_flux(
        species_list,
        tot_rop_dct,
    )
    for i, sortlist in enumerate(sortlists):
        # assign flux and process according to selected criteria
        sortdf = fluxbyclass.sort_and_filter(
            sortlist,
            filter_dct=filter_dcts[i],
            thresh=threshs[i],
            weigh=weigh,
            dropunsorted=False,
        )

        sortdfs.append(sortdf)

    return sortdfs


# cumulative rates:
# dictionary with dataframes of cumulative reaction rates for species
# then ready to plot


def cumulative_rates(
    simul_fld,
    kin_xml_fld,
    species_list,
    rate_type,
    x_axis,
    n_of_rxns=100,
    mass_ropa=False,
    threshold=0.01,
    heterogeneous_reactions = False,
    pp=None,
):
    # pp -- for ropa
    if pp is None:
        pp = PostProcessor(kin_xml_fld, simul_fld)
    # output - for x axis
    output = OpenSMOKEppXMLFile(simul_fld, kin_xml_fld)
    try:
        x = getattr(output, x_axis)
    except AttributeError:
        print(' * Warning: attribute {} not found. using "time" as default'.format(x_axis))
        x = output.time  # x coordinate
    # ROPA for each species - if species_list contains dictionary, extract flux for each
    cum_df_dct = dict.fromkeys(species_list)
    for species in species_list:
        if not pp.isHeterogeneous:
            tot_rop_dct = pp.RateOfProductionAnalysis(
                species, ropa_type="global", number_of_reactions=n_of_rxns, mass_ropa=mass_ropa )
        else:
            tot_rop_dct = pp.RateOfProductionAnalysis_Surface(
                species, ropa_type="global", number_of_reactions=n_of_rxns, heterogeneous_reactions=heterogeneous_reactions )
        cum_df_dct[species] = pp.cumulativerates(x, tot_rop_dct, rate_type=rate_type, threshold=threshold)

    return cum_df_dct


# reaction rates by class / subclass / rxn type


def reactionrates_byclasses(
    simul_fld,
    kin_xml_fld,
    rxns_sorted_df,
    sortlists,
    x_axis,
    filter_by_species=[],
    filter_dcts=None,
    threshs=None,
    mass_ropa=False,
    heterogeneous_reactions = False,
    pp=None,
):
    sortdfs = []

    if filter_dcts is None:
        filter_dcts = [{}] * len(sortlists)
    if threshs is None:
        threshs = [1e-3] * len(sortlists)

    # pp -- for reactionrates
    if pp is None:
        pp = PostProcessor(kin_xml_fld, simul_fld)
    # output - for x axis
    output = OpenSMOKEppXMLFile(simul_fld, kin_xml_fld)
    try:
        x = getattr(output, x_axis)
    except AttributeError:
        print(' * Warning: attribute {} not found. using "time" as default'.format(x_axis))
        x = output.time  # x coordinate

    for i, sortlist in enumerate(sortlists):
        # filter
        allcoeffs = []
        for species in filter_by_species:
            if not pp.isHeterogeneous:
                tot_rop_dct = pp.RateOfProductionAnalysis(
                    species, ropa_type="global", number_of_reactions=100, mass_ropa=mass_ropa )
            else:
                tot_rop_dct = pp.RateOfProductionAnalysis_Surface(
                    species, ropa_type="global", number_of_reactions=100, heterogeneous_reactions=heterogeneous_reactions )
            allcoeffs.extend(tot_rop_dct["reaction_indices"])
        if len(allcoeffs) > 0:
            allcoeffs = list(set(allcoeffs))
            rxns_sorted_df = rxns_sorted_df.loc[allcoeffs]

        rxn_class_df = filter_class0(rxns_sorted_df, filter_dcts[i])
        rxn_class_series = sortby0(rxn_class_df, sortlist, dropunsorted=False)

        # rates
        rates_df = pp.reactionrategroups(rxn_class_series, x, threshold=threshs[i])
        sortdfs.append(rates_df)

    return sortdfs