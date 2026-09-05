'''
SCRIPT: rxnclasses
@Authors:
    Luna Pratali Maffei [1]
    [1]: CRECK Modeling Lab, Department of Chemistry, Materials, and Chemical Engineering, Politecnico di Milano
@Contacts:
    luna.pratali@polimi.it
@Additional notes:
    This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
    Please report any bug to: luna.pratali@polimi.it
'''
import sys
import pandas as pd
import numpy as np

# TODO Think to move the following two lists in a separate place
# maybe also the function check bimol type
R = ['R', 'H', 'OH', 'O2', 'O', 'CH3', 'HO2', 'HCO', 'C2H3']
RSR = ['C5H5', 'C7H7', 'C6H5O']


def check_bimol_type(speciestype, reactiontype):
    """
    check if a reaction is of a type M+M, RSR+M, R+M, R+R, RSR+RSR
    e.g. TODO
    """
    try:
        species_type_1 = speciestype.split('-')[-1]
    except AttributeError:
        return 'UNSORTED'

    if reactiontype == 'ENLARGE':
        return 'UNIMOL'

    # CHECK SPECIES TYPE 2
    if reactiontype.split('_')[0] in ['HABS', 'ADD', 'REC', 'ENLARGE'] and '-' not in reactiontype:
        if len(reactiontype.split('_')) == 1:
            return 'UNSORTED'
        if reactiontype.split('_')[1] in R:
            species_type_2 = 'R'
        elif reactiontype.split('_')[1] in RSR:
            species_type_2 = 'RSR'
        else:
            species_type_2 = 'NA'
    elif reactiontype.split('_')[0] in ['ADD', 'REC', 'ENLARGE'] and '-' in reactiontype:
        species_type_2 = reactiontype.split('-')[-1]
    else:
        return 'UNIMOL'

    type_list = [species_type_1, species_type_2]
    type_list.sort()

    return '+'.join(type_list)


def filter_class0(rxn_class_df, filter_dct: dict = {}):
    """ only keep classes according to criteria listed in filter_dct """
    if len(filter_dct) > 0:
        dict_indexes = dict.fromkeys(filter_dct.keys())
        for criterion, values in filter_dct.items():
            dict_indexes[criterion] = np.array(
                list(
                    rxn_class_df[
                        [any(v in val for v in values)
                            for val in rxn_class_df[criterion]]
                    ].index
                )
            )

        indexes_filter = np.array(list(set.intersection(
            *[set(val) for val in dict_indexes.values()])))
        rxn_class_df = rxn_class_df.loc[np.array(
            list(set(indexes_filter)))]

    return rxn_class_df


def sortby0(rxn_class_df,
            sortlist: list,
            dropunsorted: bool = True):
    """ group rxn dataframe by listed criteria and  """
    # check that all criteria are columns
    if not all(criterion in rxn_class_df.columns for criterion in sortlist):
        print(' * Error: criteria not all present in dataframe columns - exiting')
        sys.exit()

    # group
    labels = []
    listofnames = []
    if len(sortlist) == 1:
        sortlist = sortlist[0]  # avoid pandas warning for single grouper
    for grp_idx, grp_df in rxn_class_df.groupby(sortlist):
        if isinstance(grp_idx, str):
            label = grp_idx
        else:
            label = '['+']['.join(grp_idx)+']'
        labels.append(label)
        listofnames.append(grp_df['name'].values)

    new_sort_df = pd.Series(listofnames, labels)
    # drop unsorted indexes
    if dropunsorted:
        new_sort_df = new_sort_df[~new_sort_df.index.str.contains('UNSORTED')]

    return new_sort_df


class reaction_classes_assign:

    def __init__(self, reactions, verbose: bool):
        """ allocate dataframe """
        # turn reactions into a dataframe
        # index = rxn index, columns = [classtype, class, reactiontype, flux]
        self.rxn_class_df = pd.DataFrame(index=np.arange(1, len(reactions)+1),
                                         columns=np.array([
                                             'name',
                                             'classtype',
                                             'speciestype',
                                             'reactiontype',
                                             'bimoltype'
                                         ]),
                                         dtype=object)

        self.reactions = reactions
        self.verbose = verbose

        # --- assign class and sublcass
        for rxn in self.reactions:
            idx = rxn['index']
            # self.rxn_class_df['name'][idx] = rxn['name']
            self.rxn_class_df.loc[idx, 'name'] = rxn['name']
            self.rxn_class_df.loc[idx, 'speciestype'] = rxn['class']
            self.rxn_class_df.loc[idx, 'reactiontype'] = rxn['reactiontype']
            self.rxn_class_df.loc[idx, 'bimoltype'] = check_bimol_type(rxn['class'], rxn['reactiontype'])

        self.rxn_class_df = self.rxn_class_df.replace('None', np.nan)
        self.rxn_class_df['speciestype'] = self.rxn_class_df['speciestype'].replace(np.nan, 'UNSORTED')
        self.rxn_class_df['reactiontype'] = self.rxn_class_df['reactiontype'].replace(np.nan, 'UNSORTED')

    def assign_class_grp(self, subcl_grp_dct):
        """ assign class group if available """
        for subcl, subset in self.rxn_class_df.groupby('reactiontype'):
            # None values are automatically discarded
            rxns = subset.index
            if subcl in subcl_grp_dct.keys():
                self.rxn_class_df.loc[rxns, 'classtype'] = subcl_grp_dct[subcl]
            else:
                self.rxn_class_df.loc[rxns, 'classtype'] = 'UNSORTED'
                print(' * Warning: reactiontype {} not found in class groups'.format(subcl))


class reaction_fluxes:
    def __init__(self, rxn_class_df, verbose: bool):
        """ initialize rxn class dataframe """
        self.rxn_class_df = rxn_class_df
        self.verbose = verbose

    def assign_flux(self, tot_rop_df):
        """ assign flux from the flux analysis """
        # print(self.rxn_class_df.loc[tot_rop_df.index])
        self.rxn_class_df = pd.concat([self.rxn_class_df, tot_rop_df], axis=1)
        # self.rxn_class_df.join(tot_rop_df)
        # replace missing flux values with 0
        flux_sp_name = tot_rop_df.columns[0]
        self.rxn_class_df[flux_sp_name] = self.rxn_class_df[flux_sp_name].replace(np.nan, 0.0)

    # renormalize
    # renorm_factor = sum(abs(self.rxn_class_df[flux_sp_name]))

    # renormalize
    # rxn_class_df = rxn_class_df.sort_values(by='absflux', ascending = False)

    def finalizeFlux(self):
        """
        Prepares rxn_class_df after assign_flux, before sort_and_filter takes over.
        Updated version of the previous netfluxes function.

        Merging duplicate/reversible-reaction pairs now happens upstream in C++ 
        (ReactionClass.mergeDuplicates, called from FluxByClass.process_flux)
        using the KineticMap stoichiometric matrix).
        By the time rxn_class_df reaches here every duplicate group is
        already collapsed to one row. What's left, all bookkeeping filter_flux/sortby
        depend on afterwards:
          - self.flux_cols: which columns are per-species flux columns.
          - drop rows with no contribution assigned on any of them.
          - self.rxn_class_df_all: a frozen copy of the table at this point, before
            sort_and_filter's filter_class/filter_flux prune rxn_class_df further -
            sortby(weigh='normbyspecies') needs it as the "before filtering" baseline
            to compute what fraction of total flux survived filtering.
        """
        self.flux_cols = np.array([col for col in self.rxn_class_df.columns if 'flux' in col], dtype=str)
        self.rxn_class_df = self.rxn_class_df[~(self.rxn_class_df[self.flux_cols] == 0.).all(axis=1)]
        self.rxn_class_df_all = self.rxn_class_df.copy()

    def filter_class(self, filter_dct):
        """ only keep classes according to criteria listed in filter_dct """
        self.rxn_class_df = filter_class0(self.rxn_class_df, filter_dct)

    def filter_flux(self, threshold=1e-3):
        """ delete all reactions with contributions below a threshold """

        indexes_filter = []

        for flux_sp_name in self.flux_cols:
            col = self.rxn_class_df[flux_sp_name]

            if col.empty:   # [LG] Addded for time-integral computations (soot)
                continue

            max_val = col.abs().max()
            if max_val == 0:
                continue

            mask = col.abs() / max_val > threshold
            indexes_filter.extend(self.rxn_class_df[mask].index.tolist())

        if len(indexes_filter) == 0:
            self.rxn_class_df = self.rxn_class_df.iloc[0:0]
        else:
            self.rxn_class_df = self.rxn_class_df.loc[list(set(indexes_filter))]

    def sortby(self,
               sortlist,
               weigh: str = 'false',
               dropunsorted: bool = True):
        """ sum fluxes by criteria in sortlist """

        # check that all criteria are columns
        if not all(criterion in self.rxn_class_df.columns for criterion in sortlist):
            print(' * Error: criteria not all present in dataframe columns - exiting')
            sys.exit()

        # group and sum
        new_sort_df = pd.DataFrame(index=self.flux_cols)
        if len(sortlist) == 1:
            sortlist = sortlist[0]  # avoid pandas warning for single grouper
        for grp_idx, grp_df in self.rxn_class_df.groupby(sortlist):
            if isinstance(grp_idx, str):
                name = grp_idx
            else:
                name = '['+']['.join(grp_idx)+']'
            if self.verbose:
                print(grp_idx, '\n', grp_df, '\n')
            new_sort_df[name] = grp_df[self.flux_cols].sum()

        # drop unsorted columns
        col_names = new_sort_df.columns
        for col in col_names:
            if 'UNSORTED' in col and dropunsorted:
                new_sort_df = new_sort_df.drop(col, axis=1)
                print(' * Warning: dropping UNSORTED rxns prior to normalization')

        # renormalize by species
        if weigh == 'normbyspecies':
            for flux_sp_name in new_sort_df.index:
                renorm_factor = max(abs(new_sort_df.loc[flux_sp_name]))
                weight_factor = sum(abs(self.rxn_class_df[flux_sp_name]))/sum(abs(self.rxn_class_df_all[flux_sp_name]))
                # new_sort_df.loc[flux_sp_name] /= renorm_factor
                new_sort_df.loc[flux_sp_name] *= (weight_factor/renorm_factor)
        elif weigh == 'omegaij':
            new_sort_df = self.omegaij(new_sort_df)
        # else: # unneeded - now plot function works with any kind of axis
            # renorm_factor = abs(self.rxn_class_df_all[new_sort_df.index]).max().max()
            # new_sort_df /= renorm_factor
        elif weigh != 'false':
            print(' * Error: "weigh" can be (str) normbyspecies, omegaij, false')
            # if weigh = false, you are doing nothing

        if self.verbose:
            # questo print e probabilmente anche alcuni di quelli sopra forse hanno poco
            # senso e magari potremmo dare la possibilità di stampare su file?
            print(new_sort_df)

        return new_sort_df

    def omegaij(self, sort_df):
        """
        takes a dataframe and computes the omegaijn coefficients
        see https://doi.org/10.1016/j.combustflame.2022.112073
        """
        # new df with same indexes and columns as the original
        df_omegaijn = pd.DataFrame(index=sort_df.index, columns=sort_df.columns, dtype=np.float64)

        # omega min and omega max
        omegamin = np.min(sort_df.values)
        omegamax = np.max(sort_df.values)

        # if min and max are not as expected: print warning:
        if omegamin > 0:
            print(' * Warning: value of omega min is > 0')
        if omegamax < 0:
            print(' * Warning: value of omega max is < 0')

        # compute
        for sp in sort_df.index:
            for cl in sort_df.columns:
                val = sort_df[cl][sp]
                if val > 0.:
                    df_omegaijn.loc[sp, cl] = val/omegamax
                elif val < 0.:
                    df_omegaijn.loc[sp, cl] = val/omegamin
                else:
                    df_omegaijn.loc[sp, cl] = 0.

        return df_omegaijn
