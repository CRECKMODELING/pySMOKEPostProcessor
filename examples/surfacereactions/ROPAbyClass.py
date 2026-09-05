import os
from pySMOKEPostProcessor.postprocessor import PostProcessor
import matplotlib.pyplot as plt
from pySMOKEPostProcessor.plotting_utilities.bar_plot import plot_bars
from pySMOKEPostProcessor import plot_heatmap
from pySMOKEPostProcessor import script_utils

kineticFolder = os.path.join("..","data","Surface_Data","kinetics")
kineticFile = os.path.join(kineticFolder,"kinetics.surface.xml")
# resultsFolder = os.path.join("..","data","Surface_Data","Output")
resultsFolder ={'Output' : os.path.join('..', 'data', 'Surface_Data', 'Output')}

class_groups_file = os.path.join('..', 'data', 'Surface_Data', 'surf_rxn_class.txt')

# pp = PostProcessor(kineticFolder, resultsFolder)
# sort the reactions - same for all examples
# (get_sortedrxns needs a PostProcessor: classification is read off the same
# ProfilesDatabase every other analysis on it uses. Classification only depends on
# the mechanism, so any one of resultsFolder's simulations does - each is still
# built fresh inside the process_classes loop below for its own ROPA.)
pp = PostProcessor(kineticFolder, next(iter(resultsFolder.values())))
rxns_sorted = script_utils.get_sortedrxns(pp, class_groups_file, heterogeneous_reactions=True)
#   print(rxns_sorted.rxn_class_df)   # This works properly
species_list =['C(B)']
sortlists = [['classtype']]
ropa_type = 'global'

for simul_name, simul_fld in resultsFolder.items():
    sortdfs = script_utils.process_classes(simul_fld, kineticFolder, rxns_sorted, species_list,
                                           sortlists, ropa_type, n_of_rxns = 20, heterogeneous_reactions=True)
print(sortdfs)
# plot heatmap
fig = plot_heatmap(sortdfs[0], symmetricaxis = True, dftype = 'flux')
fig.set_size_inches(6,4)
fig.tight_layout()
plt.show()
