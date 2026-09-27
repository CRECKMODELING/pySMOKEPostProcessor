import os
# Run from examples/python folder
kineticFolder = os.path.join("..","data","Surface_Data","kinetics")
resultsFolder = os.path.join("..","data","Surface_Data","Output")

from pySMOKEPostProcessor.postprocessor import PostProcessor
import matplotlib.pyplot as plt
# from postprocessor import PostProcessor
from pySMOKEPostProcessor.plotting_utilities.bar_plot import plot_bars

pp = PostProcessor(kineticFolder, resultsFolder)

sp = 'C(B)'
phase_ropa = True

ropa = pp.RateOfProductionAnalysis(
                    species=sp,
                    ropa_type='global',
                    number_of_reactions=10,
                    heterogeneous_reactions=phase_ropa)

result = {
    "coefficients": ropa["coefficients"],
    "reaction_names": ropa["reaction_names"],
    "reaction_indices": ropa["reaction_indices"]
}

fig1, ax1 = plot_bars(result)
fig1.suptitle(f"Species: {sp}", fontweight="bold")
if phase_ropa:
    ax1.set_title("Heterogeneous")
else:
    ax1.set_title("Homogeneous")
fig1.set_size_inches(10, 6)
fig1.tight_layout()

# global_ropa_het = pp.RateOfProductionAnalysis(species='CH4',
#                                           ropa_type='global',
#                                           number_of_reactions=10,
#                                           heterogeneous_reactions=True)

# result_het = {
#     "coefficients": global_ropa_het["coefficients"],
#     "reaction_names": global_ropa_het["reaction_names"],
#     "reaction_indices": global_ropa_het["reaction_indices"]
# }

# fig_2, ax_2 = plot_bars(result_het)
# fig_2.suptitle("Heterogeneous ROPA")
# fig_2.tight_layout()

# local_ropa = pp.RateOfProductionAnalysis(species='H2',
#                                          ropa_type='local',
#                                          local_value=0.001,
#                                          number_of_reactions=10)

# fig_2, ax_2 = plot_bars(local_ropa)

# region_ropa = pp.RateOfProductionAnalysis(species='H2',
#                                           ropa_type='region',
#                                           lower_value=0.0005,
#                                           upper_value=0.0009,
#                                           number_of_reactions=15)

# fig_3, ax_3 = plot_bars(region_ropa)

plt.show()
