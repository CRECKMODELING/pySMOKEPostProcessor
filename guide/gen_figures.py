"""
Generates the PNG figures embedded in guide/UserGuide.tex, straight from the
public API using the fixture data under examples/data/. Not part of the
package; run once with:

    conda run -n pdev python guide/gen_figures.py

from the repository root. Every figure is saved at 400 dpi into guide/figures/
(this guide is a digital reference, not meant to be printed, so figures are
sized generously for on-screen zooming rather than for page real estate).
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pySMOKEPostProcessor import (
    PostProcessor,
    plot_bars,
    plot_heatmap,
    plot_class_distribution,
    plot_distribution,
    plot_areas,
    elements_balance,
    script_utils,
)
from pySMOKEPostProcessor.surfacereactions_utilities.deposition_plot import (
    CumulativeDeposition,
    CumulativeSootProduction,
)

DPI = 400

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "examples", "data")
OUT = os.path.join(ROOT, "guide", "figures")
os.makedirs(OUT, exist_ok=True)


def savefig(fig, name, **kw):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", **kw)
    plt.close(fig)
    print("wrote", path)


# ---------------------------------------------------------------------------
# 1) ROPA bar plot (global)
# ---------------------------------------------------------------------------
pp = PostProcessor(os.path.join(DATA, "ROPA", "kinetics"), os.path.join(DATA, "ROPA", "Output"))
ropa = pp.RateOfProductionAnalysis(species="H2", ropa_type="global", number_of_reactions=10)
fig, ax = plot_bars(ropa)
fig.set_size_inches(7.5, 4.5)
ax.set_title("Global ROPA on H2")
fig.tight_layout()
savefig(fig, "ropa_bars.png")

# ---------------------------------------------------------------------------
# 1b) ROPA bar plots: local vs region, side by side
# ---------------------------------------------------------------------------
local_ropa = pp.RateOfProductionAnalysis(species="H2", ropa_type="local", local_value=0.001,
                                          number_of_reactions=10)
region_ropa = pp.RateOfProductionAnalysis(species="H2", ropa_type="region", lower_value=0.0005,
                                           upper_value=0.0009, number_of_reactions=10)
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
_, axes[0] = plot_bars(local_ropa, ax=axes[0])
axes[0].set_title("local (t = 0.001 s)")
_, axes[1] = plot_bars(region_ropa, ax=axes[1])
axes[1].set_title("region (t = 0.0005-0.0009 s)")
fig.suptitle("ROPA on H2 - local vs. region", fontweight="bold")
fig.tight_layout()
savefig(fig, "ropa_local_region.png")

# ---------------------------------------------------------------------------
# 2) Sensitivity bar plot
# ---------------------------------------------------------------------------
pp_s = PostProcessor(os.path.join(DATA, "Sensitivity", "kinetics"), os.path.join(DATA, "Sensitivity", "Output"))
sens = pp_s.SensitivityAnalysis(target="NO", sensitivity_type="global", number_of_reactions=15,
                                 ordering_type="peak-values", normalization_type="local")
fig, ax = plot_bars(sens)
fig.set_size_inches(7.5, 4.5)
ax.set_title("Global sensitivity of NO")
fig.tight_layout()
savefig(fig, "sensitivity_bars.png")

# ---------------------------------------------------------------------------
# 3) Flux analysis graphs (graphviz -> png): destruction (red) and production (blue)
# ---------------------------------------------------------------------------
G_destr = pp.FluxAnalysis(species="H2", element="H", flux_analysis_type="destruction",
                           thickness="relative", thickness_log_scale=True, label_type="relative",
                           depth=1, width=3, threshold=0.01, local_value=0.001)
G_destr.format = "png"
G_destr.attr(dpi=str(DPI))
G_destr.render(filename="flux_graph_destr", directory=OUT, cleanup=True)
os.replace(os.path.join(OUT, "flux_graph_destr.png"), os.path.join(OUT, "flux_analysis.png"))
print("wrote", os.path.join(OUT, "flux_analysis.png"))

G_prod = pp.FluxAnalysis(species="H2", element="H", flux_analysis_type="production",
                          thickness="relative", thickness_log_scale=True, label_type="relative",
                          depth=1, width=3, threshold=0.01, local_value=0.001)
G_prod.format = "png"
G_prod.attr(dpi=str(DPI))
G_prod.render(filename="flux_graph_prod", directory=OUT, cleanup=True)
os.replace(os.path.join(OUT, "flux_graph_prod.png"), os.path.join(OUT, "flux_analysis_production.png"))
print("wrote", os.path.join(OUT, "flux_analysis_production.png"))

# ---------------------------------------------------------------------------
# 4) Species-class elemental distribution (stacked area)
# ---------------------------------------------------------------------------
pp_c = PostProcessor(os.path.join(DATA, "SpeciesClasses", "kinetics"), os.path.join(DATA, "SpeciesClasses", "Output"))
carbon = pp_c.ElementalDistributionByClass(element="C", normalize=True)
fig, ax = plot_class_distribution(carbon, xlabel="time [s]", ylabel="carbon fraction",
                                   title="Carbon distribution by species class")
fig.set_size_inches(7.2, 4.5)
fig.tight_layout()
savefig(fig, "species_class_distribution.png")

# ---------------------------------------------------------------------------
# 5) Flux-by-class heatmap
# ---------------------------------------------------------------------------
res = pp_c.FluxAnalysisByClass(element="C", flux_analysis_type="destruction", class_name="C2",
                                flux_per_class=True, depth=3, width=5, threshold=1.0,
                                local_value=0.5, carbon_weighted=True)
fig = plot_heatmap(res["matrix"], symmetricaxis=False, dftype="generic")
savefig(fig, "flux_by_class_heatmap.png")

# ---------------------------------------------------------------------------
# 6) Soot PSD
# ---------------------------------------------------------------------------
pp_soot = PostProcessor(os.path.join(DATA, "Soot-01", "kinetics"), os.path.join(DATA, "Soot-01", "Output"))
fig, ax = plt.subplots(figsize=(7.2, 5.0))
for diam, marker in zip(("dmob", "dpp", "dcol", "dva"), ("o", "s", "^", "d")):
    df = pp_soot.SootPSD(local_value=0.35, particle_type="all", diameter_type=diam)
    plot_distribution(df, ax=ax, label=f"all / {diam}", marker=marker)
ax.set_xlabel("particle diameter [nm]")
ax.set_title("Soot PSD - mobility vs primary vs collision vs volume-equivalent diameter")
ax.legend()
fig.tight_layout()
savefig(fig, "soot_psd.png")

# ---------------------------------------------------------------------------
# 6b) Soot volume fraction fv along the abscissa
# ---------------------------------------------------------------------------
x, fv = pp_soot.getFvSootProfile(min_section=5)
fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax.plot(x, fv, "-")
ax.set_xlabel("axial coordinate [cm]")
ax.set_ylabel(r"$f_v$ [-]")
ax.set_yscale("log")
ax.set_title("Soot volume fraction")
fig.tight_layout()
savefig(fig, "soot_fv.png")

# ---------------------------------------------------------------------------
# 7) Reaction-rate-by-class heatmap (surface deposition example)
# ---------------------------------------------------------------------------
kin_surf = os.path.join(DATA, "Surface_Data", "kinetics")
out_surf = os.path.join(DATA, "Surface_Data", "Output")
class_groups_file = os.path.join(DATA, "Surface_Data", "surf_rxn_class.txt")
pp_surf = PostProcessor(kin_surf, out_surf)
rxns_sorted = script_utils.get_sortedrxns(pp_surf, class_groups_file, heterogeneous_reactions=True)
sortdfs = script_utils.process_classes(out_surf, kin_surf, rxns_sorted, ["C(B)"], [["classtype"]],
                                        "global", n_of_rxns=20, heterogeneous_reactions=True, pp=pp_surf)
fig = plot_heatmap(sortdfs[0], symmetricaxis=True, dftype="flux")
fig.set_size_inches(6, 4)
fig.tight_layout()
savefig(fig, "ropa_by_class_heatmap.png")

# ---------------------------------------------------------------------------
# 8) Reaction-rate line plot (GetReactionRates)
# ---------------------------------------------------------------------------
rate = pp.GetReactionRates(reaction_name=["O2+H=O+OH"])[0]
_, temperature = pp.getTemperatureProfile()

names = ["O2+H=O+OH", "H2+O=H+OH", "O2+H(+M)=HO2(+M)"]
rates = pp.GetReactionRates(reaction_name=names)

fig, ax = plt.subplots(figsize=(7.2, 4.8))
for r, c, n in zip(rates, ("r", "b", "g"), names):
    ax.plot(temperature, r, c, label=n)
ax.set_xlabel("Temperature [K]", fontsize=13)
ax.set_ylabel(r"Reaction Rate $\left[ \dfrac{kmol}{m^{3}s} \right]$", fontsize=13)
ax.grid()
ax.legend()
fig.tight_layout()
savefig(fig, "reaction_rates_lines.png")

# ---------------------------------------------------------------------------
# 9) Cumulative rates (shaded-area plot)
# ---------------------------------------------------------------------------
cum = script_utils.cumulative_rates(os.path.join(DATA, "ROPA", "Output"),
                                     os.path.join(DATA, "ROPA", "kinetics"),
                                     species_list=["H2"], rate_type="PC",
                                     x_axis=None, n_of_rxns=100, threshold=0.02, pp=pp)
fig, ax = plot_areas(cum["H2"], xlabel="time [s]", ylabel=r"rate [kmol/m$^3$s]",
                      title="Cumulative net reaction rate for H2")
fig.set_size_inches(7.5, 4.8)
fig.tight_layout()
savefig(fig, "cumulative_rates_area.png")

# ---------------------------------------------------------------------------
# 10) Elements balance (elements_balance() calls plt.show() internally; on the
#     Agg backend plt.show() is a no-op that leaves the figure open, so it can
#     still be grabbed via plt.gcf() right after the call).
# ---------------------------------------------------------------------------
elements_balance(os.path.join(DATA, "ROPA", "kinetics"), os.path.join(DATA, "ROPA", "Output"),
                  elements_list=["C", "H", "O"], threshold=0.05)
fig = plt.gcf()
savefig(fig, "elements_balance.png")

# ---------------------------------------------------------------------------
# 11) Cumulative deposition + soot production (stacked area, two panels)
# ---------------------------------------------------------------------------
surface_classes = os.path.join(DATA, "Surface_Data", "surf_rxn_class.txt")
gas_classes = os.path.join(DATA, "ReactionClasses", "rxn_class_groups.txt")
sortlists = [["reactiontype"]]
lump_steps = 4

Bulk = CumulativeDeposition(kineticFolder=kin_surf, outputFolder=out_surf,
                             class_group_file=surface_classes, sort_type=sortlists,
                             area=200. / 1e4, allCarbon=True)
Soot = CumulativeSootProduction(kineticFolder=kin_surf, outputFolder=out_surf,
                                 class_group_file=gas_classes, sort_type=sortlists,
                                 volume=100. / 1e4)

fig, ax = plt.subplots(2, 1, figsize=(10.5, 12.0), sharex=True)
fig, ax[0] = Bulk.plotCumulativeDeposition(lump_steps=lump_steps, units="mass", fig=fig, ax=ax[0])
fig, ax[1] = Soot.plotSootProduction(lump_steps=lump_steps, units="mass", fig=fig, ax=ax[1])
# Deposited carbon mass and soot mass are not on comparable scales for this
# fixture (area/volume are illustrative constants, not matched to each
# other), so each panel keeps its own natural y-limits rather than being
# forced to share one - only the bottom panel is flipped, for the "meets in
# the middle" mirrored layout.
ax[1].invert_yaxis()
fig.tight_layout()
savefig(fig, "deposition_soot.png")

print("Done.")
