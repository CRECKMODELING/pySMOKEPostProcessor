import os

import matplotlib.pyplot as plt

from pySMOKEPostProcessor import PostProcessor, plot_distribution

kineticFolder = os.path.join("..", "data", "Soot-01", "kinetics")
resultsFolder = os.path.join("..", "data", "Soot-01", "Output")

pp = PostProcessor(kineticFolder, resultsFolder)

# ---------------------------------------------------------------------------
# 1) The full <SootProperties> block, one row per BIN. Column names match
#    OpenSMOKE's BinProperties.txt; Bin_index / Bin_name tie each BIN back to
#    its gas-phase species.
# ---------------------------------------------------------------------------
props = pp.dfSootProperties
print(props.head(6).to_string())
print(f"... {len(props)} BINs, sections {props['Bin_section'].min()}"
      f"-{props['Bin_section'].max()}\n")

# ---------------------------------------------------------------------------
# 2) Particle size distribution at a chosen abscissa location. One function,
#    two knobs:
#      particle_type = "all"     -> every BIN with Bin_section >= min_section
#                      "primary" -> free primary particles (numPP == 1)  [PPSD]
#      diameter_type = "dmob" -> mobility, dm = Dpp * numPP**mobility_exponent
#                      "dpp"  -> primary-particle diameter (straight from XML)
#                      "dcol" -> collision diameter (straight from XML)
#    local_value picks the profile point like local ROPA (first point whose
#    abscissa >= local_value); the gas state used is in df.attrs.
# ---------------------------------------------------------------------------
local_value = 0.35  # s - near peak soot for this batch-reactor fixture

# Same particle population ("all"), three size coordinates. Only "dmob" applies
# the aggregation transform Dpp * numPP**exp; "dpp"/"dcol" read the BIN property
# directly, so total N is identical - only the abscissa (and its binning) moves.
psd = {}
for diam in ("dmob", "dpp", "dcol"):
    psd[diam] = pp.SootPSD(local_value=local_value, particle_type="all",
                           diameter_type=diam)
    df = psd[diam]
    size_col = next(c for c in df.columns if c.endswith("[nm]") and "_" not in c)
    print(f"all / {diam:4s} @ abscissa = {df.attrs['abscissa']:.4g}  "
          f"(T = {df.attrs['T[K]']:.0f} K): {len(df):2d} sections, "
          f"{size_col} {df[size_col].min():.3g}-{df[size_col].max():.4g} nm, "
          f"total N = {df['N[#/m3]'].sum():.3e} #/m3")

print()
print("mobility PSD (all / dmob):")
print(psd["dmob"].to_string(), "\n")

# Free primary particles only (numPP == 1) - the classic PPSD, on Dpp.
ppsd = pp.SootPSD(local_value=local_value, particle_type="primary",
                  diameter_type="dpp")
print(f"PPSD (primary / dpp) @ abscissa = {ppsd.attrs['abscissa']:.4g}")
print(ppsd.to_string(), "\n")

# ---------------------------------------------------------------------------
# 3) Plot the three "all" size coordinates on one axes.
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.2, 5.0))
for diam, style in zip(("dmob", "dpp", "dcol"), ("-o", "-s", "-^")):
    plot_distribution(psd[diam], ax=ax, label=f"all / {diam}", linestyle=style[0],
                      marker=style[1])
ax.set_xlabel("particle diameter [nm]")
ax.set_title("Soot PSD - mobility vs primary vs collision diameter")
ax.legend()
fig.tight_layout()
# fig.savefig("SootPSD.png", dpi=200, bbox_inches="tight")

plt.show()
