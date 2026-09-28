import os

import matplotlib.pyplot as plt

from pySMOKEPostProcessor import PostProcessor, plot_distribution

kineticFolder = os.path.join("..", "data", "Soot-01", "kinetics")
resultsFolder = os.path.join("..", "data", "Soot-01", "Output")

pp = PostProcessor(kineticFolder, resultsFolder)

# 1) BinProperties can be accessed directly as a pandas.DataFrame
props = pp.dfSootProperties
print(props.head(6).to_string())
print(f"... {len(props)} BINs, sections {props['Bin_section'].min()}"
      f"-{props['Bin_section'].max()}\n")

# ---------------------------------------------------------------------------
# 2) Particle size distribution at a chosen location
#      particle_type = "all"     -> every BIN with Bin_section >= min_section
#                      "primary" -> the primary-particle size distribution
#                                    (PPSD, see SootPPSD.py)
#      diameter_type = "dmob" -> mobility, dm = f(Dpp, numPP); correlation_name
#                                 picks "Kelesidis"/"Sorensen"/"Rissler"
#                      "dpp"  -> primary-particle diameter
#                      "dcol" -> collision diameter
#                      "dva"  -> volume-equivalent sphere diameter
#    local_value picks the profile point like local ROPA
# ---------------------------------------------------------------------------
local_value = 0.35  # cm - near peak soot

# Same particle population ("all"). Only "dmob" applies
# the aggregation transform Dpp * numPP**exp; "dpp"/"dcol"/"dva" read the BIN
psd = {}
for diam in ("dmob", "dcol", "dva"):
    psd[diam] = pp.SootPSD(local_value=local_value, particle_type="all",
                           diameter_type=diam)
    df = psd[diam]
    size_col = next(c for c in df.columns if c.endswith("[nm]") and "_" not in c)
    print(f"all / {diam:4s} @ abscissa = {df.attrs['abscissa']:.4g}  "
          f"(T = {df.attrs['T[K]']:.0f} K): {len(df):2d} sections, "
          f"{size_col} {df[size_col].min():.3g}-{df[size_col].max():.4g} nm, "
          f"total N = {df['N[#/m3]'].sum():.3e} #/m3")

# ---------------------------------------------------------------------------
# 3) Comparison between different diameter calculation method
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.2, 5.0))
for diam, style in zip(("dmob", "dcol", "dva"), ("-o", "-s", "-^", "-d")):
    plot_distribution(psd[diam], ax=ax, label=f"{diam}", linestyle=style[0],
                      marker=style[1])
ax.set_xlabel("Particle diameter [nm]")
ax.set_title("Soot PSD")
ax.set_ylabel(r"$dN/d\log_{10}(d_m)$ [m$^{-3}$]")
ax.legend()
fig.tight_layout()
# fig.savefig("SootPSD.png", dpi=200, bbox_inches="tight")

plt.show()
