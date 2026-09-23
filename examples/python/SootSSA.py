import os

import matplotlib.pyplot as plt

from pySMOKEPostProcessor import PostProcessor

kineticFolder = os.path.join("..", "data", "Soot-01", "kinetics")
resultsFolder = os.path.join("..", "data", "Soot-01", "Output")

pp = PostProcessor(kineticFolder, resultsFolder)

# ---------------------------------------------------------------------------
# Mean SSA along the abscissa. Mass-weighted over the BINs with
# Bin_section >= min_section: total soot surface / total soot mass.
# as_dataframe=True adds the soot mass concentration.
# ---------------------------------------------------------------------------
x, ssa = pp.getSSASootProfile(min_section=5)

fig, ax1 = plt.subplots(figsize=(7.2, 5.0))
ax1.plot(x, ssa, "-")
ax1.set_xlabel("axial coordinate [cm]")
ax1.set_ylabel(r"$SSA_{avg}$ [$m^2$/$kg$]")
ax1.set_yscale("log")
ax1.set_title("Soot Surface Specific Area")
fig.tight_layout()

# ---------------------------------------------------------------------------
# Mean H/C along the abscissa over the same BINs: mass-weighted, and
# carbon-weighted (= total H atoms / total C atoms). Any per-BIN property can
# be averaged the same way: pp.soot.averagedProfile(prop, 5, "mass")
# ---------------------------------------------------------------------------
x, htoc_mass = pp.getHtoCSootProfile(min_section=5, weighting="mass")
x, htoc_carbon = pp.getHtoCSootProfile(min_section=5, weighting="carbon")

fig, ax2 = plt.subplots(figsize=(7.2, 5.0))
ax2.plot(x, htoc_mass, "-", label="mass-weighted")
ax2.plot(x, htoc_carbon, "--", label="carbon-weighted")
ax2.legend()
ax2.set_xlabel("axial coordinate [cm]")
ax2.set_ylabel(r"$H/C_{avg}$ [-]")
ax2.set_title("Soot H/C ratio")
fig.tight_layout()
plt.show()
