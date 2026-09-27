import os

import matplotlib.pyplot as plt

from pySMOKEPostProcessor import PostProcessor

kineticFolder = os.path.join("..", "data", "Soot-01", "kinetics")
resultsFolder = os.path.join("..", "data", "Soot-01", "Output")

pp = PostProcessor(kineticFolder, resultsFolder)

# ---------------------------------------------------------------------------
# Soot volume fraction fv [-] along the abscissa: getFvSootProfile() sums, at
# every abscissa point, (1/rho_BIN) * mass_BIN.
# This is the same as OpenSMOKE; the definition is VSoot/VTot
# ---------------------------------------------------------------------------
x, fv = pp.getFvSootProfile(min_section=5)

fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax.plot(x, fv, "-")
ax.set_xlabel("axial coordinate [cm]")
ax.set_ylabel(r"$f_v$ [-]")
ax.set_yscale("log")
ax.set_title("Soot volume fraction")
fig.tight_layout()
plt.show()
