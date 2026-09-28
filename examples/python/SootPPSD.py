import os

import matplotlib.pyplot as plt

from pySMOKEPostProcessor import PostProcessor, plot_distribution

kineticFolder = os.path.join("..", "data", "Soot-01", "kinetics")
resultsFolder = os.path.join("..", "data", "Soot-01", "Output")

pp = PostProcessor(kineticFolder, resultsFolder)

local_value = 0.35  # cm - near peak soot

ppsd = pp.SootPSD(local_value=local_value, particle_type="primary", diameter_type="dpp")
print(f"PPSD @ abscissa = {ppsd.attrs['abscissa']:.4g} "
      f"(T = {ppsd.attrs['T[K]']:.0f} K): {len(ppsd):2d} dpp classes, "
      f"dpp[nm] {ppsd['dpp[nm]'].min():.3g}-{ppsd['dpp[nm]'].max():.4g} nm, "
      f"total primary-particle N = {ppsd['N[#/m3]'].sum():.3e} #/m3")
print(ppsd[["dpp[nm]", "n_bins", "N[#/m3]"]].to_string(index=False))


fig, ax = plt.subplots(figsize=(8.0, 5.0))
plot_distribution(ppsd, ax=ax, label="PPSD")
ax.set_title("Primary Particle Size Distribution")

fig.tight_layout()
# fig.savefig("SootPPSD.png", dpi=200, bbox_inches="tight")

plt.show()
