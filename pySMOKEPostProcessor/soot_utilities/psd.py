"""DataFrame + plotting front-end for the C++ soot size distribution.

The maths lives in cxx_wrapper/source/Soot.hpp (``Soot::ParticleSizeDistribution``),
which reads the BIN properties *and* the state profiles straight from the same
ProfilesDatabase every other analysis uses - no Output.xml re-read. Here we only
shape the returned column dict into a DataFrame and offer a plot helper.

One distribution, a few knobs:

* ``particle_type`` - ``"all"`` keeps every BIN with ``Bin_section >= min_section``
  (numPP <= 0 nascent bins count as one spherule); ``"primary"`` keeps only the
  free primary particles (numPP == 1), and ignores ``min_section``.
* ``diameter_type`` - ``"dmob"`` mobility diameter dm = Dpp * numPP**exponent,
  ``"dpp"`` primary-particle diameter, ``"dcol"`` collision diameter.

``local_value`` selects the profile point like local ROPA does (first point whose
abscissa >= local_value); the abscissa is time for a reactor, a coordinate for a
flame. The gas state actually used lands in ``df.attrs``.
"""

import pandas as pd

# plotting lives with the other plot helpers; re-exported here for convenience.
from ..plotting_utilities.psd_plot import plot_distribution  # noqa: F401

_DIAMETER_SYMBOL = {"dmob": "dm", "dpp": "dpp", "dcol": "dcol"}

_COLUMN_ORDER = [
    "{d}[nm]", "{d}_min[nm]", "{d}_max[nm]", "N[#/m3]", "n_bins",
    "Dlog10({d}[nm])", "dN/dlog10({d}[nm])[#/m3]",
]

_RAW_TO_PRETTY = {
    "d_nm": "{d}[nm]",
    "d_min_nm": "{d}_min[nm]",
    "d_max_nm": "{d}_max[nm]",
    "N_per_m3": "N[#/m3]",
    "n_bins": "n_bins",
    "Dlog10_d": "Dlog10({d}[nm])",
    "dN_dlog10_d_per_m3": "dN/dlog10({d}[nm])[#/m3]",
}


def compute_psd(pp, local_value, *, particle_type="all", diameter_type="dmob",
                min_section=5, mobility_exponent=0.45, merge_tol=0.20):
    """Soot particle size distribution at ``local_value``. Returns a DataFrame."""
    if pp.dfSootProperties is None:
        raise ValueError("The mechanism has no <SootProperties> block")
    if diameter_type not in _DIAMETER_SYMBOL:
        raise ValueError('diameter_type must be "dmob", "dpp" or "dcol"')

    cols, meta = pp.soot.psd(
        local_value, particle_type, diameter_type, min_section, mobility_exponent,
        merge_tol,
    )

    symbol = _DIAMETER_SYMBOL[diameter_type]
    rename = {k: v.format(d=symbol) for k, v in _RAW_TO_PRETTY.items()}
    df = pd.DataFrame(cols).rename(columns=rename)
    df = df[[c.format(d=symbol) for c in _COLUMN_ORDER]]
    df["n_bins"] = df["n_bins"].astype(int)
    df.attrs.update(
        {
            "particle_type": particle_type,
            "diameter_type": diameter_type,
            "min_section": min_section,
            "local_value": float(local_value),
            "abscissa": meta["abscissa"],
            "T[K]": meta["T"],
            "P[Pa]": meta["P"],
            "rho[kg/m3]": meta["rho"],
            "MW[kg/kmol]": meta["MW"],
        }
    )
    return df
