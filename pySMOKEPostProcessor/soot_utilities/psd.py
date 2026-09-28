import pandas as pd

# plotting lives with the other plot helpers; re-exported here for convenience.
from ..plotting_utilities.plot_distributions import plot_distribution

_DIAMETER_SYMBOL = {"dmob": "dm", "dpp": "dpp", "dcol": "dcol", "dva": "dva"}

# Mobility-diameter correlations from literature.
#   "Kelesidis" -> dpp * numPP**0.45 
#       (Kelesidis et al. 2017, Carbon, 10.1016/j.carbon.2017.06.004).
#   "Sorensen"  -> dpp * numPP**0.465 
#       (Sorensen 2011, Aerosol Sci. Technol., 10.1080/02786826.2011.560909).
#   "Rissler"   -> 0.794 * dpp * numPP**0.51 
#       (Rissler et al. 2013, Aerosol Sci. Technol., 10.1080/02786826.2013.791381).
_MOBILITY_CORRELATIONS = {"Kelesidis", "Sorensen", "Rissler"}

_COLUMN_ORDER = [
    "{d}[nm]", "{d}_min[nm]", "{d}_max[nm]", "N[#/m3]", "n_bins",
    "Dlog10({d}[nm])", "dN/dlog10({d}[nm])[#/m3]",
]

_RAW_TO_PRETTY = {
    "x": "{d}[nm]",
    "x_min": "{d}_min[nm]",
    "x_max": "{d}_max[nm]",
    "N_per_m3": "N[#/m3]",
    "n_bins": "n_bins",
    "Dlog10_x": "Dlog10({d}[nm])",
    "dN_dlog10_x_per_m3": "dN/dlog10({d}[nm])[#/m3]",
}


def compute_psd(pp, local_value, *, particle_type="all", diameter_type="dmob",
                min_section=1, correlation_name="Kelesidis", merge_tol=0.20):
    """Soot particle size distribution at ``local_value``. Returns a DataFrame.    
    particle_type == "primary" toggles the PPSD calculation rather than PSD,
    requiring diameter_type = "dpp"
    """
    if pp.dfSootProperties is None:
        raise ValueError("The mechanism has no <SootProperties> block")
    if diameter_type not in _DIAMETER_SYMBOL:
        raise ValueError('diameter_type must be "dmob", "dpp", "dcol" or "dva"')
    if correlation_name not in _MOBILITY_CORRELATIONS:
        raise ValueError(
            f"correlation_name must be one of {sorted(_MOBILITY_CORRELATIONS)}"
        )
    if particle_type == "primary" and diameter_type != "dpp":
        raise ValueError(
            'particle_type == "primary" (PPSD) is only meaningful with '
            'diameter_type == "dpp"'
        )

    cols, meta = pp.soot.psd(
        local_value, particle_type, diameter_type, min_section, correlation_name,
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
            "correlation_name": correlation_name,
            "local_value": float(local_value),
            "abscissa": meta["abscissa"],
            "T[K]": meta["T"],
            "P[Pa]": meta["P"],
            "rho[kg/m3]": meta["rho"],
            "MW[kg/kmol]": meta["MW"],
        }
    )
    return df
