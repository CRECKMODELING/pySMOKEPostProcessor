import matplotlib.pyplot as plt
import numpy as np

_LABELS = {
    "dm[nm]": ("Mobility diameter $d_m$ [nm]", "d_m"),
    "dpp[nm]": ("Primary-particle diameter $d_{pp}$ [nm]", "d_{pp}"),
    "dcol[nm]": ("Collision diameter $d_{col}$ [nm]", "d_{col}"),
    "dva[nm]": ("Volume-equivalent diameter $d_{va}$ [nm]", "d_{va}"),
}


def plot_distribution(df, ax=None, label=None, **plot_kw):
    """Plot for a SootPSD DataFrame; ``d`` vs either ``dN/dlog10(d)`` or ``N``.

    ``df`` is what ``PostProcessor.SootPSD`` returns: a size column named
    ``dm[nm]`` / ``dpp[nm]`` / ``dcol[nm]`` / ``dva[nm]`` plus the matching
    ``dN/dlog10(<d>[nm])[#/m3]`` and ``N[#/m3]``. Points with a non-finite or
    non-positive coordinate/ordinate are dropped. Returns the Axes.

    In case ``df.attrs["particle_type"] == "primary"`` (PPSD), 
    the plotted distribution is ``N(dpp)``.
    Otherwise, the plotted distribution is ``dN/dlog10(dx)`` 
    where x may be ``dva``, ``dmob``, ``dcol``
    """
    size_col = next(c for c in df.columns if c in _LABELS)
    xlabel, sym = _LABELS[size_col]
    x = df[size_col].to_numpy(dtype=float)

    if ax is None:
        _, ax = plt.subplots(figsize=(7.2, 5.0))
    if label is None:
        label = f"abscissa = {df.attrs.get('abscissa', float('nan')):.3g}"

    if df.attrs.get("particle_type") == "primary":
        y = df["N[#/m3]"].to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
        plot_kw.setdefault("markerfmt", "o")
        plot_kw.setdefault("basefmt", " ")
        markerline, stemlines, _ = ax.stem(x[mask], y[mask], label=label, **plot_kw)
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$N$ [m$^{-3}$]")
        ax.grid(True, which="both", alpha=0.25)
        return ax

    y_col = f"dN/dlog10({size_col})[#/m3]"
    y = df[y_col].to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)

    plot_kw.setdefault("marker", "o")
    plot_kw.setdefault("markersize", 4)
    plot_kw.setdefault("linewidth", 1.5)
    ax.plot(x[mask], y[mask], label=label, **plot_kw)

    ax.set_yscale("log")
    ax.set_xlim(left=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(rf"$dN/d\log_{{10}}({sym})$ [m$^{{-3}}$]")
    ax.grid(True, which="both", alpha=0.25)
    return ax
