import glob
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .postprocessor import PostProcessor

# From Adrien's Python implementation
# This is now operated in C++ with Python as interface only.


def elements_balance(kineticFolder, resultsFolder, elements_list, threshold=0.05):
    """Plot, for each element in `elements_list`, a stackplot of which species carries
    it across an OpenSMOKE run.

    resultsFolder layout picks the mode automatically:
      - a folder of Case*/ subfolders (each holding an Output.xml, e.g. a temperature
        sweep): one steady-state point per Case (last row of its profile), plotted
        against that Case's final temperature.
      - a single run folder (its own Output.xml directly inside, e.g. a Flame1D run):
        every point of that run's own independent variable (time, or the spatial
        coordinate for a flame) is used.

    kineticFolder    -- folder containing kinetics.xml (element composition per species)
    resultsFolder    -- folder containing Case*/Output.xml (a sweep), or a single run's
                         output folder (its own Output.xml directly inside)
    elements_list    -- element symbols to plot, e.g. ["C", "H", "O"]; any case is accepted
                         ("he", "HE", "He" all match); elements absent from every species in
                         the mechanism are skipped with a warning rather than raising.
    threshold        -- a species gets its own band on an element's plot if it holds more
                         than this fraction of that element's total at some point in the
                         run; everything else is lumped into "Others"
    """
    O2_break_threshold = 0.30    # trigger a broken y-axis on the O plot if O2's peak share of total O exceeds this
    O2_break_margin = 0.03       # gap (as a fraction of total O) left on each side of the cut
    O2_break_bottom_fraction = 0.03  # bottom (0-margin) panel's share of the row's vertical space

    # matplotlib's tab20 hues, reordered as the 10 saturated colors followed by their
    # 10 light variants (instead of tab20's default saturated/light interleaving), so
    # the first 10 species assigned stay maximally distinct from each other.
    # Species keep the same slot across all element plots (see assign_colors).
    palette = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
        "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    ]
    others_colour = "silver"  # reserved for the "Others" lumped bucket to make it different from standard grey

    title_fontsize = 14
    legend_fontsize = 9

    def normalize_element(symbol):
        """'HE' / 'he' / 'He' => 'He', so the mechanism's all-caps element names
        (and however the caller capitalized elements_list) always match up."""
        return symbol[0].upper() + symbol[1:].lower() if len(symbol) > 1 else symbol.upper()

    def lump_minor_species(elem_df, threshold):
        """Keep species whose peak fraction of the element total exceeds "threshold";
        sum the rest into an 'Others' column. Returns DataFrame sorted by mean share, descending."""
        total = elem_df.sum(axis=1)
        share = elem_df.div(total.replace(0, np.nan), axis=0).fillna(0)
        keep = share.columns[share.max(axis=0) >= threshold]
        kept = elem_df[keep]
        others = elem_df.drop(columns=keep).sum(axis=1)
        out = kept.copy()
        if others.abs().sum() > 0:
            out['Others'] = others
        return out[out.mean(axis=0).sort_values(ascending=False).index]

    def assign_colors(lumped_by_element, elements):
        """Map each species name to a fixed palette slot, shared across every element's
        plot so the same species always reads as the same color. First seen order,
        scanning elements in 'elements' order; 'Others' always gets others_colour."""
        color_of = {}
        next_slot = 0
        for el in elements:
            for sname in lumped_by_element[el].columns:
                if sname == 'Others' or sname in color_of:
                    continue
                color_of[sname] = palette[next_slot % len(palette)]
                next_slot += 1
        if next_slot > len(palette):
            print(f"[warn] {next_slot} distinct species across plots, only {len(palette)} "
                  f"palette colors some species share a color", file=sys.stderr)
        color_of['Others'] = others_colour
        return color_of

    def add_percent_axis(ax, total_ref, label=False):
        secax = ax.secondary_yaxis(
            'right',
            functions=(lambda y: y / total_ref * 100, lambda p: p / 100 * total_ref)
        )
        if label:
            secax.set_ylabel("% of total")
        return secax

    def draw_break_marks(ax_top, ax_bottom):
        """ Its break an axis in case of a species taking too much space (O2 mostly)
        """
        ax_top.spines.bottom.set_visible(False)
        ax_bottom.spines.top.set_visible(False)
        ax_top.tick_params(bottom=False, labelbottom=False)
        d = 0.5
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=10, linestyle="none",
                      color='k', mec='k', mew=1, clip_on=False)
        ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **kwargs)
        ax_bottom.plot([0, 1], [1, 1], transform=ax_bottom.transAxes, **kwargs)

    def draw_plot(elem_dfs, x_label, elements, threshold):
        x = next(iter(elem_dfs.values())).index.values
        x_short = x_label.split('[')[0].strip()

        lumped_by_element = {el: lump_minor_species(elem_dfs[el], threshold) for el in elements}
        color_of = assign_colors(lumped_by_element, elements)

        # Decide which elements need a broken axis: an O2 band that eats most of the plot.
        breaks = {}
        for el in elements:
            edf = elem_dfs[el]
            if 'O2' not in edf.columns:
                continue
            total = edf.sum(axis=1)
            share_o2 = (edf['O2'] / total.replace(0, np.nan)).fillna(0)
            if share_o2.max() <= O2_break_threshold:
                continue
            total_ref = total.mean()
            low = O2_break_margin * total_ref
            high = (share_o2.min() - O2_break_margin) * total_ref
            if high > low:
                breaks[el] = (low, high, total.max(), total_ref)

        fig = plt.figure(figsize=(10.5, 3.6 * len(elements)))
        # Outer grid: one loose-spaced row per element, so titles never crowd the plot
        # above. Elements needing a broken axis get a tight nested 2-row grid of their own.
        outer = fig.add_gridspec(len(elements), 1, hspace=0.55, left=0.08, right=0.78, top=0.95, bottom=0.07)

        bottom_ax = None
        for i, el in enumerate(elements):
            lumped = lumped_by_element[el]
            colors = [color_of[c] for c in lumped.columns]

            if el in breaks:
                low, high, top, total_ref = breaks[el]
                inner = outer[i].subgridspec(2, 1, height_ratios=(1.0 - O2_break_bottom_fraction,
                                                                    O2_break_bottom_fraction), hspace=0.08)
                ax_top = fig.add_subplot(inner[0])
                ax_bottom = fig.add_subplot(inner[1], sharex=ax_top)
                ax_top.stackplot(x, [lumped[c].values for c in lumped.columns], colors=colors,
                                  labels=lumped.columns, alpha=0.9)
                ax_bottom.stackplot(x, [lumped[c].values for c in lumped.columns], colors=colors, alpha=0.9)
                ax_top.set_ylim(high, top * 1.02)
                ax_bottom.set_ylim(0, low)
                draw_break_marks(ax_top, ax_bottom)
                add_percent_axis(ax_top, total_ref, label=True)
                secax_bottom = add_percent_axis(ax_bottom, total_ref)
                # bottom panel is a thin sliver (just the 0-margin strip) — too little
                # room for readable moles tick labels there, so the left axis just gets
                # tick marks and the graduation (0, 3%) lives on the right (%) axis
                margin_pct = O2_break_margin * 100
                ax_bottom.set_yticks([0, low])
                ax_bottom.set_yticklabels([])
                secax_bottom.set_yticks([0, margin_pct])
                ax_top.set_title(f"{el} speciation vs {x_short}  (species >{threshold:.0%} of {el} total; "
                                  f"axis broken {low / total_ref:.0%}-{high / total_ref:.0%})",
                                  fontsize=title_fontsize, pad=10)
                ax_top.set_ylabel(f"{el} moles")
                # stackplot draws the first column at the bottom of the stack; reverse
                # the legend so it reads top-to-bottom in the same order as the stack
                handles, labels = ax_top.get_legend_handles_labels()
                ax_top.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1.15, 0.5),
                              fontsize=legend_fontsize)
                bottom_ax = ax_bottom
            else:
                ax = fig.add_subplot(outer[i], sharex=bottom_ax)
                ax.stackplot(x, [lumped[c].values for c in lumped.columns], colors=colors, labels=lumped.columns,
                             alpha=0.9)
                total_ref = elem_dfs[el].sum(axis=1).mean()
                add_percent_axis(ax, total_ref, label=True)
                ax.set_ylabel(f"{el} moles")
                ax.set_title(f"{el} speciation vs {x_short}  (species >{threshold:.0%} of {el} total)",
                             fontsize=title_fontsize, pad=10)
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1.15, 0.5),
                          fontsize=legend_fontsize)
                bottom_ax = ax

        bottom_ax.set_xlabel(x_label)

    elements_list = [normalize_element(e) for e in elements_list]

    case_folders = sorted(glob.glob(os.path.join(resultsFolder, "Case*")))
    is_sweep = bool(case_folders)
    # One PostProcessor built once for the whole mechanism; the sweep branch below
    # repoints it at each Case's Output.xml via updateOutput() instead of re-parsing
    # kinetics.xml per Case - the motivating use case for updateOutput().
    pp = PostProcessor(kineticFolder, case_folders[0] if is_sweep else resultsFolder)

    available_elements = []
    for el in elements_list:
        try:
            pp.ElementMolesBySpecies(el)
        except ValueError:
            warnings.warn(
                f"element '{el}' is not present in any species of this mechanism, skipped"
            )
            continue
        available_elements.append(el)
    elements_list = available_elements
    if not elements_list:
        raise ValueError("None of the requested elements are present in this mechanism")

    if is_sweep:
        rows = {el: {} for el in elements_list}
        for i, case_folder in enumerate(case_folders):
            if i > 0:
                pp.updateOutput(case_folder)
            _, temperature = pp.getTemperatureProfile()
            T_last = float(temperature[-1])
            for el in elements_list:
                rows[el][T_last] = pp.ElementMolesBySpecies(el).iloc[-1]
        elem_dfs = {
            el: pd.DataFrame.from_dict(case_rows, orient="index").sort_index().fillna(0.0)
            for el, case_rows in rows.items()
        }
        x_label = "Temperature [K]"
    else:
        elem_dfs = {el: pp.ElementMolesBySpecies(el) for el in elements_list}
        x_label = "Independent variable"

    draw_plot(elem_dfs, x_label, elements_list, threshold)
    plt.show()
