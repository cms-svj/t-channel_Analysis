#!/usr/bin/env python3

# Combines what used to be plot_variations.py (N_SVJ tagged-jet histograms +
# raw Lund weight distributions) and systematics_compare_MET_stitched.py
# (preselection-level MET spectra split by n_SVJ category). Both parts read
# from the same input skims (same YEARS/SAMPLE, same TAGGER_CONFIGS, same
# tagging functions and variation list), so they're merged here into one
# script with one shared event-loading setup and one combined ROOT output.
# The old two scripts are left in place for now (not deleted).

import os
import argparse
import uproot
import awkward as ak
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullLocator, AutoMinorLocator

# Set the style for matplotlib
plt.style.use(hep.style.CMS)

plt.rcParams['legend.title_fontsize'] = 'xx-small'
plt.rcParams["axes.labelsize"] = 22
plt.rcParams["xtick.labelsize"] = 18
plt.rcParams["ytick.labelsize"] = 18
plt.rcParams["legend.fontsize"] = 16

# Combined ROOT output for every histogram plotted below. Recreated fresh on
# every run (rather than uproot.update(), which appends a new ROOT "cycle"
# for each overwritten key instead of truly replacing it) so the file never
# accumulates stale cycles across repeated runs.
ROOT_OUTPUT_PATH = "histograms.root"


def open_output_root():
    return uproot.recreate(ROOT_OUTPUT_PATH)


# The preselection MET spectra split by n_SVJ category aren't needed for
# now — set True to also generate those plots/histograms. This does NOT
# control the tagger-agnostic "All" preselection MET histogram (that one
# isn't split by n_SVJ at all, so it's generated unconditionally below).
MAKE_NSVJ_PLOTS = False

# The raw Lund weight distribution plots (lund_weights_all_variations.png +
# one plot per variation) aren't needed for now either — set True to
# generate those plots/histograms.
MAKE_LUND_WEIGHT_PLOTS = False


# CMS Petroff-6 color-blind-friendly palette; nominal/uncorrected stays black,
# nominal Lund-corrected is blue (matches "Lund Nominal" in COLOR_MAP below);
# Up/Down variations share color+style. Petroff purple ("#964a8b") was too
# close to the violet next to it (the last two variations), so it's swapped
# for green.
NOMINAL_COLOR = "black"
NOMINAL_LUND_COLOR = "#118AB2"
VARIATION_COLORS = ['#3f90da', '#ffa90e', '#94a4a2', '#832db6', '#a96b59', '#e76300', '#b9ac70', '#717581', '#92dad']
VARIATION_LINESTYLE = "--"


def line_legend_handle(color, linestyle='-'):
    """A plain line legend swatch, in place of matplotlib's default step-histogram box."""
    return Line2D([], [], color=color, linestyle=linestyle, linewidth=2)


def get_number_of_tagged_svjs(df, wp):

    def __is_tagged(
        wnae_pt_0_200,
        wnae_pt_200_300,
        wnae_pt_300_400,
        wnae_pt_400_500,
        wnae_pt_500_inf,
        pt,
        wp,
    ):
        # High stat WP 17/11/2024 - 2018 but should work for all years
        if wp == 10:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 26.259)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 20.065)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 23.472)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 22.599)
               + (pt >= 500) * (wnae_pt_500_inf > 18.747)
            )
        elif wp == 20:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 25.156)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 18.284)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 20.383)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.941)
               + (pt >= 500) * (wnae_pt_500_inf > 16.370)
            )
        elif wp == 25:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 24.843)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 17.802)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 19.525)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.739)
               + (pt >= 500) * (wnae_pt_500_inf > 15.925)
            )
        elif wp == 30:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 24.594)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 17.440)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 18.776)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.558)
               + (pt >= 500) * (wnae_pt_500_inf > 15.468)
            )
        elif wp == 35:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 24.373)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 17.107)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 18.100)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.389)
               + (pt >= 500) * (wnae_pt_500_inf > 14.981)
            )
        elif wp == 40:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 24.166)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 16.781)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 17.493)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.218)
               + (pt >= 500) * (wnae_pt_500_inf > 14.526)
            )
        elif wp == 45:
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 23.962)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 16.459)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 16.883)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 21.041)
               + (pt >= 500) * (wnae_pt_500_inf > 14.079)
            )
        #
        # High stat WP 17/11/2024 - 2018 data
        #
        elif wp == "20_data":
            is_tagged = (
               (pt < 200) * (wnae_pt_0_200 > 24.974)
               + (pt >= 200) * (pt < 300) * (wnae_pt_200_300 > 18.445)
               + (pt >= 300) * (pt < 400) * (wnae_pt_300_400 > 21.031)
               + (pt >= 400) * (pt < 500) * (wnae_pt_400_500 > 22.131)
               + (pt >= 500) * (wnae_pt_500_inf > 16.535)
            )
        else:
            print(f"Non valid working point: {wp}")
            exit(1)

        # is_tagged = is_tagged.astype(int)
        is_tagged = is_tagged.to_numpy().astype(int)

        return is_tagged

    is_good_jet = df["JetsAK8_isGood"].array(library="ak")
    pt = ak.pad_none(df["JetsAK8_/.fPt"].array(library="ak")[is_good_jet], 4)

    scores_pt_0_200 = ak.pad_none(df["JetsAK8_WNAEPt0To200Loss"].array(library="ak")[is_good_jet], 4)
    scores_pt_0_200 = ak.fill_none(scores_pt_0_200, -1.)
    scores_pt_200_300 = ak.pad_none(df["JetsAK8_WNAEPt200To300Loss"].array(library="ak")[is_good_jet], 4)
    scores_pt_300_400 = ak.pad_none(df["JetsAK8_WNAEPt300To400Loss"].array(library="ak")[is_good_jet], 4)
    scores_pt_400_500 = ak.pad_none(df["JetsAK8_WNAEPt400To500Loss"].array(library="ak")[is_good_jet], 4)
    scores_pt_500_inf = ak.pad_none(df["JetsAK8_WNAEPt500ToInfLoss"].array(library="ak")[is_good_jet], 4)

    j0_wnae_pt_0_200 = scores_pt_0_200[:, 0]
    j1_wnae_pt_0_200 = scores_pt_0_200[:, 1]
    j2_wnae_pt_0_200 = scores_pt_0_200[:, 2]
    j3_wnae_pt_0_200 = scores_pt_0_200[:, 3]

    j0_wnae_pt_200_300 = scores_pt_200_300[:, 0]
    j1_wnae_pt_200_300 = scores_pt_200_300[:, 1]
    j2_wnae_pt_200_300 = scores_pt_200_300[:, 2]
    j3_wnae_pt_200_300 = scores_pt_200_300[:, 3]

    j0_wnae_pt_300_400 = scores_pt_300_400[:, 0]
    j1_wnae_pt_300_400 = scores_pt_300_400[:, 1]
    j2_wnae_pt_300_400 = scores_pt_300_400[:, 2]
    j3_wnae_pt_300_400 = scores_pt_300_400[:, 3]

    j0_wnae_pt_400_500 = scores_pt_400_500[:, 0]
    j1_wnae_pt_400_500 = scores_pt_400_500[:, 1]
    j2_wnae_pt_400_500 = scores_pt_400_500[:, 2]
    j3_wnae_pt_400_500 = scores_pt_400_500[:, 3]

    j0_wnae_pt_500_inf = scores_pt_500_inf[:, 0]
    j1_wnae_pt_500_inf = scores_pt_500_inf[:, 1]
    j2_wnae_pt_500_inf = scores_pt_500_inf[:, 2]
    j3_wnae_pt_500_inf = scores_pt_500_inf[:, 3]

    j0_pt = pt[:, 0]
    j1_pt = pt[:, 1]
    j2_pt = pt[:, 2]
    j3_pt = pt[:, 3]

    svj0 = __is_tagged(
        j0_wnae_pt_0_200,
        j0_wnae_pt_200_300,
        j0_wnae_pt_300_400,
        j0_wnae_pt_400_500,
        j0_wnae_pt_500_inf,
        j0_pt,
        wp,
    )

    svj1 = __is_tagged(
        j1_wnae_pt_0_200,
        j1_wnae_pt_200_300,
        j1_wnae_pt_300_400,
        j1_wnae_pt_400_500,
        j1_wnae_pt_500_inf,
        j1_pt,
        wp,
    )

    svj2 = __is_tagged(
        j2_wnae_pt_0_200,
        j2_wnae_pt_200_300,
        j2_wnae_pt_300_400,
        j2_wnae_pt_400_500,
        j2_wnae_pt_500_inf,
        j2_pt,
        wp,
    )
    svj3 = __is_tagged(
        j3_wnae_pt_0_200,
        j3_wnae_pt_200_300,
        j3_wnae_pt_300_400,
        j3_wnae_pt_400_500,
        j3_wnae_pt_500_inf,
        j3_pt,
        wp,
    )

    n_svjs = svj0 + svj1 + svj2 + svj3

    return n_svjs


def get_number_of_tagged_svjs_pnet(df, wp):
    # Unlike WNAE, PNet has a single score per jet with no pt-dependent
    # working points, so tagging is just a flat threshold on the leading 4 jets.
    is_good_jet = df["JetsAK8_isGood"].array(library="ak")
    score = ak.pad_none(df["JetsAK8_pNetJetTaggerScore"].array(library="ak")[is_good_jet], 4, clip=True)
    score = ak.fill_none(score, -1.)
    score = score.to_numpy()
    return (score > wp).sum(axis=1).astype(int)


YEARS = ["2016", "2017", "2018"]
SAMPLE = "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1"

# Each tagger pulls from a different skim path (PNet only lives in the
# t_channel_pre_selection skim, WNAE only in signals_WNAE — see
# project_lundplotspaper_pnet_plot memory) and has its own tagging function,
# working point, and labels. Shared by both the tagged-SVJ plots (xlabel) and
# the preselection MET plots (cat_superscript/corner_label) below.
TAGGER_CONFIGS = [
    {
        "name": "wnae",
        "subdir_template": "signals_WNAE/{year}/nominal",
        "tag_func": get_number_of_tagged_svjs,
        "wp": 20,
        "xlabel": r'$n_{SVJ}^{WNAE}$',
        "cat_superscript": "WNAE",
        "corner_label": "WNAE",
        "suffix": "",
    },
    {
        "name": "pnet",
        "subdir_template": "{year}/t_channel_pre_selection/nominal",
        "tag_func": get_number_of_tagged_svjs_pnet,
        "wp": 0.9,
        "xlabel": r'$n_{SVJ}^{PNet}$',
        "cat_superscript": "PNet",
        "corner_label": "ParticleNet",
        "suffix": "_pnet",
    },
]


def year_file(year, subdir_template):
    return (
        f"/eos/uscms/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/"
        f"{subdir_template.format(year=year)}/{SAMPLE}/part-0.root"
    )


def get_variation_label(variation_name):
    return (
        variation_name.replace('weight', ' ')
        .replace('sys', 'correction systematics')
        .replace('stat', 'correction statistical')
        .replace('pt', 'correction $p_{T}$ extrapolation')
        .replace('distortion', 'correction distortion')
    )


multiplicative_variations = [
    # 'lundWeight_nom',
    #'lundWeightProngs',
    'lundWeightSys',
    'lundWeightDistortion',
    #'lundWeightUnclust',
    #'lundWeightBquark',
    # 'nlundWeight_pt_vars',
    # 'lundWeight_pt_vars',
    # 'nlundWeight_stat_vars',
    # 'lundWeight_stat_vars',
    ]

alternative_variations = [
    'lundWeightPt',
    'lundWeightStat',
]

VARIATIONS = multiplicative_variations + alternative_variations

# Paired subsets for the tagged-SVJ plots: Nominal + Lund correction plus
# just these two systematics shown together (as opposed to all four at once).
TAGGED_SVJ_GROUPS = {
    "SysDistortion": ["lundWeightSys", "lundWeightDistortion"],
    "PtStat": ["lundWeightPt", "lundWeightStat"],
}


# ---------------------------------------------------------------------------
# Tagged-SVJ (N_SVJ) histograms + raw Lund weight distributions
# (formerly plot_variations.py)
# ---------------------------------------------------------------------------

def load_tagged_svj_data(tagger_cfg):
    # Each year's cutflow normalizes that year's own events (different Initial /
    # InitialLund* counts per year), so every weight must be computed per-year
    # BEFORE concatenating across years — never mix cutflow factors across years.
    n_tagged_jets_parts = []
    pu_weights_parts = []
    weights_nominal_lund_parts = []
    variation_weighted_parts = {var: {'up': [], 'down': []} for var in VARIATIONS}

    for year in YEARS:
        file_path = year_file(year, tagger_cfg["subdir_template"])
        year_data = uproot.open(file_path)["Events"]
        year_cutflow = uproot.open(file_path)["CutFlow"]

        year_pu_weights = year_data["puWeight"].array(library="np")
        year_lund_weights = year_data["lundWeightNom"].array(library="np")
        year_initial = year_cutflow["Initial"].array(library="np")[0]

        n_tagged_jets_parts.append(tagger_cfg["tag_func"](year_data, tagger_cfg["wp"]))
        pu_weights_parts.append(year_pu_weights)
        weights_nominal_lund_parts.append(
            year_pu_weights * year_lund_weights
            / year_cutflow["InitialLundNominal"].array(library="np")[0] * year_initial
        )

        for var in VARIATIONS:
            variation_weighted_parts[var]['up'].append(
                year_pu_weights * year_data[f"{var}Up"].array(library="np")
                / year_cutflow[f"InitialLund{var.capitalize()}Up"].array(library="np")[0] * year_initial
            )
            variation_weighted_parts[var]['down'].append(
                year_pu_weights * year_data[f"{var}Down"].array(library="np")
                / year_cutflow[f"InitialLund{var.capitalize()}Down"].array(library="np")[0] * year_initial
            )

    n_tagged_jets = np.concatenate(n_tagged_jets_parts)
    pu_weights = np.concatenate(pu_weights_parts)
    weights_nominal_lund = np.concatenate(weights_nominal_lund_parts)
    variation_weights = {
        var: {
            'up': np.concatenate(variation_weighted_parts[var]['up']),
            'down': np.concatenate(variation_weighted_parts[var]['down']),
        }
        for var in VARIATIONS
    }
    return n_tagged_jets, pu_weights, weights_nominal_lund, variation_weights


def write_tagged_svj_histograms(root_file, tagger_cfg, n_tagged_jets, pu_weights, weights_nominal_lund,
                                 variation_weights):
    # Raw (un-normalized) histograms only — the normalized version is just
    # these divided by their sum, so it doesn't need its own copy in the
    # ROOT file. Written directly from the loaded data rather than from
    # make_tagged_svj_plot, since not every plot variant (e.g. the paired
    # subsets) covers the full set of variations.
    edges = np.arange(0, 5, dtype=float)
    prefix = f"tagged_svj/{tagger_cfg['name']}"
    nominal_counts, _ = np.histogram(n_tagged_jets, bins=edges, weights=pu_weights)
    lund_nominal_counts, _ = np.histogram(n_tagged_jets, bins=edges, weights=weights_nominal_lund)
    root_file[f"{prefix}/Nominal"] = (nominal_counts, edges)
    root_file[f"{prefix}/LundNominal"] = (lund_nominal_counts, edges)
    for var in VARIATIONS:
        up_counts, _ = np.histogram(n_tagged_jets, bins=edges, weights=variation_weights[var]['up'])
        down_counts, _ = np.histogram(n_tagged_jets, bins=edges, weights=variation_weights[var]['down'])
        root_file[f"{prefix}/{var}Up"] = (up_counts, edges)
        root_file[f"{prefix}/{var}Down"] = (down_counts, edges)


def make_tagged_svj_plot(tagger_cfg, n_tagged_jets, pu_weights, weights_nominal_lund, variation_weights,
                          normalize=True, variations=None, tag=None):
    # normalize only rescales the already-loaded histograms (density=True in
    # plt.hist below) — it doesn't need its own data load, so callers should
    # load once via load_tagged_svj_data() and reuse it for both variants.
    # Styled to match plot_systematic_comparison (the preselection MET
    # plots): same colors/labels (COLOR_MAP/DISPLAY_NAME/legend_label_for),
    # shaded uncertainty bands instead of separate Up/Down lines, and the
    # same legend/corner-label box conventions.
    # variations restricts which systematics are drawn (default: all).
    if variations is None:
        variations = VARIATIONS

    edges = np.arange(0, 5, dtype=float)

    fig, ax = plt.subplots()

    nominal_counts, _ = np.histogram(n_tagged_jets, bins=edges, weights=pu_weights, density=normalize)
    lund_nominal_counts, _ = np.histogram(n_tagged_jets, bins=edges, weights=weights_nominal_lund,
                                           density=normalize)

    ax.stairs(nominal_counts, edges, linewidth=2.5, color=COLOR_MAP["Nominal"], label="Nominal")
    ax.stairs(lund_nominal_counts, edges, linewidth=2, color=COLOR_MAP["Lund Nominal"], label="Lund Nominal")

    for var in variations:
        up_counts, _ = np.histogram(n_tagged_jets, bins=edges, weights=variation_weights[var]['up'],
                                     density=normalize)
        down_counts, _ = np.histogram(n_tagged_jets, bins=edges, weights=variation_weights[var]['down'],
                                       density=normalize)

        disp = DISPLAY_NAME[var]
        lower = np.minimum(up_counts, down_counts)
        upper = np.maximum(up_counts, down_counts)
        ax.fill_between(
            edges,
            np.append(lower, lower[-1]),
            np.append(upper, upper[-1]),
            step="post",
            color=COLOR_MAP[f"{disp}Up"],
            alpha=0.25,
            linewidth=0,
            label=legend_label_for(f"{disp}Up"),
        )

    ax.set_xlabel(tagger_cfg["xlabel"])
    ax.set_ylabel('Arbitrary Units' if normalize else 'Events')
    ax.set_xticks(range(0, 5))
    # Widen x-range to open a blank column to the right of the data (bins only go
    # to 4) for the legend, instead of stacking vertical headroom above the data
    ax.set_xlim(-0.5, 6)
    _, ymax = ax.get_ylim()
    ax.set_ylim(0, ymax * 1.15)

    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(
        handles, labels,
        fontsize=17,
        ncol=max(2, len(labels) // 3),
        loc="upper right",
        frameon=True,
        edgecolor="white",
        framealpha=1,
    )
    legend.set_zorder(20)

    bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="white", alpha=1)
    ax.text(
        0.02, 0.96, tagger_cfg["corner_label"],
        transform=ax.transAxes,
        va="top", ha="left",
        fontsize=21,
        bbox=bbox_props
    )
    # Stacked (not one long line) and kept in the blank column to the right
    # of the data (bins only go to 4, xlim extends to 6) so it doesn't run
    # across the tall bins in the middle of the plot — the legend itself
    # spans that same column at the top, so x is fixed at the edge of the
    # data rather than aligned to the legend's own (much wider) left edge.
    # y tracks the legend's bottom edge so it sits right below it regardless
    # of how many rows the legend ends up with.
    fig.canvas.draw()
    legend_bbox = ax.transAxes.inverted().transform_bbox(legend.get_window_extent(fig.canvas.get_renderer()))
    ax.text(
        0.70, legend_bbox.y0 - 0.03,
        "$m_{\\Phi} = 2000$ GeV\n$m_{Dark} = 20$ GeV\n$r_{inv} = 0.3$\n$\\lambda = 1$",
        transform=ax.transAxes,
        fontsize=17,
        ha='left',
        va='top',
        linespacing=1.6,
    )

    plt.tight_layout()
    hep.cms.label("", data=False, ax=ax)
    norm_part = "_normalized" if normalize else ""
    tag_part = f"_{tag}" if tag else ""
    fig.savefig(f"lund_weights_plots/tagged_svjs_all_variations{norm_part}{tagger_cfg['suffix']}{tag_part}.png",
                bbox_inches="tight")
    plt.close(fig)


def load_lund_weights_data():
    # This data only depends on lundWeightNom/lundWeight*Up/Down, never on a
    # tagging score, so it's loaded once from a single file set rather than
    # once per tagger (the two skims share the same underlying events and
    # weights for this sample, so per-tagger copies would be identical).
    lund_weights_parts = []
    variation_raw_parts = {var: {'up': [], 'down': []} for var in VARIATIONS}

    for year in YEARS:
        file_path = year_file(year, TAGGER_CONFIGS[0]["subdir_template"])
        year_data = uproot.open(file_path)["Events"]

        lund_weights_parts.append(year_data["lundWeightNom"].array(library="np"))
        for var in VARIATIONS:
            variation_raw_parts[var]['up'].append(year_data[f"{var}Up"].array(library="np"))
            variation_raw_parts[var]['down'].append(year_data[f"{var}Down"].array(library="np"))

    lund_weights = np.concatenate(lund_weights_parts)
    variation_raw = {
        var: {
            'up': np.concatenate(variation_raw_parts[var]['up']),
            'down': np.concatenate(variation_raw_parts[var]['down']),
        }
        for var in VARIATIONS
    }
    lund_weights_range = (0, np.quantile(lund_weights, 0.95))
    return lund_weights, variation_raw, lund_weights_range


def make_lund_weights_plot(lund_weights, variation_raw, lund_weights_range, root_file=None):
    plt.figure()
    nominal_counts, edges, _ = plt.hist(
        lund_weights,
        bins=100,
        range=lund_weights_range,
        label='Lund Weights',
        histtype='step',
        linewidth=2,
        edgecolor=NOMINAL_COLOR,
        linestyle='-',
    )
    variation_labels = []
    variation_hist_data = {}
    for var, color in zip(VARIATIONS, VARIATION_COLORS):
        label = get_variation_label(var.replace("lundWeight_", "").replace("_", " ").title())
        up_counts, _, _ = plt.hist(
            variation_raw[var]['up'],
            bins=100,
            range=lund_weights_range,
            label=label,
            histtype='step',
            linewidth=2,
            edgecolor=color,
            linestyle=VARIATION_LINESTYLE,
        )
        down_counts, _, _ = plt.hist(
            variation_raw[var]['down'],
            bins=100,
            range=lund_weights_range,
            label=None,
            histtype='step',
            linewidth=2,
            edgecolor=color,
            linestyle=VARIATION_LINESTYLE,
        )
        variation_labels.append((label, color))
        variation_hist_data[var] = (up_counts, down_counts)

    if root_file is not None:
        root_file["lund_weights/LundWeightNom"] = (nominal_counts, edges)
        for var, (up_counts, down_counts) in variation_hist_data.items():
            root_file[f"lund_weights/{var}Up"] = (up_counts, edges)
            root_file[f"lund_weights/{var}Down"] = (down_counts, edges)

    plt.xlabel('Lund Weights')
    plt.ylabel('Number of Events')
    # Widen x-range to open a blank column to the right of the data (the histograms
    # stop at lund_weights_range[1]) for the legend, instead of vertical headroom
    plt.xlim(0, lund_weights_range[1] * 1.1)
    _, ymax = plt.ylim()
    plt.ylim(0, ymax * 1.1)
    legend_handles = [line_legend_handle(NOMINAL_COLOR, '-')] + [
        line_legend_handle(color, VARIATION_LINESTYLE) for _, color in variation_labels
    ]
    legend_labels = ['Lund Weights'] + [label for label, _ in variation_labels]
    plt.legend(
        handles=legend_handles,
        labels=legend_labels,
        fontsize=17,
        ncol=1,
        loc='upper left',
        bbox_to_anchor=(0.45, 0.95),
        borderaxespad=0,
        framealpha=0.9,
    )
    hep.cms.label(data=False, loc=0, fontsize=20)
    plt.savefig("lund_weights_plots/lund_weights_all_variations.png", bbox_inches="tight")
    plt.close()


def make_lund_weight_variation_plot(var, color, lund_weights, variation_raw, lund_weights_range):
    label = get_variation_label(var.replace("lundWeight_", "").replace("_", " ").title())

    plt.figure()
    plt.hist(
        lund_weights,
        bins=100,
        range=lund_weights_range,
        label='Lund Weights',
        histtype='step',
        linewidth=2,
        edgecolor=NOMINAL_COLOR,
        linestyle='-',
    )
    plt.hist(
        variation_raw[var]['up'],
        bins=100,
        range=lund_weights_range,
        label=label,
        histtype='step',
        linewidth=2,
        edgecolor=color,
        linestyle=VARIATION_LINESTYLE,
    )
    plt.hist(
        variation_raw[var]['down'],
        bins=100,
        range=lund_weights_range,
        label=None,
        histtype='step',
        linewidth=2,
        edgecolor=color,
        linestyle=VARIATION_LINESTYLE,
    )
    plt.xlabel('Lund Weights')
    plt.ylabel('Number of Events')
    plt.xlim(0, lund_weights_range[1] * 1.1)
    _, ymax = plt.ylim()
    plt.ylim(0, ymax * 1.1)
    legend_handles = [
        line_legend_handle(NOMINAL_COLOR, '-'),
        line_legend_handle(color, VARIATION_LINESTYLE),
    ]
    legend_labels = ['Lund Weights', label]
    plt.legend(
        handles=legend_handles,
        labels=legend_labels,
        fontsize=17,
        ncol=1,
        loc='upper left',
        bbox_to_anchor=(0.45, 0.95),
        borderaxespad=0,
        framealpha=0.9,
    )
    hep.cms.label(data=False, loc=0, fontsize=20)
    slug = var.replace('lundWeight', '').lower()
    plt.savefig(f"lund_weights_plots/lund_weights_{slug}.png", bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Preselection-level MET spectra, split by n_SVJ category
# (formerly systematics_compare_MET_stitched.py)
# ---------------------------------------------------------------------------

DISPLAY_NAME = {
    "lundWeightStat": "Lund Stat",
    "lundWeightSys": "Lund Sys",
    "lundWeightPt": "Lund pT",
    "lundWeightDistortion": "Lund Distortion",
}

# Reverse of DISPLAY_NAME, used to give ROOT histogram keys the same
# branch-name-based convention as the tagged_svj/lund_weights histograms
# above (e.g. "Lund SysUp" -> "lundWeightSysUp") instead of the
# display-only "Lund Sys" style names.
KEY_NAME = {"Nominal": "Nominal", "Lund Nominal": "LundNominal"}
for _var, _disp in DISPLAY_NAME.items():
    KEY_NAME[f"{_disp}Up"] = f"{_var}Up"
    KEY_NAME[f"{_disp}Down"] = f"{_var}Down"

CATEGORIES = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]

# Placeholder binning — not yet reconciled with the ABCD framework's own MET
# binning, just enough to look at shape differences at preselection level.
MET_EDGES = np.arange(200.0, 1250.0, 50.0)
N_BINS = len(MET_EDGES) - 1

COLOR_MAP = {
    "Nominal": "black",
    "Lund Nominal": "#118AB2",
    "Lund StatUp": "#5790FC",
    "Lund StatDown": "#5790FC",
    "Lund SysUp": "#E42536",
    "Lund SysDown": "#E42536",
    "Lund pTUp": "#F89C20",
    "Lund pTDown": "#F89C20",
    "Lund DistortionUp": "#7A21DD",
    "Lund DistortionDown": "#7A21DD",
}

LINESTYLE_MAP = {name: "-" for name in COLOR_MAP}

LEGEND_LABEL_OVERRIDES = {
    "Lund Stat": "Lund statistical",
    "Lund Sys": "Lund systematic",
    "Lund pT": r"Lund $p_{T}$ extrapolation",
    "Lund Distortion": "Lund distortion",
}

# Paired subsets, mirroring TAGGED_SVJ_GROUPS: Nominal + Lund correction plus
# two systematics shown together (as bands) instead of just one at a time.
SYSTEMATIC_PAIR_GROUPS = {
    "SysDistortion": ["Lund SysUp", "Lund SysDown", "Lund DistortionUp", "Lund DistortionDown"],
    "PtStat": ["Lund pTUp", "Lund pTDown", "Lund StatUp", "Lund StatDown"],
}

LUND_CORRECTION_PT_LABELS = {
    1: r"$200 < p_{\mathrm{T}} < 300$ GeV",
    2: r"$300 < p_{\mathrm{T}} < 400$ GeV",
    3: r"$400 < p_{\mathrm{T}} < 500$ GeV",
    4: r"$500 < p_{\mathrm{T}} < 600$ GeV",
    5: r"$600 < p_{\mathrm{T}} < 700$ GeV",
    6: r"$p_{\mathrm{T}} > 700$ GeV",
}


def draw_order(names):
    return [n for n in names if n not in ("Nominal", "Lund Nominal")] + \
           [n for n in ("Nominal", "Lund Nominal") if n in names]


def legend_label_for(name):
    if name.endswith("Down"):
        return None
    if name.endswith("Up"):
        base = name[:-2]
        return LEGEND_LABEL_OVERRIDES.get(base, base)
    return name


def load_events(subdir_template, tagger_cfg=None):
    # tagger_cfg is only needed to compute n_tagged for the per-category
    # (split) plots. The "All"-category plots don't mask on n_tagged at all,
    # so they can skip tagging entirely and just load MET/weights once from
    # a single canonical skim (see project_lundplotspaper_pnet_plot memory:
    # the WNAE and PNet skims contain the same underlying events/weights,
    # just with different tagger branches attached).
    met_parts = []
    n_tagged_parts = [] if tagger_cfg is not None else None
    weight_parts = {"Nominal": [], "Lund Nominal": []}
    for var in VARIATIONS:
        weight_parts[f"{DISPLAY_NAME[var]}Up"] = []
        weight_parts[f"{DISPLAY_NAME[var]}Down"] = []

    for year in YEARS:
        file_path = year_file(year, subdir_template)
        data = uproot.open(file_path)["Events"]
        cutflow = uproot.open(file_path)["CutFlow"]

        met = data["MET"].array(library="np")
        pu_weight = data["puWeight"].array(library="np")
        lund_nom = data["lundWeightNom"].array(library="np")
        initial = cutflow["Initial"].array(library="np")[0]
        initial_lund_nom = cutflow["InitialLundNominal"].array(library="np")[0]

        met_parts.append(met)
        if tagger_cfg is not None:
            n_tagged_parts.append(tagger_cfg["tag_func"](data, tagger_cfg["wp"]))
        weight_parts["Nominal"].append(pu_weight)
        weight_parts["Lund Nominal"].append(
            pu_weight * lund_nom / initial_lund_nom * initial
        )

        for var in VARIATIONS:
            for direction in ("Up", "Down"):
                raw = data[f"{var}{direction}"].array(library="np")
                initial_var = cutflow[f"InitialLund{var.capitalize()}{direction}"].array(library="np")[0]
                weight_parts[f"{DISPLAY_NAME[var]}{direction}"].append(
                    pu_weight * raw / initial_var * initial
                )

    met = np.concatenate(met_parts)
    n_tagged = np.concatenate(n_tagged_parts) if tagger_cfg is not None else None
    weights = {name: np.concatenate(parts) for name, parts in weight_parts.items()}
    return met, n_tagged, weights


def category_mask(n_tagged, category):
    if category == "3PSVJ":
        return n_tagged >= 3
    return n_tagged == int(category[0])


def weighted_hist(values, weights, edges):
    counts, _ = np.histogram(values, bins=edges, weights=weights)
    sumw2, _ = np.histogram(values, bins=edges, weights=weights ** 2)
    return counts, np.sqrt(sumw2)


def write_root_histograms(root_file, prefix, met, weights, mask=None):
    m = met if mask is None else met[mask]
    for name, w in weights.items():
        wm = w if mask is None else w[mask]
        counts, _ = weighted_hist(m, wm, MET_EDGES)
        key = KEY_NAME.get(name, name)
        root_file[f"{prefix}/{key}"] = (counts, MET_EDGES)


def build_xlabels(n_segments, skip=4):
    labels = []
    for _ in range(n_segments):
        for i in range(N_BINS):
            labels.append(str(int(MET_EDGES[i])) if i % skip == 0 else "")
    labels.append(str(int(MET_EDGES[-1])))
    return labels


def build_axis_edges(n_segments):
    return np.arange(0, n_segments * N_BINS + 1, dtype=float)


def plot_systematic_comparison(met, n_tagged, weights, output_dir, tagger_cfg, with_ratio=True, names=None, tag=None,
                                combine_categories=False, band_systematics=False):
    os.makedirs(output_dir, exist_ok=True)

    if names is None:
        names = list(COLOR_MAP.keys())

    if band_systematics:
        line_names = [n for n in names if n in ("Nominal", "Lund Nominal")]
        band_bases = sorted({n[:-2] for n in names if n.endswith("Up") and n not in ("Nominal", "Lund Nominal")})
    else:
        line_names = names
        band_bases = []

    if combine_categories:
        categories = ["All"]
        masks = {"All": np.ones(len(met), dtype=bool)}
    else:
        categories = CATEGORIES
        masks = {cat: category_mask(n_tagged, cat) for cat in categories}

    counts = {name: {} for name in names}
    errors = {name: {} for name in names}
    for cat in categories:
        mask = masks[cat]
        cat_met = met[mask]
        for name in names:
            c, e = weighted_hist(cat_met, weights[name][mask], MET_EDGES)
            counts[name][cat] = c
            errors[name][cat] = e

    stitched = {name: np.concatenate([counts[name][cat] for cat in categories]) for name in names}

    axis_edges = build_axis_edges(len(categories))
    x_labels = build_xlabels(len(categories))

    # Keep per-bin width roughly constant across category counts (14" was
    # tuned for the 4-category split plot) so the single-category "All" plot
    # isn't blown up to 4x-wide bins; floor it so the legend/CMS label still
    # have room to fit.
    fig_width = max(10.0, 14.0 * len(categories) / len(CATEGORIES))

    if with_ratio:
        fig, (ax, rax) = plt.subplots(
            2, 1, figsize=(fig_width, 10), sharex=True,
            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05}
        )
    else:
        fig, ax = plt.subplots(figsize=(fig_width, 8))
        rax = None

    hmax = max(np.max(stitched[name] + np.concatenate([errors[name][cat] for cat in categories])) for name in names)

    for name in draw_order(line_names):
        ax.stairs(
            stitched[name],
            axis_edges,
            linewidth=2.5 if name == "Nominal" else 2,
            linestyle=LINESTYLE_MAP[name],
            color=COLOR_MAP[name],
            label=legend_label_for(name),
        )

    for base in band_bases:
        lower = np.minimum(stitched[f"{base}Up"], stitched[f"{base}Down"])
        upper = np.maximum(stitched[f"{base}Up"], stitched[f"{base}Down"])
        ax.fill_between(
            axis_edges,
            np.append(lower, lower[-1]),
            np.append(upper, upper[-1]),
            step="post",
            color=COLOR_MAP[f"{base}Up"],
            alpha=0.25,
            linewidth=0,
            label=legend_label_for(f"{base}Up"),
        )

    if with_ratio:
        lund_nominal = stitched.get("Lund Nominal")
        if lund_nominal is not None:
            for name in draw_order(line_names):
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = np.where(lund_nominal != 0, stitched[name] / lund_nominal, 0.0)
                rax.stairs(
                    ratio,
                    axis_edges,
                    linewidth=2.5 if name == "Nominal" else 2,
                    linestyle=LINESTYLE_MAP[name],
                    color=COLOR_MAP[name],
                )

            for base in band_bases:
                lower = np.minimum(stitched[f"{base}Up"], stitched[f"{base}Down"])
                upper = np.maximum(stitched[f"{base}Up"], stitched[f"{base}Down"])
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio_lower = np.where(lund_nominal != 0, lower / lund_nominal, 0.0)
                    ratio_upper = np.where(lund_nominal != 0, upper / lund_nominal, 0.0)
                rax.fill_between(
                    axis_edges,
                    np.append(ratio_lower, ratio_lower[-1]),
                    np.append(ratio_upper, ratio_upper[-1]),
                    step="post",
                    color=COLOR_MAP[f"{base}Up"],
                    alpha=0.25,
                    linewidth=0,
                )

    ax.set_ylabel("Events")
    ax.set_yscale("log")
    ax.set_ylim(0.1, max(1, hmax) * 200)

    ax.yaxis.set_minor_locator(
        LogLocator(base=10, subs=np.arange(2, 10), numticks=100)
    )
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(NullLocator())

    bottom_ax = rax if with_ratio else ax

    if with_ratio:
        ax.tick_params(labelbottom=False)
        rax.axhline(1.0, color="black", ls="--", lw=1)
        rax.set_ylim(0.5, 1.5)
        rax.set_ylabel("Relative Uncertainty", fontsize=14)
        rax.yaxis.set_minor_locator(AutoMinorLocator(5))
        rax.minorticks_on()
        rax.xaxis.set_minor_locator(NullLocator())

    bottom_ax.set_xlabel(r"$p_{\mathrm{T}}^{\mathrm{miss}}$ [GeV]")

    cat_superscript = tagger_cfg["cat_superscript"]
    cat_lbl_map = {
        "0SVJ": f"$n_{{SVJ}}^{{{cat_superscript}}}=0$",
        "1SVJ": f"$n_{{SVJ}}^{{{cat_superscript}}}=1$",
        "2SVJ": f"$n_{{SVJ}}^{{{cat_superscript}}}=2$",
        "3PSVJ": f"$n_{{SVJ}}^{{{cat_superscript}}}\\geq3$"
    }

    for i, cat in enumerate(categories):
        idx = i * N_BINS
        if i > 0:
            ax.axvline(axis_edges[idx], color="black", linestyle=":", alpha=0.7, zorder=10)
            if with_ratio:
                rax.axvline(axis_edges[idx], color="black", linestyle=":", alpha=0.7, zorder=10)

        if cat in cat_lbl_map:
            local = axis_edges[idx:idx + N_BINS + 1]
            ax.text(
                (local[0] + local[-1]) / 2,
                0.72,
                cat_lbl_map[cat],
                transform=ax.get_xaxis_transform(),
                ha="center", va="top",
                fontsize=17, fontweight="bold"
            )

    bottom_ax.set_xlim(axis_edges[0], axis_edges[-1])
    bottom_ax.set_xticks(axis_edges)
    bottom_ax.set_xticklabels(x_labels, fontsize=21)

    tick_axes = [ax, rax] if with_ratio else [ax]
    for tick_ax in tick_axes:
        for label, tick in zip(x_labels, tick_ax.xaxis.get_major_ticks()):
            length = 16 if label != "" else 8
            tick.tick1line.set_markersize(length)
            tick.tick2line.set_markersize(length)

    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(
        handles, labels,
        fontsize=17,
        ncol=max(2, len(labels) // 3),
        loc="upper right",
        frameon=True,
        edgecolor="white",
        framealpha=1
    )
    legend.set_zorder(20)

    bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="white", alpha=1)
    ax.text(
        0.02, 0.96, tagger_cfg["corner_label"],
        transform=ax.transAxes,
        va="top", ha="left",
        fontsize=21,
        bbox=bbox_props
    )

    # Signal-parameter text, same convention as make_tagged_svj_plot: stacked
    # and positioned from the legend's own bounding box (kept at the edge of
    # the data rather than the legend's own left edge, which spans too wide)
    # so it tracks the legend regardless of how many rows it ends up with.
    # Note: with a 4-category split plot (MAKE_NSVJ_PLOTS=True) and the full
    # systematics legend (3 rows), this can land close to the category
    # labels (fixed at y=0.72) — not an issue for the current default output
    # (combine_categories=True only), but worth another look if that's
    # re-enabled.
    fig.canvas.draw()
    legend_bbox = ax.transAxes.inverted().transform_bbox(legend.get_window_extent(fig.canvas.get_renderer()))
    ax.text(
        0.70, legend_bbox.y0 - 0.03,
        "$m_{\\Phi} = 2000$ GeV\n$m_{Dark} = 20$ GeV\n$r_{inv} = 0.3$\n$\\lambda = 1$",
        transform=ax.transAxes,
        fontsize=17,
        ha='left',
        va='top',
        linespacing=1.6,
    )

    hep.cms.label("", data=False, ax=ax)

    plt.tight_layout()

    ratio_suffix = "" if with_ratio else "_no_ratio"
    tag_part = f"_{tag}" if tag else ""
    cat_part = "_allNSVJ" if combine_categories else ""
    tagger_suffix = tagger_cfg["suffix"]
    outfile = os.path.join(
        output_dir,
        f"preselection_MET_lund_weights{tagger_suffix}{cat_part}{tag_part}{ratio_suffix}.pdf"
    )

    fig.savefig(outfile, bbox_inches="tight")
    print(f"Saved {outfile}")

    plt.close(fig)


def _hist1(root_file, key):
    values, edges = root_file[key].to_numpy()
    return np.asarray(values, dtype=float), np.asarray(edges, dtype=float)


def _hist2(root_file, key):
    values, xedges, yedges = root_file[key].to_numpy()
    return np.asarray(values, dtype=float), np.asarray(xedges, dtype=float), np.asarray(yedges, dtype=float)


def _key_object(root_file, key_name):
    for key in root_file.file.root_directory._keys:
        if key.name() == key_name:
            return key.get()
    raise KeyError(key_name)


def _cms_label(ax, simulation=True, fontsize=22, lumi_fontsize=20):
    hep.cms.text("Simulation" if simulation else "", ax=ax, fontsize=fontsize)
    ax.text(1.0, 1.005, "(13 TeV)", transform=ax.transAxes, ha="right", va="bottom", fontsize=lumi_fontsize)


def _style_axis(ax):
    ax.tick_params(axis="both", which="major", labelsize=18, length=8, width=1.5, direction="in", top=True, right=True)
    ax.tick_params(axis="both", which="minor", length=4, width=1.2, direction="in", top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)


def _style_colorbar(cax, label):
    cax.tick_params(axis="y", which="major", labelsize=16, length=7, width=1.4, direction="in")
    cax.tick_params(axis="y", which="minor", length=3.5, width=1.1, direction="in")
    for spine in cax.spines.values():
        spine.set_linewidth(1.5)
    cax.set_ylabel(label, fontsize=18)


def _save(fig, output_dir, stem):
    os.makedirs(output_dir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(output_dir, f"{stem}.{ext}"), bbox_inches="tight")
    plt.close(fig)


def plot_tagged_svj_from_root(lund_root, output_dir, tagger, group_name, variations, normalize=True):
    prefix = f"tagged_svj/{tagger}"
    edges = _hist1(lund_root, f"{prefix}/Nominal")[1]
    nominal = _hist1(lund_root, f"{prefix}/Nominal")[0]
    lund_nom = _hist1(lund_root, f"{prefix}/LundNominal")[0]
    if normalize:
        nominal = nominal / nominal.sum() if nominal.sum() else nominal
        lund_nom = lund_nom / lund_nom.sum() if lund_nom.sum() else lund_nom

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.stairs(nominal, edges, color=COLOR_MAP["Nominal"], linewidth=2.5, label="Nominal")
    ax.stairs(lund_nom, edges, color=COLOR_MAP["Lund Nominal"], linewidth=2.5, label="Lund nominal")

    for var in variations:
        disp = DISPLAY_NAME[var]
        up = _hist1(lund_root, f"{prefix}/{var}Up")[0]
        down = _hist1(lund_root, f"{prefix}/{var}Down")[0]
        if normalize:
            up = up / up.sum() if up.sum() else up
            down = down / down.sum() if down.sum() else down
        lower = np.minimum(up, down)
        upper = np.maximum(up, down)
        ax.fill_between(
            edges,
            np.append(lower, lower[-1]),
            np.append(upper, upper[-1]),
            step="post",
            color=COLOR_MAP[f"{disp}Up"],
            alpha=0.28,
            linewidth=0,
            label=legend_label_for(f"{disp}Up"),
        )

    tagger_label = "WNAE" if tagger == "wnae" else "ParticleNet"
    ax.set_xlabel(rf"$n_{{\mathrm{{SVJ}}}}^{{\mathrm{{{tagger_label}}}}}$")
    ax.set_ylabel("Arbitrary units" if normalize else "Events")
    ax.set_xticks(np.arange(0.5, 4.5, 1.0), ["0", "1", "2", "3+"])
    ax.set_xlim(0, 4)
    tagged_ymax = max(ax.get_ylim()[1], np.max(lund_nom) * 1.35 if len(lund_nom) else 1.0)
    if normalize:
        tagged_ymax = max(tagged_ymax, 0.6)
    ax.set_ylim(0, tagged_ymax)
    ax.legend(fontsize=16, loc="upper right", frameon=False)
    _cms_label(ax)
    ax.text(
        0.04,
        0.84,
        r"$m_{\Phi}=2000$ GeV," "\n" r"$m_{\mathrm{dark}}=20$ GeV," "\n" r"$r_{\mathrm{inv}}=0.3,\ \lambda=1$",
        transform=ax.transAxes,
        fontsize=15,
        va="top",
    )
    fig.tight_layout()
    norm = "normalized" if normalize else "events"
    _save(fig, output_dir, f"lund_tagged_svj_{tagger}_{group_name}_{norm}")


def plot_preselection_met_from_root(lund_root, output_dir, group_name, variations):
    prefix = "preselection_met"
    nominal, edges = _hist1(lund_root, f"{prefix}/Nominal")
    lund_nom = _hist1(lund_root, f"{prefix}/LundNominal")[0]

    fig, (ax, rax) = plt.subplots(
        2, 1, figsize=(10, 9), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )
    ax.stairs(nominal, edges, color=COLOR_MAP["Nominal"], linewidth=2.5, label="Nominal")
    ax.stairs(lund_nom, edges, color=COLOR_MAP["Lund Nominal"], linewidth=2.5, label="Lund nominal")

    band_max = np.maximum(nominal, lund_nom)
    for var in variations:
        disp = DISPLAY_NAME[var]
        up = _hist1(lund_root, f"{prefix}/{var}Up")[0]
        down = _hist1(lund_root, f"{prefix}/{var}Down")[0]
        lower = np.minimum(up, down)
        upper = np.maximum(up, down)
        band_max = np.maximum(band_max, upper)
        ax.fill_between(
            edges,
            np.append(lower, lower[-1]),
            np.append(upper, upper[-1]),
            step="post",
            color=COLOR_MAP[f"{disp}Up"],
            alpha=0.28,
            linewidth=0,
            label=legend_label_for(f"{disp}Up"),
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_lower = np.where(lund_nom != 0, lower / lund_nom, 0.0)
            ratio_upper = np.where(lund_nom != 0, upper / lund_nom, 0.0)
        rax.fill_between(
            edges,
            np.append(ratio_lower, ratio_lower[-1]),
            np.append(ratio_upper, ratio_upper[-1]),
            step="post",
            color=COLOR_MAP[f"{disp}Up"],
            alpha=0.28,
            linewidth=0,
        )

    for values, color, width in ((nominal, COLOR_MAP["Nominal"], 2.5), (lund_nom, COLOR_MAP["Lund Nominal"], 2.5)):
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(lund_nom != 0, values / lund_nom, 0.0)
        rax.stairs(ratio, edges, color=color, linewidth=width)

    ax.set_yscale("log")
    ax.set_ylim(0.1, max(1.0, np.max(band_max)) * 80.0)
    ax.set_ylabel("Events")
    rax.axhline(1.0, color="black", linestyle="--", linewidth=1)
    rax.set_ylim(0.5, 1.5)
    rax.set_ylabel("Ratio to\nLund nominal", fontsize=14)
    rax.set_xlabel(r"$p_{\mathrm{T}}^{\mathrm{miss}}$ [GeV]")
    ax.legend(fontsize=15, ncol=2, loc="upper right", frameon=False)
    _cms_label(ax)
    ax.text(0.04, 0.82, "Preselection", transform=ax.transAxes, fontsize=18, va="top")
    fig.tight_layout()
    _save(fig, output_dir, f"lund_preselection_met_{group_name}")


def plot_lund_plane_from_root(distortion_root, output_dir):
    for idx in range(1, 7):
        pt_label = LUND_CORRECTION_PT_LABELS.get(idx, rf"correction $p_{{\mathrm{{T}}}}$ interval {idx}")
        ratio, xedges, yedges = _hist2(distortion_root, f"lundPlane_bin{idx}_ratio")
        unc_obj = _key_object(distortion_root, f"lundplots/distortion/lundPlane_bin{idx}_ratio_unc")
        unc, _, _ = unc_obj.to_numpy()

        fig, (ax, cax) = plt.subplots(
            1, 2, figsize=(8.6, 7.2),
            gridspec_kw={"width_ratios": [1, 0.055], "wspace": 0.08},
        )
        mesh = ax.pcolormesh(xedges, yedges, ratio.T, cmap="viridis", shading="auto")
        cb = fig.colorbar(mesh, cax=cax)
        _style_colorbar(cax, "Lund-plane signal / W+jets ratio")
        ax.set_xlabel(r"$\ln(1/\Delta R)$", fontsize=22)
        ax.set_ylabel(r"$\ln(k_{\mathrm{T}}/\mathrm{GeV})$", fontsize=22)
        _style_axis(ax)
        _cms_label(ax, fontsize=19, lumi_fontsize=18)
        ax.text(
            0.03,
            0.95,
            pt_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=16,
            color="white",
            fontweight="bold",
        )
        fig.tight_layout()
        _save(fig, output_dir, f"lund_plane_ratio_bin{idx}")

        fig, (ax, cax) = plt.subplots(
            1, 2, figsize=(8.6, 7.2),
            gridspec_kw={"width_ratios": [1, 0.055], "wspace": 0.08},
        )
        mesh = ax.pcolormesh(xedges, yedges, unc.T, cmap="viridis", shading="auto")
        cb = fig.colorbar(mesh, cax=cax)
        _style_colorbar(cax, "Lund-plane distortion relative uncertainty")
        ax.set_xlabel(r"$\ln(1/\Delta R)$", fontsize=22)
        ax.set_ylabel(r"$\ln(k_{\mathrm{T}}/\mathrm{GeV})$", fontsize=22)
        _style_axis(ax)
        _cms_label(ax, fontsize=19, lumi_fontsize=18)
        ax.text(
            0.03,
            0.95,
            pt_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=16,
            color="white",
            fontweight="bold",
        )
        fig.tight_layout()
        _save(fig, output_dir, f"lund_plane_distortion_unc_bin{idx}")


def plot_from_root_inputs(lund_root_path, distortion_root_path, output_dir):
    with uproot.open(lund_root_path) as lund_root:
        for tagger in ("wnae", "pnet"):
            for group_name, variations in TAGGED_SVJ_GROUPS.items():
                plot_tagged_svj_from_root(lund_root, output_dir, tagger, group_name, variations, normalize=True)
        plot_preselection_met_from_root(lund_root, output_dir, "SysDistortion", TAGGED_SVJ_GROUPS["SysDistortion"])
        plot_preselection_met_from_root(lund_root, output_dir, "PtStat", TAGGED_SVJ_GROUPS["PtStat"])
    with uproot.open(distortion_root_path) as distortion_root:
        plot_lund_plane_from_root(distortion_root, output_dir)
    print(f"Saved Lund plots to {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--from-root", action="store_true", help="draw plots from existing Lund ROOT caches")
    parser.add_argument("--lund-root", default="lund-histograms.root")
    parser.add_argument("--distortion-root", default="distortion-LJP.root")
    parser.add_argument("--output-dir", default="assets_v2/lund_systematics")
    args, _ = parser.parse_known_args()
    if args.from_root:
        plot_from_root_inputs(args.lund_root, args.distortion_root, args.output_dir)
        raise SystemExit(0)

    systematics_output_dir = "./signal_systematic_plots"

    PRESELECTION_CFG = {
        "cat_superscript": "",
        "corner_label": "Preselection",
        "suffix": "",
    }

    with open_output_root() as root_file:
        # --- tagged-SVJ (N_SVJ) histograms, both taggers ---
        # Loaded once per tagger and reused for both plot variants — normalize
        # only rescales the histogram at plot time, it doesn't need its own
        # load. Only the paired-systematic plots are generated (no single
        # "all four at once" plot); the full ROOT histograms are still
        # written directly from the loaded data regardless.
        for tagger_cfg in TAGGER_CONFIGS:
            n_tagged_jets, pu_weights, weights_nominal_lund, variation_weights = load_tagged_svj_data(tagger_cfg)
            write_tagged_svj_histograms(root_file, tagger_cfg, n_tagged_jets, pu_weights, weights_nominal_lund,
                                         variation_weights)

            for tag, group_vars in TAGGED_SVJ_GROUPS.items():
                make_tagged_svj_plot(tagger_cfg, n_tagged_jets, pu_weights, weights_nominal_lund, variation_weights,
                                      normalize=True, variations=group_vars, tag=tag)
                make_tagged_svj_plot(tagger_cfg, n_tagged_jets, pu_weights, weights_nominal_lund, variation_weights,
                                      normalize=False, variations=group_vars, tag=tag)

        # --- raw Lund weight distributions (combined + one plot per variation) ---
        if MAKE_LUND_WEIGHT_PLOTS:
            lund_weights, variation_raw, lund_weights_range = load_lund_weights_data()
            make_lund_weights_plot(lund_weights, variation_raw, lund_weights_range, root_file=root_file)
            for var, color in zip(VARIATIONS, VARIATION_COLORS):
                make_lund_weight_variation_plot(var, color, lund_weights, variation_raw, lund_weights_range)

        # --- preselection MET spectrum, tagger-agnostic (not split by n_SVJ) ---
        # Doesn't mask on n_tagged at all, so it's generated unconditionally
        # (independent of MAKE_NSVJ_PLOTS) from a single canonical skim
        # rather than once per tagger, and kept at its own top-level ROOT
        # prefix ("preselection_met") separate from the n_SVJ-split
        # histograms below ("nsvj_met").
        met_all, _, weights_all = load_events(TAGGER_CONFIGS[0]["subdir_template"])
        write_root_histograms(root_file, "preselection_met", met_all, weights_all)

        for with_ratio in (True, False):
            plot_systematic_comparison(met_all, None, weights_all, systematics_output_dir, PRESELECTION_CFG,
                                        with_ratio=with_ratio, combine_categories=True, band_systematics=True)
            for tag, group_names in SYSTEMATIC_PAIR_GROUPS.items():
                names = ["Nominal", "Lund Nominal"] + group_names
                plot_systematic_comparison(met_all, None, weights_all, systematics_output_dir, PRESELECTION_CFG,
                                            with_ratio=with_ratio, names=names, tag=tag, combine_categories=True,
                                            band_systematics=True)

        # --- preselection MET spectra, split by n_SVJ category ---
        if MAKE_NSVJ_PLOTS:
            for tagger_cfg in TAGGER_CONFIGS:
                met, n_tagged, weights = load_events(tagger_cfg["subdir_template"], tagger_cfg)

                for cat in CATEGORIES:
                    mask = category_mask(n_tagged, cat)
                    write_root_histograms(root_file, f"nsvj_met/{tagger_cfg['name']}/{cat}", met, weights, mask)

                for with_ratio in (True, False):
                    plot_systematic_comparison(met, n_tagged, weights, systematics_output_dir, tagger_cfg,
                                                with_ratio=with_ratio, combine_categories=False, band_systematics=True)
                    for tag, group_names in SYSTEMATIC_PAIR_GROUPS.items():
                        names = ["Nominal", "Lund Nominal"] + group_names
                        plot_systematic_comparison(met, n_tagged, weights, systematics_output_dir, tagger_cfg,
                                                    with_ratio=with_ratio, names=names, tag=tag,
                                                    combine_categories=False, band_systematics=True)
