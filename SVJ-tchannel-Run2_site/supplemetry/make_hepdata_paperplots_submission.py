#!/usr/bin/env python3
"""Build the HEPData submission for the figures in ``PaperPLots``.

The numerical values are read from the same ROOT caches, AUC arrays, limit
CSVs, and FitDiagnostics files used to make the paper figures.  No values are
digitized from the PDFs.

Usage::

    source condor/initCondor.sh
    python3 SVJ-tchannel-Run2_site/supplemetry/make_hepdata_paperplots_submission.py

The output is written next to this script in ``hepdata_paperplots_submission``.
"""

import contextlib
import csv
import importlib.util
import io
import math
import os
import sys
from types import SimpleNamespace

import numpy as np
from hepdata_lib import RootFileReader, Submission, Table, Uncertainty, Variable


HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.dirname(HERE)
BASE = os.path.dirname(SITE)
PAPER_DIR = os.path.join(SITE, "PaperPLots")
OUTPUT_DIR = os.path.join(HERE, "hepdata_paperplots_submission")

SUPPLEMENTARY_CACHE = os.path.join(HERE, "supplemetry_cache.root")
FIGURE2_CACHE = os.path.join(BASE, "Figure2_plots", "WNAE_alljets_dataMC_vars_Run2_wlundcorrection_v4.root")
DNN_CACHE = os.path.join(BASE, "EventDNNScore_plots", "dnn_scores_wlundcorrection.root")
LIMITS_DIR = os.path.join(BASE, "t-channel_plotting_scripts", "limits")
POSTFIT_DIR = os.path.join(BASE, "t-channel_plotting_scripts", "postfits")

PROC_LABEL = {
    "QCD": "QCD multijet",
    "TTJets": r"$t\bar{t}$+jets",
    "WJetsToLNu": r"$W\rightarrow\ell\nu$+jets",
    "ZJetsToNuNu": r"$Z\rightarrow\nu\nu$+jets",
    "ST": "Single top",
}
PROC_ORDER = ["QCD", "TTJets", "WJetsToLNu", "ZJetsToNuNu", "ST"]

SIGNALS = [
    (
        "t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
        r"$m_{\Phi}=600$ GeV, $r_{\mathrm{inv}}=0.3$",
        10.0,
    ),
    (
        "t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1",
        r"$m_{\Phi}=2000$ GeV, $r_{\mathrm{inv}}=0.1$",
        500.0,
    ),
    (
        "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
        r"$m_{\Phi}=2000$ GeV, $r_{\mathrm{inv}}=0.3$",
        500.0,
    ),
]

THEORY_XSEC_PB = {
    500: 41.76,
    600: 19.65,
    700: 10.38,
    800: 6.006,
    900: 3.717,
    1000: 2.413,
    1500: 0.4832,
    2000: 0.1634,
    2500: 0.07178,
    3000: 0.03603,
    3500: 0.01999,
    4000: 0.01185,
}


def require_files(paths):
    missing = [path for path in paths if not os.path.isfile(path)]
    if missing:
        raise SystemExit("Missing required input(s):\n  " + "\n  ".join(missing))


def paper_table(number, description, image, observables):
    table = Table(f"Figure {number}")
    table.description = description
    table.location = f"Data from Figure {number}."
    table.keywords["observables"] = observables
    table.keywords["reactions"] = ["P P --> JET JET"]
    table.keywords["cmenergies"] = [13000.0]
    table.add_image(os.path.join(PAPER_DIR, image))
    return table


def add_stat(variable, errors, label="stat"):
    uncertainty = Uncertainty(label, is_symmetric=True)
    uncertainty.values = [float(value) for value in errors]
    variable.add_uncertainty(uncertainty)


def hist_copy(hist):
    return {
        "x_edges": [(float(low), float(high)) for low, high in hist["x_edges"]],
        "y": np.asarray(hist["y"], dtype=float),
        "dy": np.asarray(hist["dy"], dtype=float),
    }


def wnae_overflow_display(hist):
    """Apply the exact display binning used by Figure2_makerusingskims.py."""
    source = hist_copy(hist)
    xmin, xmax, width, working_point = 5.0, 60.0, 5.0, 18.284
    start = working_point - math.ceil((working_point - xmin) / width) * width
    stop = start + math.ceil((xmax - start) / width) * width
    edges = np.arange(start, stop + 0.5 * width, width)
    edges[np.argmin(np.abs(edges - working_point))] = working_point
    edges = np.asarray([round(float(edge), 6) for edge in edges])
    values = np.zeros(len(edges) - 1)
    variances = np.zeros(len(edges) - 1)
    for (low, high), value, error in zip(source["x_edges"], source["y"], source["dy"]):
        center = 0.5 * (low + high)
        if center < xmin:
            index = 0
        elif center >= xmax:
            index = len(values) - 1
        else:
            index = int(np.searchsorted(edges, center, side="right") - 1)
            index = min(max(index, 0), len(values) - 1)
        values[index] += value
        variances[index] += error * error
    return {
        "x_edges": list(zip(edges[:-1], edges[1:])),
        "y": values,
        "dy": np.sqrt(variances),
    }


def add_distribution_table(
    submission,
    number,
    description,
    image,
    cache,
    background_path,
    signal_path,
    x_title,
    x_units="",
    data_path=None,
    signal_scales=False,
    transform=None,
    y_units="Events",
):
    reader = RootFileReader(cache)
    table = paper_table(number, description, image, ["N"])
    transform = transform or hist_copy

    backgrounds = []
    for process in PROC_ORDER:
        backgrounds.append(transform(reader.read_hist_1d(background_path.format(proc=process))))
    x = Variable(x_title, is_independent=True, is_binned=True, units=x_units)
    x.values = backgrounds[0]["x_edges"]
    table.add_variable(x)

    total = np.zeros_like(backgrounds[0]["y"])
    total_variance = np.zeros_like(total)
    for process, hist in zip(PROC_ORDER, backgrounds):
        variable = Variable(
            PROC_LABEL[process], is_independent=False, is_binned=False, units=y_units,
            zero_uncertainties_warning=False,
        )
        variable.values = hist["y"].tolist()
        add_stat(variable, hist["dy"])
        table.add_variable(variable)
        total += hist["y"]
        total_variance += hist["dy"] ** 2

    total_variable = Variable(
        "Total simulation", is_independent=False, is_binned=False, units=y_units,
        zero_uncertainties_warning=False,
    )
    total_variable.values = total.tolist()
    add_stat(total_variable, np.sqrt(total_variance), "MC stat")
    table.add_variable(total_variable)

    if data_path:
        hist = transform(reader.read_hist_1d(data_path))
        variable = Variable(
            "Data", is_independent=False, is_binned=False, units=y_units,
            zero_uncertainties_warning=False,
        )
        variable.values = hist["y"].tolist()
        add_stat(variable, hist["dy"])
        table.add_variable(variable)

    for signal_name, signal_label, display_scale in SIGNALS:
        hist = transform(reader.read_hist_1d(signal_path.format(signal=signal_name)))
        scale = display_scale if signal_scales else 1.0
        label = signal_label + (rf" ($\times {display_scale:g}$)" if signal_scales else "")
        variable = Variable(
            label, is_independent=False, is_binned=False, units=y_units,
            zero_uncertainties_warning=False,
        )
        variable.values = (hist["y"] * scale).tolist()
        add_stat(variable, hist["dy"] * scale)
        table.add_variable(variable)

    submission.add_table(table)


def add_nminus1_tables(submission, start=1):
    specs = [
        (
            "h_dPhiMinjMETAK8_pre__dPhiMin",
            r"$\min\Delta\phi(J,p_{T}^{\mathrm{miss}})$",
            "",
            "Nminus1/nminus1_dPhiMinjMETAK8.pdf",
            "Normalized minimum azimuthal separation between an AK8 jet and missing transverse momentum. "
            "The requirement on the plotted variable is omitted while the other preselection requirements are applied.",
        ),
        (
            "h_ST_pre__stcut",
            r"$S_{T}$",
            "GeV",
            "Nminus1/nminus1_ST.pdf",
            "Normalized scalar sum of AK8 jet transverse momenta. The requirement on the plotted variable is "
            "omitted while the other preselection requirements are applied.",
        ),
        (
            "h_MET_pre__metcut",
            r"$p_{T}^{\mathrm{miss}}$",
            "GeV",
            "Nminus1/nminus1_MET.pdf",
            "Normalized missing transverse momentum distribution. The requirement on the plotted variable is "
            "omitted while the other preselection requirements are applied.",
        ),
    ]
    for offset, (key, title, units, image, description) in enumerate(specs):
        add_distribution_table(
            submission,
            start + offset,
            description,
            image,
            SUPPLEMENTARY_CACHE,
            f"nminus1/shapes/groups/{{proc}}/{key}",
            f"nminus1/shapes/signals/{{signal}}/{key}",
            title,
            units,
            y_units="Arbitrary units",
        )


def add_score_tables(submission, start=4):
    specs = [
        (
            DNN_CACHE,
            "DNNScore",
            "Event classifier score",
            "Scores/event_DNN_score_Run2.pdf",
            "Distribution of the event-level classifier score after preselection for the Run 2 data, simulated "
            "backgrounds, and benchmark signals.",
            None,
        ),
        (
            FIGURE2_CACHE,
            "allWNAEPt200To300LossAK8",
            r"WNAE ($200 < p_{T} < 300$ GeV) score ($J$)",
            "Scores/WNAE_score_200to300_Run2.pdf",
            "Distribution of the WNAE score for AK8 jets with transverse momentum between 200 and 300 GeV "
            "after event preselection.",
            wnae_overflow_display,
        ),
        (
            FIGURE2_CACHE,
            "allPNetScoreAK8",
            r"ParticleNet score ($J$)",
            "Scores/ParticleNet_score_Run2.pdf",
            "Distribution of the ParticleNet score for all selected AK8 jets after event preselection.",
            None,
        ),
    ]
    for offset, (cache, key, title, image, description, transform) in enumerate(specs):
        add_distribution_table(
            submission,
            start + offset,
            description,
            image,
            cache,
            f"raw/groups/Run2/{{proc}}/{key}",
            f"raw/signals/Run2/{{signal}}/{key}",
            title,
            data_path=f"raw/data/Run2/observed/{key}",
            signal_scales=True,
            transform=transform,
        )


def load_auc_module():
    path = os.path.join(BASE, "makeAUCtables_pnet.py")
    spec = importlib.util.spec_from_file_location("paper_auc_values", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def add_auc_table(submission, number, image, description, z_title, mmed, rinv, values):
    values = np.asarray(values, dtype=float)
    if values.shape != (len(mmed), len(rinv)):
        raise ValueError(f"AUC grid has shape {values.shape}; expected {(len(mmed), len(rinv))}")
    table = paper_table(number, description, image, ["AUC"])
    x = Variable(r"$m_{\Phi}$", is_independent=True, is_binned=False, units="GeV")
    y = Variable(r"$r_{\mathrm{inv}}$", is_independent=True, is_binned=False)
    z = Variable(z_title, is_independent=False, is_binned=False)
    x.values = [float(mx) for mx in mmed for _ in rinv]
    y.values = [float(ry) for _ in mmed for ry in rinv]
    z.values = values.reshape(-1).tolist()
    table.add_variable(x)
    table.add_variable(y)
    table.add_variable(z)
    submission.add_table(table)


def add_auc_tables(submission, start=7):
    auc = load_auc_module()
    add_auc_table(
        submission,
        start,
        "AUC_QCD/event_DNN_AUC_vs_QCD_mMed_rinv.pdf",
        "Area under the receiver operating characteristic curve for the event-level classifier against QCD multijet background.",
        "Event classifier AUC vs. QCD",
        auc.DNN_BACKGROUND_MMED,
        auc.DNN_BACKGROUND_RINV,
        auc.DNN_BACKGROUND_GRIDS["qcd_from_dnn"]["raw"],
    )
    add_auc_table(
        submission,
        start + 1,
        "AUC_QCD/WNAE_AUC_vs_QCD_mMed_rinv.pdf",
        "Area under the receiver operating characteristic curve for the WNAE tagger against QCD multijet background.",
        "WNAE AUC vs. QCD",
        auc.MMED_WNAE,
        auc.RINV_WNAE,
        auc.WNAE_RAW,
    )
    add_auc_table(
        submission,
        start + 2,
        "AUC_QCD/ParticleNet_AUC_vs_QCD_mMed_rinv.pdf",
        "Area under the receiver operating characteristic curve for the ParticleNet tagger against QCD multijet background.",
        "ParticleNet AUC vs. QCD",
        auc.MMED_PNET,
        auc.RINV_PNET,
        auc.PNET_RAW,
    )


def load_postfit_module():
    path = os.path.join(POSTFIT_DIR, "postfit_plotter.py")
    spec = importlib.util.spec_from_file_location("paper_postfit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.ROOT.gErrorIgnoreLevel = module.ROOT.kError
    return module


def root_hist_arrays(hist):
    edges = [(float(hist.GetBinLowEdge(i)), float(hist.GetBinLowEdge(i + 1))) for i in range(1, hist.GetNbinsX() + 1)]
    values = [float(hist.GetBinContent(i)) for i in range(1, hist.GetNbinsX() + 1)]
    errors = [float(hist.GetBinError(i)) for i in range(1, hist.GetNbinsX() + 1)]
    return edges, values, errors


def postfit_payload(tagger):
    postfit = load_postfit_module()
    tag_dir = "PNET" if tagger == "ParticleNet" else "WNAE"
    root_path = os.path.join(
        POSTFIT_DIR,
        tag_dir,
        "fitDiagnosticscombinedFitDiagnostics_mMed2000_mDark20_rinv0p3_yukawa1_Run2.root",
    )
    card_path = os.path.join(
        POSTFIT_DIR,
        tag_dir,
        "combinedFitDiagnostics_mMed2000_mDark20_rinv0p3_yukawa1_Run2.txt",
    )
    mapping = postfit._cli_mapping(card_path)
    abcd = SimpleNamespace(remap_back_postfit=False, map_regions={r: r for r in "ABCD"})
    by_year = {}
    # The plotting helpers are intentionally reused so RooFit density-bin and
    # overflow corrections exactly match the published postfit figures.
    with contextlib.redirect_stdout(io.StringIO()):
        for category, channel_map in mapping.items():
            hists = {
                "fit_b_bkg": postfit.__get_bkg_histograms(root_path, "fit_b", mapping, category, abcd),
                "fit_b_data": postfit.__get_data_histograms(root_path, None, "fit_b", mapping, category, abcd),
                "prefit_signal": postfit.__get_signal_histograms(root_path, "prefit", mapping, category, abcd),
            }
            base, year = (category.split("Y", 1) + ["Run2"])[:2] if "Y" in category else (category, "Run2")
            by_year.setdefault(base, {})[year] = hists

        combined = {}
        for category, years in by_year.items():
            combined[category] = {}
            for fit, flavor in (("fit_b", "bkg"), ("fit_b", "data"), ("prefit", "signal")):
                combined[category][f"{fit}_{flavor}"] = postfit.__combine_histograms_years(years, fit, flavor)
        order = [category for category in ("0SVJ", "1SVJ", "2SVJ", "3PSVJ") if category in combined]
        return postfit.__collect_stitched_region_payload(
            combined,
            "A",
            "fit_b_bkg",
            "fit_b_data",
            {"Signal": "prefit_signal"},
            order,
            2 if tagger == "ParticleNet" else 4,
        )


def add_postfit_table(submission, number, tagger):
    payload = postfit_payload(tagger)
    image = (
        "Postfit_ARegion_BOnly/ParticleNet_postfit_bonly_A_Run2.pdf"
        if tagger == "ParticleNet"
        else "Postfit_ARegion_BOnly/WNAE_postfit_bonly_A_Run2.pdf"
    )
    table = paper_table(
        number,
        f"Background-only postfit missing transverse momentum distributions in region A, split by {tagger} "
        "semivisible-jet multiplicity. The signal curve is shown at its nominal prefit normalization.",
        image,
        ["N"],
    )
    categories, pt_bins, backgrounds, background_errors, data, data_errors, signal, signal_errors = ([] for _ in range(8))
    for category, h_bkg, h_data, h_signal in zip(
        payload["categories"], payload["bkg_hists"], payload["data_hists"], payload["sig_hists"]["Signal"]
    ):
        edges, b_values, b_errors = root_hist_arrays(h_bkg)
        _, d_values, d_errors = root_hist_arrays(h_data)
        _, s_values, s_errors = root_hist_arrays(h_signal)
        categories.extend([category.replace("3PSVJ", "3+SVJ")] * len(edges))
        pt_bins.extend(edges)
        backgrounds.extend(b_values)
        background_errors.extend(b_errors)
        data.extend(d_values)
        data_errors.extend(d_errors)
        signal.extend(s_values)
        signal_errors.extend(s_errors)

    category_var = Variable(rf"$n_{{\mathrm{{SVJ}}}}^{{{'PN' if tagger == 'ParticleNet' else 'WNAE'}}}$ category", True, False)
    category_var.values = categories
    pt_var = Variable(r"$p_{T}^{\mathrm{miss}}$", True, True, "GeV")
    pt_var.values = pt_bins
    table.add_variable(category_var)
    table.add_variable(pt_var)
    for label, values, errors, uncertainty_label in (
        ("Postfit background", backgrounds, background_errors, "postfit"),
        ("Data", data, data_errors, "stat"),
        (r"Signal ($m_{\Phi}=2000$ GeV, $m_{\mathrm{dark}}=20$ GeV, $r_{\mathrm{inv}}=0.3$, $\lambda=1$)", signal, signal_errors, "MC stat"),
    ):
        variable = Variable(label, False, False, "Events", zero_uncertainties_warning=False)
        variable.values = values
        add_stat(variable, errors, uncertainty_label)
        table.add_variable(variable)
    submission.add_table(table)


def token_float(value):
    return float(str(value).replace("p", "."))


def read_limit_points(path, axes, fixed, min_y=None):
    points = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            numeric = {key: token_float(row[key]) for key in row}
            if any(not np.isclose(numeric[key], value) for key, value in fixed.items()):
                continue
            if min_y is not None and numeric[axes[1]] < min_y:
                continue
            points.append(numeric)
    points.sort(key=lambda row: (row[axes[0]], row[axes[1]]))
    return points


def add_limit_table(submission, number, tagger, axes, fixed):
    tag_lower = tagger.lower()
    axis_slug = "_".join(axes)
    csv_path = os.path.join(LIMITS_DIR, f"unblinded_gapVeto_perCateYearNonclosure_{tagger.upper()}-2.csv")
    points = read_limit_points(csv_path, axes, fixed, min_y=0.5 if axes[1] == "yukawa" else None)
    image = f"Limits2D/limits2d_{tag_lower}_{axis_slug}.pdf"
    table = paper_table(
        number,
        f"Observed and expected 95% confidence level upper limits on the signal production cross section for the "
        f"{tagger} analysis in the {axes[0]}-{axes[1]} plane.",
        image,
        ["SIG"],
    )
    labels = {
        "mMed": (r"$m_{\Phi}$", "GeV"),
        "mDark": (r"$m_{\mathrm{dark}}$", "GeV"),
        "rinv": (r"$r_{\mathrm{inv}}$", ""),
        "yukawa": (r"$\lambda$", ""),
    }
    for axis in axes:
        name, units = labels[axis]
        variable = Variable(name, True, False, units)
        variable.values = [row[axis] for row in points]
        table.add_variable(variable)

    columns = [
        ("Observed 95% CL upper limit", "obs_lim"),
        ("Median expected 95% CL upper limit", "expected"),
        ("Expected 95% CL upper limit, -2 s.d.", "expected_m2sigma"),
        ("Expected 95% CL upper limit, -1 s.d.", "expected_m1sigma"),
        ("Expected 95% CL upper limit, +1 s.d.", "expected_p1sigma"),
        ("Expected 95% CL upper limit, +2 s.d.", "expected_p2sigma"),
    ]
    for label, column in columns:
        variable = Variable(label, False, False, "pb")
        variable.values = [row[column] * THEORY_XSEC_PB[int(round(row["mMed"]))] for row in points]
        table.add_variable(variable)
    theory = Variable("Theory cross section", False, False, "pb")
    theory.values = [THEORY_XSEC_PB[int(round(row["mMed"]))] for row in points]
    table.add_variable(theory)
    submission.add_table(table)


def add_limit_tables(submission, start=12):
    specs = [
        ("PNET", ("mMed", "mDark"), {"rinv": 0.3, "yukawa": 1.0}),
        ("PNET", ("mMed", "rinv"), {"mDark": 20.0, "yukawa": 1.0}),
        ("PNET", ("mMed", "yukawa"), {"mDark": 20.0, "rinv": 0.3}),
        ("WNAE", ("mMed", "mDark"), {"rinv": 0.3, "yukawa": 1.0}),
        ("WNAE", ("mMed", "rinv"), {"mDark": 20.0, "yukawa": 1.0}),
        ("WNAE", ("mMed", "yukawa"), {"mDark": 20.0, "rinv": 0.3}),
    ]
    for offset, (tagger, axes, fixed) in enumerate(specs):
        add_limit_table(submission, start + offset, tagger, axes, fixed)


def main():
    require_files(
        [
            SUPPLEMENTARY_CACHE,
            FIGURE2_CACHE,
            DNN_CACHE,
            os.path.join(BASE, "makeAUCtables_pnet.py"),
            os.path.join(LIMITS_DIR, "unblinded_gapVeto_perCateYearNonclosure_PNET-2.csv"),
            os.path.join(LIMITS_DIR, "unblinded_gapVeto_perCateYearNonclosure_WNAE-2.csv"),
        ]
    )
    submission = Submission()
    submission.comment = (
        "Numerical data for the t-channel semivisible-jet paper figures. Signal benchmark distributions include "
        "the display scale factors stated in their variable names."
    )

    print("[1/5] N-1 distributions")
    add_nminus1_tables(submission)
    print("[2/5] classifier-score distributions")
    add_score_tables(submission)
    print("[3/5] AUC maps")
    add_auc_tables(submission)
    print("[4/5] postfit region-A distributions")
    add_postfit_table(submission, 10, "ParticleNet")
    add_postfit_table(submission, 11, "WNAE")
    print("[5/5] two-dimensional limits")
    add_limit_tables(submission)

    submission.create_files(OUTPUT_DIR, remove_old=True)
    print(f"[DONE] wrote 17 paper-figure tables to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
