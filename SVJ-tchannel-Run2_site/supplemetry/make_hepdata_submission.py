#!/usr/bin/env python3
"""
Build a HEPData submission (submission.yaml + one data YAML + one thumbnail
per table) for the N-1, ParticleNet-input, WNAE-input, and DNN-input plots in
supplementryv2.pdf, reading directly from the same merged cache
(supplemetry_cache.root) that make_supplementary_plots.py draws from.

Uses hepdata_lib (https://github.com/HEPData/hepdata_lib), the library
HEPData itself recommends, rather than hand-writing the YAML schema.

Usage:
  source condor/initCondor.sh
  python3 make_supplementary_plots.py     # if supplemetry_cache.root doesn't exist yet
  python3 make_hepdata_submission.py
"""
import os

from hepdata_lib import Submission, Table, Variable, Uncertainty, RootFileReader

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "supplemetry_cache.root")
OUTPUT_DIR = os.path.join(HERE, "hepdata_submission")

PROC_LABEL = {
    "QCD": "QCD multijet",
    "TTJets": r"$t\bar{t}$+jets",
    "WJetsToLNu": r"$W\rightarrow\ell\nu$+jets",
    "ZJetsToNuNu": r"$Z\rightarrow\nu\nu$+jets",
    "ST": "Single top",
}
PROC_ORDER = ["QCD", "TTJets", "WJetsToLNu", "ZJetsToNuNu", "ST"]

# Matches Figure2_makerusingskims.DEFAULT_SIGNAL_DIRS -- the three reference
# signal points overlaid throughout this document.
SIGNALS = [
    ("t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1", r"$m_{\Phi}=600$ GeV, $r_{\mathrm{inv}}=0.3$"),
    ("t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1", r"$m_{\Phi}=2000$ GeV, $r_{\mathrm{inv}}=0.1$"),
    ("t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", r"$m_{\Phi}=2000$ GeV, $r_{\mathrm{inv}}=0.3$"),
]

# One entry per plot in supplementryv2.pdf, in the same order the tex draws
# them (N-1, then ParticleNet inputs, then WNAE inputs, then DNN inputs) so
# TABLE_SPECS[i]'s position gives the right "Additional Figure N" number.
#
# fields:
#   key           cache variable name
#   section       top-level namespace in supplemetry_cache.root
#   image         PDF under assets_v2/<section>/ used for the thumbnail
#   x_title/units header for the independent (x-axis) variable
#   description   exact caption text from supplementryv2.tex
#   crop          optional (xmin, xmax) to match the plot's display crop --
#                 without this, variables like soft-drop mass would publish
#                 hundreds of empty high-mass bins that were never shown
#   has_data      whether to include an observed-data dependent variable
TABLE_SPECS = [
    dict(key="h_j1PtAK8_pre__stcut", section="nminus1", image="nminus1/h_j1PtAK8_pre__stcut_normalized.pdf",
         x_title=r"$p_{T}(J_{1})$", x_units="GeV", has_data=False,
         description="The normalized distribution of leading AK8 jet $p_{\\mathrm{T}}$ for simulated background "
                      "and signal events. The requirement on the plotted variable is omitted, while all other "
                      "preselection requirements are applied. The vertical dotted line indicates the preselection "
                      "(final selection) requirement on the variable."),
    dict(key="h_j2PtAK8_pre__stcut", section="nminus1", image="nminus1/h_j2PtAK8_pre__stcut_normalized.pdf",
         x_title=r"$p_{T}(J_{2})$", x_units="GeV", has_data=False,
         description="The normalized distribution of subleading AK8 jet $p_{\\mathrm{T}}$ for simulated background "
                      "and signal events. The requirement on the plotted variable is omitted, while all other "
                      "preselection requirements are applied. The vertical dotted line indicates the preselection "
                      "(final selection) requirement on the variable."),
    dict(key="h_jPtAK8_pre__stcut", section="nminus1", image="nminus1/h_jPtAK8_pre__stcut_normalized.pdf",
         x_title=r"$p_{T}(J)$", x_units="GeV", has_data=False,
         description="The normalized distribution of AK8 jet $p_{\\mathrm{T}}$ for simulated background and "
                      "signal events. The requirement on the plotted variable is omitted, while all other "
                      "preselection requirements are applied. The vertical dotted line indicates the preselection "
                      "(final selection) requirement on the variable."),
    dict(key="h_njetsAK8_pre__2jetsAK8", section="nminus1", image="nminus1/h_njetsAK8_pre__2jetsAK8_normalized.pdf",
         x_title="Number of AK8 jets", x_units="", has_data=False,
         description="The normalized distribution of AK8 jet multiplicity for simulated background and signal "
                      "events. The requirement on the plotted variable is omitted, while all other preselection "
                      "requirements are applied. The vertical dotted line indicates the preselection (final "
                      "selection) requirement on the variable."),
    dict(key="del_r", section="pnet", image="pnet_inputs/del_r_linear_raw_ratio.pdf",
         x_title=r"$\Delta R$", x_units="", has_data=True,
         description="The distribution of the ParticleNet input variable $\\Delta R$ of jet constituents for Run 2 "
                      "events, comparing signal and background samples after the Lund-plane correction."),
    dict(key="log_e_jete", section="pnet", image="pnet_inputs/log_e_jete_linear_raw_ratio.pdf",
         x_title=r"$\log(E/E_{J})$", x_units="", has_data=True, crop=(-8.0, -1.0),
         description="The distribution of the ParticleNet input variable constituent energy fraction for Run 2 "
                      "events, comparing signal and background samples after the Lund-plane correction."),
    dict(key="log_e", section="pnet", image="pnet_inputs/log_e_linear_raw_ratio.pdf",
         x_title=r"$\log(E)$", x_units="", has_data=True, crop=(-1.0, 6.0),
         description="The distribution of the ParticleNet input variable constituent energy for Run 2 events, "
                      "comparing signal and background samples after the Lund-plane correction."),
    dict(key="log_pt_jetpt", section="pnet", image="pnet_inputs/log_pt_jetpt_linear_raw_ratio.pdf",
         x_title=r"$\log(p_{T}/p_{T}^{J})$", x_units="", has_data=True, crop=(-8.0, -1.0),
         description="The distribution of the ParticleNet input variable constituent $p_{\\mathrm{T}}$ fraction "
                      "for Run 2 events, comparing signal and background samples after the Lund-plane correction."),
    dict(key="log_pt", section="pnet", image="pnet_inputs/log_pt_linear_raw_ratio.pdf",
         x_title=r"$\log(p_{T})$", x_units="", has_data=True, crop=(-1.0, 5.0),
         description="The distribution of the ParticleNet input variable constituent $p_{\\mathrm{T}}$ for Run 2 "
                      "events, comparing signal and background samples after the Lund-plane correction."),
    dict(key="allPtDAK8", section="wnae", image="wnae_inputs/allPtDAK8_log_raw_ratio.pdf",
         x_title=r"$D_{p_{T}}(J)$", x_units="", has_data=True,
         description="The distribution of the WNAE tagger input variable $D_{p_{\\mathrm{T}}}$ ($p_{\\mathrm{T}}$ "
                      "dispersion) for Run 2 events, comparing signal and background samples after the "
                      "Lund-plane correction."),
    dict(key="allSoftDropMassAK8", section="wnae", image="wnae_inputs/allSoftDropMassAK8_log_raw_ratio.pdf",
         x_title=r"$m_{SD}(J)$", x_units="GeV", has_data=True, crop=(0.0, 40.0),
         description="The distribution of the WNAE tagger input variable soft-drop mass for Run 2 events, "
                      "comparing signal and background samples after the Lund-plane correction."),
    dict(key="allTau43AK8", section="wnae", image="wnae_inputs/allTau43AK8_log_raw_ratio.pdf",
         x_title=r"$\tau_{43}(J)$", x_units="", has_data=True,
         description="The distribution of the WNAE tagger input variable $\\tau_{43}$ for Run 2 events, comparing "
                      "signal and background samples after the Lund-plane correction."),
    dict(key="allTau32AK8", section="wnae", image="wnae_inputs/allTau32AK8_log_raw_ratio.pdf",
         x_title=r"$\tau_{32}(J)$", x_units="", has_data=True,
         description="The distribution of the WNAE tagger input variable $\\tau_{32}$ for Run 2 events, comparing "
                      "signal and background samples after the Lund-plane correction."),
    dict(key="allTau21AK8", section="wnae", image="wnae_inputs/allTau21AK8_log_raw_ratio.pdf",
         x_title=r"$\tau_{21}(J)$", x_units="", has_data=True,
         description="The distribution of the WNAE tagger input variable $\\tau_{21}$ for Run 2 events, comparing "
                      "signal and background samples after the Lund-plane correction."),
    dict(key="allEcfD2b1AK8", section="wnae", image="wnae_inputs/allEcfD2b1AK8_log_raw_ratio.pdf",
         x_title=r"$D_{2}^{\beta=1}(J)$", x_units="", has_data=True,
         description="The distribution of the WNAE tagger input variable $D_{2}^{\\beta=1}$ for Run 2 events, "
                      "comparing signal and background samples after the Lund-plane correction."),
    dict(key="allEcfC2b1AK8", section="wnae", image="wnae_inputs/allEcfC2b1AK8_log_raw_ratio.pdf",
         x_title=r"$C_{2}^{\beta=1}(J)$", x_units="", has_data=True,
         description="The distribution of the WNAE tagger input variable $C_{2}^{\\beta=1}$ for Run 2 events, "
                      "comparing signal and background samples after the Lund-plane correction."),
    dict(key="allEcfC2b2AK8", section="wnae", image="wnae_inputs/allEcfC2b2AK8_log_raw_ratio.pdf",
         x_title=r"$C_{2}^{\beta=2}(J)$", x_units="", has_data=True,
         description="The distribution of the WNAE tagger input variable $C_{2}^{\\beta=2}$ for Run 2 events, "
                      "comparing signal and background samples after the Lund-plane correction."),
    dict(key="allEcfN2b1AK8", section="wnae", image="wnae_inputs/allEcfN2b1AK8_log_raw_ratio.pdf",
         x_title=r"$N_{2}^{\beta=1}(J)$", x_units="", has_data=True,
         description="The distribution of the WNAE tagger input variable $N_{2}^{\\beta=1}$ for Run 2 events, "
                      "comparing signal and background samples after the Lund-plane correction."),
    dict(key="DeltaEta01GoodJetsAK8", section="dnn", image="dnn_inputs/DeltaEta01GoodJetsAK8_log_raw_wdata_ratio.pdf",
         x_title=r"$\Delta\eta(J_{1},J_{2})$", x_units="", has_data=True,
         description="The distribution of the event-level DNN input variable $\\Delta\\eta(J_{1},J_{2})$ for "
                      "Run 2 events, comparing signal and background samples."),
    dict(key="GoodJetsAK80_MTMETLog", section="dnn", image="dnn_inputs/GoodJetsAK80_MTMETLog_log_raw_wdata_ratio.pdf",
         x_title=r"$m_{T}(J_{1},p_{T}^{miss})$", x_units="GeV", has_data=True,
         description="The distribution of the event-level DNN input variable transverse mass of the leading AK8 "
                      "jet and $p_{\\mathrm{T}}^{\\mathrm{miss}}$ for Run 2 events, comparing signal and "
                      "background samples."),
    dict(key="DeltaR01GoodJetsAK8", section="dnn", image="dnn_inputs/DeltaR01GoodJetsAK8_log_raw_wdata_ratio.pdf",
         x_title=r"$\Delta R(J_{1},J_{2})$", x_units="", has_data=True,
         description="The distribution of the event-level DNN input variable $\\Delta R(J_{1},J_{2})$ for Run 2 "
                      "events, comparing signal and background samples."),
    dict(key="LundJetPlaneZ01GoodJetsAK8", section="dnn", image="dnn_inputs/LundJetPlaneZ01GoodJetsAK8_log_raw_wdata_ratio.pdf",
         x_title=r"$z(J_{1},J_{2})$", x_units="", has_data=True,
         description="The distribution of the event-level DNN input variable Lund-plane $z$ observable for the "
                      "two leading AK8 jets for Run 2 events, comparing signal and background samples."),
]


def cache_paths(section: str, key: str):
    """(background_group_path_fmt, data_path, signal_path_fmt) for one section."""
    if section == "nminus1":
        return (
            f"nminus1/shapes/groups/{{proc}}/{key}",
            None,
            f"nminus1/shapes/signals/{{signal}}/{key}",
        )
    return (
        f"{section}/raw/groups/Run2/{{proc}}/{key}",
        f"{section}/raw/data/Run2/observed/{key}",
        f"{section}/raw/signals/Run2/{{signal}}/{key}",
    )


def crop_hist(hist: dict, crop) -> dict:
    """Keep only bins whose center falls within [xmin, xmax]."""
    if crop is None:
        return hist
    xmin, xmax = crop
    keep = [i for i, (lo, hi) in enumerate(hist["x_edges"]) if xmin <= 0.5 * (lo + hi) <= xmax]
    out = dict(hist)
    out["x_edges"] = [hist["x_edges"][i] for i in keep]
    out["y"] = [hist["y"][i] for i in keep]
    out["dy"] = [hist["dy"][i] for i in keep]
    return out


def build_table(reader: RootFileReader, figure_number: int, spec: dict) -> Table:
    key = spec["key"]
    bkg_fmt, data_path, sig_fmt = cache_paths(spec["section"], key)

    table = Table(f"Additional Figure {figure_number}")
    table.description = spec["description"]
    table.location = f"Data from Additional Figure {figure_number}."
    table.keywords["observables"] = ["N"]
    table.keywords["reactions"] = ["P P --> JET JET"]
    table.keywords["cmenergies"] = [13000.0]
    table.add_image(os.path.join(HERE, "assets_v2", spec["image"]))

    x_var = None
    for proc in PROC_ORDER:
        h = crop_hist(reader.read_hist_1d(bkg_fmt.format(proc=proc)), spec.get("crop"))
        if x_var is None:
            x_var = Variable(spec["x_title"], is_independent=True, is_binned=True, units=spec["x_units"])
            x_var.values = h["x_edges"]
            table.add_variable(x_var)
        v = Variable(PROC_LABEL[proc], is_independent=False, is_binned=False, units="Events")
        v.values = h["y"]
        unc = Uncertainty("stat", is_symmetric=True)
        unc.values = h["dy"]
        v.add_uncertainty(unc)
        table.add_variable(v)

    if spec["has_data"] and data_path is not None:
        h = crop_hist(reader.read_hist_1d(data_path), spec.get("crop"))
        v = Variable("Data", is_independent=False, is_binned=False, units="Events")
        v.values = h["y"]
        unc = Uncertainty("stat", is_symmetric=True)
        unc.values = h["dy"]
        v.add_uncertainty(unc)
        table.add_variable(v)

    for signal_name, signal_label in SIGNALS:
        h = crop_hist(reader.read_hist_1d(sig_fmt.format(signal=signal_name)), spec.get("crop"))
        v = Variable(signal_label, is_independent=False, is_binned=False, units="Events")
        v.values = h["y"]
        table.add_variable(v)

    return table


def main() -> None:
    if not os.path.isfile(CACHE):
        raise SystemExit(f"Merged cache not found: {CACHE}\nRun make_supplementary_plots.py first.")

    reader = RootFileReader(CACHE)
    submission = Submission()
    for i, spec in enumerate(TABLE_SPECS, start=1):
        print(f"[{i}/{len(TABLE_SPECS)}] {spec['section']}/{spec['key']}")
        table = build_table(reader, i, spec)
        submission.add_table(table)

    submission.create_files(OUTPUT_DIR, remove_old=True)
    print(f"\n[DONE] wrote HEPData submission under {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
