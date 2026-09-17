#!/usr/bin/env python3
"""Plot WNAE AUC heatmaps from the JSON tables in supplemetry/AUCs."""

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, BASE)

import makeAUCtables_pnet as auc_mod

AUC_DIR = os.path.join(HERE, "AUCs")
OUT_DIR = os.path.join(HERE, "assets_v2", "auc")


AXIS_LABELS = {
    "r_inv": r"$r_{\mathrm{inv}}$",
    "m_dark_gev": r"$m_{\mathrm{dark}}\ [\mathrm{GeV}]$",
    "lambda": r"$\lambda$",
}


BACKGROUND_LABELS = {
    "qcd_multijet": "QCD",
    "ttbar": "ttbar",
    "single_top": "single top",
    "w_to_lnu": r"W$\rightarrow l\nu$",
    "z_to_nunu": r"Z$\rightarrow\nu\nu$",
}


def tick_fmt(value):
    return f"{value:g}"


def output_name(path):
    stem = os.path.splitext(os.path.basename(path))[0]
    stem = re.sub(r"[^A-Za-z0-9_]+", "_", stem)
    return os.path.join(OUT_DIR, f"wnae_auc_{stem}")


def plot_one(path):
    with open(path) as fin:
        data = json.load(fin)

    rows = data["axes"]["rows"]["values"]
    cols = data["axes"]["columns"]["values"]
    col_name = data["axes"]["columns"]["name"]
    raw = data["auc"]
    bkg_key = os.path.basename(path).split("_pt_")[0]
    bkg_label = BACKGROUND_LABELS.get(bkg_key, bkg_key.replace("_", " "))
    pt = data["jet_pt_bin_gev"]
    pt_label = f"{pt['min_inclusive']:g}-{pt['max_exclusive']:g} GeV" if pt["max_exclusive"] is not None else f"{pt['min_inclusive']:g}+ GeV"
    z_label = f"WNAE AUC vs. {bkg_label}"

    auc_mod.plot_auc_grid(
        mMed_vals=[int(v) if float(v).is_integer() else v for v in rows],
        rinv_vals=cols,
        raw=raw,
        z_label=z_label,
        out_prefix=output_name(path),
        vmin=0.5,
        vmax=1.0,
        cmap="viridis",
        training_samples=None,
        y_label=AXIS_LABELS.get(col_name, col_name),
        y_tick_fmt=tick_fmt,
    )


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    paths = sorted(
        os.path.join(AUC_DIR, name)
        for name in os.listdir(AUC_DIR)
        if name.endswith(".json")
    )
    for path in paths:
        plot_one(path)
    print(f"[OK] wrote {len(paths)} WNAE AUC heatmap sets under {OUT_DIR}")


if __name__ == "__main__":
    main()
