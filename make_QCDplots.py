#!/usr/bin/env python3
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep  # Added for CMS styling
from matplotlib.backends.backend_pdf import PdfPages

# Apply CMS Style globally
plt.style.use(hep.style.CMS)

# ==========================
# CONFIG
# ==========================
YEARS = [2016, 2017, 2018]
BASE_PATH = "/uscms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis/output"
MODEL_TAG = "Current_Model_wp90_w2016/{year}_QCD.root"
OUTDIR = "./ABCD_hists3/plots_QCD_pdf"
os.makedirs(OUTDIR, exist_ok=True)

PNET_HISTS = ['h_MET_pre', "h_MET_pre_0SVJ", "h_MET_pre_1SVJ", "h_MET_pre_2SVJ", "h_MET_pre_3PSVJ"]
WNAE_HISTS = ['h_MET_pre', "h_MET_pre_WNAE_0SVJ", "h_MET_pre_WNAE_1SVJ", "h_MET_pre_WNAE_2SVJ", "h_MET_pre_WNAE_3PSVJ"]

LUMI_BY_YEAR = {"2016": "36.31", "2017": "41.48", "2018": "59.83"}
XLIM = (0, 2000)
YFLOOR = 1e-4 

# ==========================
# HELPER FUNCTIONS
# ==========================

def get_hist(rootfile, histname):
    if histname not in rootfile:
        print(f"WARNING: {histname} not found")
        return None, None
    h = rootfile[histname]
    values = h.values()
    edges = h.axes[0].edges()
    return edges, values

def apply_cms_decorations(ax, year, is_data=False):
    """Applies the CMS labels and standard grid styling."""
    LUMI_BY_YEAR = {"2016": "36.31", "2017": "41.48", "2018": "59.83"}
    lumi = LUMI_BY_YEAR.get(str(year), "")
    
    # 1. CMS Label (Top left/right)
    hep.cms.label(ax=ax, label="Preliminary", data=is_data, lumi=lumi, year=year)
    
    # 2. Scale and Limits
    ax.set_yscale("log")
    ax.set_ylim(2e-4, 1e5)
    ax.set_xlim(*XLIM)
    
    # 3. Label Alignment (CMS style: labels at the end of axes)
    ax.set_xlabel("MET [GeV]", loc='right')
    ax.set_ylabel("Events", loc='top')
    
    # 4. Grid Styling
    ax.grid(True, which='major', linestyle='-', alpha=0.3)
    ax.grid(True, which='minor', linestyle=':', alpha=0.2)
    
    # 5. Ticks (Ensuring they point inward and appear on all sides)
    ax.tick_params(which='both', direction='in', top=True, right=True, length=6)
    ax.tick_params(which='minor', length=3)
    ax.grid(True, which='both', linestyle=':', alpha=0.6)

def add_hist_page(pdf, edges, values, title, year):
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Use hep.histplot for that high-energy physics look
    hep.histplot(values, edges, ax=ax, histtype='step', linewidth=2, color='blue', label=title)
    #clean_label = name.replace("h_MET_pre_", "").replace("h_MET_pre", "Inclusive")
    ax.legend(loc='upper right', fontsize=12, frameon=True)
    #ax.set_title(title, fontsize=16, pad=20)
    apply_cms_decorations(ax, year)
    
    pdf.savefig(fig)
    plt.close(fig)

def add_overlay_page(pdf, hists_dict, year, label):
    fig, ax = plt.subplots(figsize=(14, 8))
    
    for name, (edges, values) in hists_dict.items():
        if values is None: continue
        # Clean up label for legend
        clean_label = name.replace("h_MET_pre_", "").replace("h_MET_pre", "Inclusive")
        hep.histplot(values, edges, ax=ax, histtype='step', label=clean_label, linewidth=1.5)
        
    apply_cms_decorations(ax, year)
    ax.legend(loc='upper right', fontsize=12, frameon=True)
    
    pdf.savefig(fig)
    plt.close(fig)

# ==========================
# MAIN LOOP
# ==========================
for year in YEARS:
    filepath = os.path.join(BASE_PATH, MODEL_TAG.format(year=year))
    print(f"Processing Year {year}: {filepath}")
    
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        continue

    pdf_path = os.path.join(OUTDIR, f"{year}_MET_plots.pdf")
    with uproot.open(filepath) as f, PdfPages(pdf_path) as pdf:

        # ---- PNET ----
        pnet_data = {}
        for hname in PNET_HISTS:
            edges, values = get_hist(f, hname)
            if edges is not None:
                add_hist_page(pdf, edges, values, f"PNET Selection: {hname}", year)
                pnet_data[hname] = (edges, values)
        add_overlay_page(pdf, pnet_data, year, "PNET")

        # ---- WNAE ----
        wnae_data = {}
        for hname in WNAE_HISTS:
            edges, values = get_hist(f, hname)
            if edges is not None:
                add_hist_page(pdf, edges, values, f"WNAE Selection: {hname}", year)
                wnae_data[hname] = (edges, values)
        add_overlay_page(pdf, wnae_data, year, "WNAE")

    print(f"Saved CMS-style PDF: {pdf_path}")

print("\nAll done.")