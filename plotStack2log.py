#!/usr/bin/env python3
import os
import math
import optparse
import glob

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gStyle.SetOptStat(0)

import utils.DataSetInfo as info
# import utils.CMS_lumi as CMS_lumi


'''
Example:
python3 plotStack2.py \
    -d /uscms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis/output/Current_Model_wp90_w2016 \
    -o paper_stacks \
    --years 2016,2017,2018 \
    --cuts _pre_,_pre_WNAE_ \
    --ratio

New signal options:
    --signals "SVJ,signal"
    --signal-scale 10


    --signal-norm

python3 plotStack2log.py    -d /uscms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis/output/Current_Model_wp90_w2016   -o paper_stackslog    --years 2016,2017,2018   --cuts _pre_   --signals "m600_d20_r0p3,m2000_d20_r0p1,m2000_d20_r0p3,m2000_d20_r0p7,m4000_d20_r0p3"
'''

LUMI_BY_YEAR = {"2016": "36.31", "2017": "42.07", "2018": "59.56"}

PROC_LABEL = {
    "QCD": "QCD",
    "TTJets": "t#bar{t}+jets",
    "ZJetsToNuNu": "Z#rightarrow#nu#nu+jets",
    "WJetsToLNu": "W#rightarrowl#nu+jets",
    "ST": "Single top",
}

PROC_COLOR = {
    "QCD": ROOT.TColor.GetColor("#9c9ca1"),
    "TTJets": ROOT.TColor.GetColor("#7a21dd"),
    "WJetsToLNu": ROOT.TColor.GetColor("#e42536"),
    "ZJetsToNuNu": ROOT.TColor.GetColor("#f89c20"),
    "ST": ROOT.TColor.GetColor("#5790fc"),
}

PROC_HATCH = {
    "QCD": 1001,
    "TTJets": 1001,
    "WJetsToLNu": 1001,
    "ZJetsToNuNu": 1001,
    "ST": 1001,
}

LEGEND_ORDER = ["QCD", "TTJets", "WJetsToLNu", "ZJetsToNuNu", "ST"]
STACK_ORDER = ["ST", "ZJetsToNuNu", "WJetsToLNu", "TTJets", "QCD"]

SVJ_BINS = ["0SVJ", "1SVJ", "2SVJ", "3SVJ", "3PSVJ"]

# ---------------------------------------------------------
# Signal overlay config
# ---------------------------------------------------------
SIGNAL_LINE_COLORS = [
    ROOT.TColor.GetColor("#92dadd"),  # cyan
    ROOT.TColor.GetColor("#6b3e26"),  # dark brown
    ROOT.TColor.GetColor("#0b3d02"),  # dark green
    ROOT.TColor.GetColor("#228833"),  # green dashed
    ROOT.TColor.GetColor("#1f4e79"),  # dark blue
]

SIGNAL_LINE_STYLES = [
    1,  # solid
    3,  # dotted
    1,  # solid
    2,  # dashed
    1,  # solid
]


def style_signal_hist(h, idx):
    h.SetLineColor(SIGNAL_LINE_COLORS[idx % len(SIGNAL_LINE_COLORS)])
    h.SetLineWidth(2)
    h.SetLineStyle(SIGNAL_LINE_STYLES[idx % len(SIGNAL_LINE_STYLES)])
    h.SetFillStyle(0)
    h.SetMarkerSize(0)
    return h

BKG_FILE_TOKENS = [
    "_Data.root",
    "_ST.root",
    "_TTJets.root",
    "_ZJets.root",
    "_WJets.root",
    "_QCD.root",
]
def get_signal_label_from_filename(fname, year="", hem_period=""):
    base = os.path.basename(fname).replace(".root", "")

    prefix = f"{year}{hem_period}_"
    if base.startswith(prefix):
        base = base[len(prefix):]

    mphi = None
    rinv = None

    for part in base.split("_"):
        if part.startswith("m") and part[1:].isdigit():
            mphi = part[1:]

        elif part.startswith("r"):
            rinv = part[1:].replace("p", ".")

    if mphi is not None and rinv is not None:
        return rf"m_{{#phi}} = {mphi} GeV, r_{{inv}} = {rinv}"

    return base

def discover_signal_files(path, year="2018", HEMPeriod="", patterns_csv=""):
    """
    If --signals is provided, use those comma-separated substrings to filter files.
    Otherwise auto-pick all YEAR[_HEM]_*.root that are not standard bg/data files.
    """
    pattern = os.path.join(path, f"{year}{HEMPeriod}_*.root")
    all_files = sorted(glob.glob(pattern))

    requested_patterns = [x.strip() for x in patterns_csv.split(",") if x.strip()]
    out = []

    for fullpath in all_files:
        base = os.path.basename(fullpath)

        if any(tok in base for tok in BKG_FILE_TOKENS):
            continue

        if requested_patterns and not any(pat in base for pat in requested_patterns):
            continue

        out.append(base)

    return out


def get_signal_data(path, scale=1.0, year="2018", HEMPeriod="", signal_patterns=""):
    sgData = []
    signal_files = discover_signal_files(
        path,
        year=year,
        HEMPeriod=HEMPeriod,
        patterns_csv=signal_patterns,
    )

    for fname in signal_files:
        try:
            sgData.append(
                info.DataSetInfo(
                    basedir=path,
                    fileName=fname,
                    label=get_signal_label_from_filename(fname, year=year, hem_period=HEMPeriod),
                    scale=scale,
                )
            )
        except Exception as e:
            print(f"[WARN] Could not build signal dataset for {fname}: {e}")

    if sgData:
        print("[INFO] Signals found:")
        for s in sgData:
            print(f"   - {s.fileName}")
    else:
        print("[INFO] No signal files found for overlay.")

    return sgData


def getData(path, scale=1.0, year="2018", HEMPeriod=""):
    Data = [
        info.DataSetInfo(
            basedir=path,
            fileName=f"{year}{HEMPeriod}_Data.root",
            sys=-1.0,
            label="Data",
            scale=scale,
        ),
    ]

    bgData = [
        info.DataSetInfo(
            basedir=path,
            fileName=f"{year}{HEMPeriod}_ST.root",
            label="Single top",
            scale=scale,
            color=(ROOT.kRed + 1),
        ),
        info.DataSetInfo(
            basedir=path,
            fileName=f"{year}{HEMPeriod}_TTJets.root",
            label="t#bar{t}+jets",
            scale=scale,
            color=(ROOT.kBlue - 6),
        ),
        info.DataSetInfo(
            basedir=path,
            fileName=f"{year}{HEMPeriod}_ZJets.root",
            label="Z#rightarrow#nu#nu+jets",
            scale=scale,
            color=(ROOT.kGray + 1),
        ),
        info.DataSetInfo(
            basedir=path,
            fileName=f"{year}{HEMPeriod}_WJets.root",
            label="W#rightarrowl#nu+jets",
            scale=scale,
            color=(ROOT.kYellow + 1),
        ),
        info.DataSetInfo(
            basedir=path,
            fileName=f"{year}{HEMPeriod}_QCD.root",
            label="QCD",
            scale=scale,
            color=(ROOT.kGreen + 1),
            ),
    ]

    sgData = []
    return Data, sgData, bgData


def dataset_exists(basedir, filename):
    return os.path.exists(os.path.join(basedir, filename))


def getData_with_missing_data_ok(path, scale=1.0, year="2018", HEMPeriod="", signal_patterns=""):
    Data, _, bgData = getData(path, scale=scale, year=year, HEMPeriod=HEMPeriod)

    if not dataset_exists(path, f"{year}{HEMPeriod}_Data.root"):
        print(f"[INFO] Missing data file for {year}{HEMPeriod}; proceeding with MC only.")
        Data = []

    existing_bg = []
    for sample in bgData:
        if dataset_exists(path, sample.fileName):
            existing_bg.append(sample)
        else:
            print(f"[INFO] Missing background file: {sample.fileName}")

    sgData = get_signal_data(
        path,
        scale=scale,
        year=year,
        HEMPeriod=HEMPeriod,
        signal_patterns=signal_patterns,
    )

    existing_sg = []
    for sample in sgData:
        if dataset_exists(path, sample.fileName):
            existing_sg.append(sample)
        else:
            print(f"[INFO] Missing signal file: {sample.fileName}")

    return Data, existing_sg, existing_bg


def ensure_dataset_labels(datasets):
    for d in datasets:
        if not hasattr(d, "label"):
            try:
                d.label = d.legEntry()
            except Exception:
                d.label = getattr(d, "fileName", "sample")


def makeDirs(outdir, cut):
    target = os.path.join(outdir, cut[1:] if cut.startswith("_") else cut)
    os.makedirs(target, exist_ok=True)
    return target


def divisorGenerator(n):
    large_divisors = []
    for i in range(1, int(math.sqrt(n) + 1)):
        if n % i == 0:
            yield i
            if i * i != n:
                large_divisors.append(n // i)
    for divisor in reversed(large_divisors):
        yield divisor


def find_nearest(values, target):
    values = list(values)
    return min(values, key=lambda x: abs(x - target))


def rebinCalc(nBins, target):
    if nBins <= 0:
        return 1
    rebinFloat = nBins / float(target)
    allDivs = list(divisorGenerator(nBins))
    if not allDivs:
        return 1
    return max(1, int(find_nearest(allDivs, rebinFloat)))


def clone_hist(h, name):
    out = h.Clone(name)
    out.SetDirectory(0)
    return out

def setup_canvas(isRatio=False, isLogY=True):
    if isRatio:
        c = ROOT.TCanvas("c", "c", 800, 800)

        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.28, 1.0, 1.0)
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.00, 1.0, 0.28)

        pad1.SetBottomMargin(0.02)
        pad1.SetLeftMargin(0.16)
        pad1.SetRightMargin(0.05)
        pad1.SetTopMargin(0.08)
        pad1.SetTicks(1, 1)
        pad1.SetLogy(isLogY)
        pad1.Draw()

        pad2.SetTopMargin(0.03)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.16)
        pad2.SetRightMargin(0.05)
        pad2.SetTicks(1, 1)
        pad2.SetGridy()
        pad2.Draw()

        return c, pad1, pad2

    c = ROOT.TCanvas("c", "c", 800, 800)
    c.cd()

    ROOT.gPad.SetLeftMargin(0.16)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1, 1)
    ROOT.gPad.SetLogy(isLogY)

    return c, None, None

def make_legend(with_data=True, n_signal=0):
    """
    MC-only style is left untouched.
    For data plots, this function is not used directly anymore.
    """

    # Original MC-only style untouched
    leg = ROOT.TLegend(0.19, 0.730, 0.93, 0.895)
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.030)
    leg.SetTextFont(42)
    leg.SetColumnSeparation(0.04)
    leg.SetEntrySeparation(0.028)

    return leg

def make_data_legends():
    """
    Data plot only:
    Use two independent legends so entries do not overlap.

    Left legend:
      Data, QCD, ttbar, W, Z, Single top, MC unc.

    Right legend:
      signal overlays
    """

    # Higher left legend
    leg_bkg = ROOT.TLegend(0.19, 0.665, 0.50, 0.900)
    leg_bkg.SetNColumns(1)
    leg_bkg.SetBorderSize(0)
    leg_bkg.SetFillStyle(0)
    leg_bkg.SetTextSize(0.030)
    leg_bkg.SetTextFont(42)
    leg_bkg.SetEntrySeparation(0.028)

    # Higher right legend
    leg_sig = ROOT.TLegend(0.58, 0.720, 0.93, 0.900)
    leg_sig.SetNColumns(1)
    leg_sig.SetBorderSize(0)
    leg_sig.SetFillStyle(0)
    leg_sig.SetTextSize(0.030)
    leg_sig.SetTextFont(42)
    leg_sig.SetEntrySeparation(0.033)

    return leg_bkg, leg_sig
    
def get_paper_vars():
    return {
        "MET": ("p_{T}^{miss} [GeV]", 500, 200.0, 2000.0),
        "HT": ("H_{T} [GeV]", 500, 0.0, 5000.0),
        "ST": ("S_{T} [GeV]", 500, 1300.0, 5000.0),
        "mT": ("m_{T} [GeV]", 500, 0.0, 6000.0),
        "dPhiMinjMETAK8": ("#Delta#phi_{min}(J, p_{T}^{miss})", 100, 0.0, 2.0),
        "nsvjJetsAK8": ("Number of SVJ AK8 jets", 20, 0.0, 5.0),
        "dnnEventClassScore": ("Event classifier score", 100, 0.0, 1.0),
        "j1PtAK8": ("p_{T}(J_{1}) [GeV]", 400, 0.0, 4000.0),
        "j1SoftDropMassAK8": ("m_{SD}(J_{1}) [GeV]", 200, 0.0, 900.0),
        "j1Tau21AK8": ("#tau_{21}(J_{1})", 40, 0.0, 1.0),
    }


def get_proc_key(dataset):
    label = dataset.legEntry().lower() if hasattr(dataset, "legEntry") else str(dataset).lower()
    fname = getattr(dataset, "fileName", "").lower()
    full = f"{label} {fname}"

    if "qcd" in full:
        return "QCD"
    if "ttjets" in full or "ttbar" in full or "t#bar{t}" in full or full.endswith("_tt") or " tt" in full:
        return "TTJets"
    if "zjets" in full or "znunu" in full or "z#rightarrow#nu#nu" in full or "zto" in full:
        return "ZJetsToNuNu"
    if "wjets" in full or "w+jets" in full or "wtolnu" in full or "wto" in full:
        return "WJetsToLNu"
    if "single top" in full or "_st" in full or full.startswith("st") or "/st" in full:
        return "ST"
    return "QCD"


def style_mc_hist(h, proc_key):
    color = PROC_COLOR.get(proc_key, ROOT.kGray)
    h.SetFillColor(color)
    h.SetFillStyle(1001)

    # Cleaner paper-style stack: no black outline on each component
    h.SetLineColor(color)
    h.SetLineWidth(0)

    h.SetMarkerSize(0)
    return h

def style_data_hist(h):
    h.SetLineColor(ROOT.kBlack)
    h.SetMarkerColor(ROOT.kBlack)
    h.SetMarkerStyle(20)
    h.SetMarkerSize(1.0)
    h.SetLineWidth(2)
    return h


def get_summed_hist(dset, histoName, rebinx, xmin, xmax, fill=False):
    h = None

    for svj in SVJ_BINS:
        full_name = f"{histoName}{svj}"

        try:
            h_tmp = dset.getHisto(
                full_name,
                rebinx=rebinx,
                xmin=xmin,
                xmax=xmax,
                fill=fill,
                showEvents=False,
            )
        except Exception as e:
            print(f"[WARN] Failed to get {full_name} from {dset.fileName}: {e}")
            continue

        if not h_tmp or not hasattr(h_tmp, "GetNbinsX"):
            continue

        if h is None:
            h = clone_hist(h_tmp, f"{full_name}_sum")
        else:
            h.Add(h_tmp)

    return h


def build_mc_stack(bgData, histoName, rebinx, xmin, xmax):
    stack = ROOT.THStack("hs", "hs")
    total_mc = None
    proc_hists = {}

    grouped = {key: [] for key in STACK_ORDER}
    for dset in bgData:
        grouped.setdefault(get_proc_key(dset), []).append(dset)

    for proc_key in STACK_ORDER:
        proc_hist = None

        for i, dset in enumerate(grouped.get(proc_key, [])):
            h = get_summed_hist(
                dset,
                histoName=histoName,
                rebinx=rebinx,
                xmin=xmin,
                xmax=xmax,
                fill=True,
            )
            if h is None:
                continue

            h = clone_hist(h, f"mc_{proc_key}_{i}_{histoName}")
            if proc_hist is None:
                proc_hist = h
            else:
                proc_hist.Add(h)

        if proc_hist is None:
            continue

        style_mc_hist(proc_hist, proc_key)
        stack.Add(proc_hist)
        proc_hists[proc_key] = proc_hist

        if total_mc is None:
            total_mc = clone_hist(proc_hist, f"totalmc_{histoName}")
        else:
            total_mc.Add(proc_hist)

    if total_mc is not None:
        total_mc.SetFillStyle(0)
        total_mc.SetLineColor(ROOT.kBlack)
        total_mc.SetLineWidth(2)

    return stack, total_mc, proc_hists


def normalize_stack_to_unity(proc_hists, total_mc):
    """
    Normalize the stacked MC components by the total MC integral.
    This makes the stack a shape distribution instead of event yield.
    """
    if total_mc is None:
        return

    integral = total_mc.Integral()
    if integral <= 0:
        return

    for h in proc_hists.values():
        h.Scale(1.0 / integral)

    total_mc.Scale(1.0 / integral)

def build_data_hist(dataList, histoName, rebinx, xmin, xmax):
    if not dataList:
        return None

    data_hist = None
    for i, dset in enumerate(dataList):
        h = get_summed_hist(
            dset,
            histoName=histoName,
            rebinx=rebinx,
            xmin=xmin,
            xmax=xmax,
            fill=False,
        )
        if h is None:
            continue

        h = clone_hist(h, f"data_{i}_{histoName}")
        if data_hist is None:
            data_hist = h
        else:
            data_hist.Add(h)

    if data_hist is not None:
        style_data_hist(data_hist)

    return data_hist


def build_signal_hists(sgData, histoName, rebinx, xmin, xmax, signal_scale=1.0, normalize=False):
    signal_hists = []

    if not sgData:
        return signal_hists

    for i, dset in enumerate(sgData):
        h = get_summed_hist(
            dset,
            histoName=histoName,
            rebinx=rebinx,
            xmin=xmin,
            xmax=xmax,
            fill=False,
        )
        if h is None:
            continue

        h = clone_hist(h, f"sig_{i}_{histoName}")

        if signal_scale != 1.0:
            h.Scale(signal_scale)

        if normalize:
            integral = h.Integral()
            if integral > 0:
                h.Scale(1.0 / integral)

        style_signal_hist(h, i)
        signal_hists.append((dset, h))

    return signal_hists

def build_signal_hists_combined_by_label(
    sgData,
    histoName,
    rebinx,
    xmin,
    xmax,
    signal_scale=1.0,
    normalize=False,
):
    """
    For combined-years MC-only plots.

    If the same signal exists in 2016, 2017, and 2018, this sums those
    histograms into one signal curve before optional normalization.
    """
    if not sgData:
        return []

    grouped = {}

    for dset in sgData:
        label = getattr(dset, "label", getattr(dset, "fileName", "Signal"))

        h = get_summed_hist(
            dset,
            histoName=histoName,
            rebinx=rebinx,
            xmin=xmin,
            xmax=xmax,
            fill=False,
        )

        if h is None:
            continue

        h = clone_hist(h, f"sig_combined_tmp_{label}_{histoName}")

        if label not in grouped:
            grouped[label] = h
        else:
            grouped[label].Add(h)

    signal_hists = []

    for i, label in enumerate(sorted(grouped.keys())):
        h = grouped[label]
        h = clone_hist(h, f"sig_combined_{i}_{histoName}")

        if signal_scale != 1.0:
            h.Scale(signal_scale)

        if normalize:
            integral = h.Integral()
            if integral > 0:
                h.Scale(1.0 / integral)

        style_signal_hist(h, i)

        # Fake lightweight object so existing legend code works
        class SignalLabel:
            pass

        dset_label = SignalLabel()
        dset_label.label = label

        signal_hists.append((dset_label, h))

    return signal_hists

def make_ratio(data_hist, mc_hist, xtitle):
    ratio = clone_hist(data_hist, f"ratio_{data_hist.GetName()}")
    ratio.Divide(mc_hist)
    ratio.SetTitle("")
    ratio.SetMinimum(0.0)
    ratio.SetMaximum(2.0)
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetXaxis().SetTitle(xtitle)
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetNdivisions(505)

    ratio.GetXaxis().SetTitleSize(0.13)
    ratio.GetXaxis().SetLabelSize(0.11)
    ratio.GetXaxis().SetTitleOffset(1.05)

    ratio.GetYaxis().SetLabelSize(0.10)
    ratio.GetYaxis().SetTitleOffset(0.40)
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetLabelOffset(0.01)

    return ratio


def make_mc_unc_band(mc_hist):
    band = clone_hist(mc_hist, f"band_{mc_hist.GetName()}")
    band.SetFillColor(ROOT.kGray + 2)
    band.SetFillStyle(3004)
    band.SetMarkerSize(0)
    band.SetLineWidth(0)
    return band

def get_ymax(stack, mc_hist, data_hist=None, signal_hists=None, isLogY=True):
    ymax = 0.0

    if mc_hist:
        ymax = max(ymax, mc_hist.GetMaximum())

    if data_hist:
        ymax = max(ymax, data_hist.GetMaximum())

    if signal_hists:
        for _, h in signal_hists:
            ymax = max(ymax, h.GetMaximum())

    if ymax <= 0:
        ymax = 1.0

    if isLogY:
        if data_hist is not None:
            return ymax * 900.0   # data only: raise roof more
        else:
            return ymax * 60.0    # MC-only unchanged

    if data_hist is not None:
        return ymax * 2.2

    return ymax * 1.5



def draw_common_signal_text():
    """
    Common signal parameters shown below the legend.
    Keep this short so it does not run off the canvas.
    """
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.030) # Matches legend text size precisely
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextAlign(12)  # Left-aligned, vertically centered

    latex.DrawLatex(
        0.605,  # Shifted left to align cleanly with the m_phi labels
        0.700,  # Moved higher up to close the gap right below the legend matrix
        "m_{dark} = 20 GeV, #lambda = 1"
    )
def draw_cms_label(year=None, extra_text="Preliminary", show_lumi=True):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextColor(ROOT.kBlack)

    # CMS
    latex.SetTextFont(61)
    latex.SetTextSize(0.060)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.16, 0.94, "CMS")

    # Preliminary should be close, Simulation needs more space
    latex.SetTextFont(52)
    latex.SetTextSize(0.040)
    latex.SetTextAlign(11)

    if extra_text == "Simulation":
        latex.DrawLatex(0.290, 0.94, extra_text)   # farther from CMS
    else:
        latex.DrawLatex(0.280, 0.94, extra_text)   # close to CMS for Preliminary

    # Upper-right label
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(31)

    if show_lumi and year is not None:
        lumi_text = f"{LUMI_BY_YEAR.get(year, '')} fb^{{-1}} (13 TeV)"
    else:
        lumi_text = "(13 TeV)"

    latex.DrawLatex(0.95, 0.94, lumi_text)

def save_canvas(canvas, outbase):
    canvas.SaveAs(outbase + ".pdf")
    canvas.SaveAs(outbase + ".png")


def plot_paper_stack(data_tuple, histoName, totalBin, outputPath, xTitle, yTitle="Events",
                     isLogY=True, xmin=999.9, xmax=-999.9, year="2018", do_ratio=True,
                     signal_scale=1.0, signal_norm=False,
                     combine_signal_years=False):
    ROOT.TH1.AddDirectory(False)

    dataList, bgData, sgData = data_tuple
    if not bgData:
        print(f"[WARN] No background datasets available for {histoName}")
        return

    do_ratio = bool(do_ratio and dataList)

    rebinx = rebinCalc(totalBin, 40)
    canvas, pad1, pad2 = setup_canvas(isRatio=do_ratio, isLogY=isLogY)

    if do_ratio:
        pad1.cd()

    stack, mc_hist, proc_hists = build_mc_stack(bgData, histoName, rebinx, xmin, xmax)
    data_hist = build_data_hist(dataList, histoName, rebinx, xmin, xmax)
    if combine_signal_years:
        signal_hists = build_signal_hists_combined_by_label(
            sgData,
            histoName=histoName,
            rebinx=rebinx,
            xmin=xmin,
            xmax=xmax,
            signal_scale=signal_scale,
            normalize=signal_norm,
        )
    else:
        signal_hists = build_signal_hists(
            sgData,
            histoName=histoName,
            rebinx=rebinx,
            xmin=xmin,
            xmax=xmax,
            signal_scale=signal_scale,
            normalize=signal_norm,
        )

    # For MC-only shape plots:
    # normalize the stacked MC to unit area if signal_norm=True.
    if data_hist is None and signal_norm:
        normalize_stack_to_unity(proc_hists, mc_hist)

        # Rebuild the stack after scaling the individual components
        stack = ROOT.THStack("hs_norm", "hs_norm")
        for proc_key in STACK_ORDER:
            if proc_key in proc_hists:
                stack.Add(proc_hists[proc_key])

    if mc_hist is None:
        print(f"[WARN] Could not build MC histogram for {histoName}")
        canvas.Close()
        return

    dummy = clone_hist(mc_hist, f"dummy_{histoName}")
    dummy.Reset()
    dummy.SetStats(0)
    dummy.SetTitle("")
    #dummy.GetYaxis().SetTitle(yTitle)
    if data_hist is None and signal_norm:
        dummy.GetYaxis().SetTitle("Arbitrary units")
    else:
        dummy.GetYaxis().SetTitle(yTitle)
    dummy.GetYaxis().SetTitleSize(0.05)

    if do_ratio:
        dummy.GetXaxis().SetLabelSize(0)
        dummy.GetXaxis().SetTitleSize(0)
    else:
        dummy.GetXaxis().SetTitle(xTitle)
        dummy.GetXaxis().SetTitleSize(0.05)
        dummy.GetXaxis().SetLabelSize(0.04)

    dummy.GetYaxis().SetTitleSize(0.05)
    dummy.GetYaxis().SetLabelSize(0.04)
    dummy.GetYaxis().SetLabelOffset(0.01)
    #dummy.GetYaxis().SetTitleOffset(0.95)
    dummy.GetYaxis().SetTitleOffset(1.3) # Increased from 0.95

    
    if isLogY:
        if data_hist is not None:
            dummy.SetMinimum(1e-3)
        else:
            dummy.SetMinimum(1e-5)
    else:
        dummy.SetMinimum(0.0)
    #dummy.SetMinimum(1e-4 if isLogY else 0.0)
    dummy.SetMaximum(get_ymax(stack, mc_hist, data_hist, signal_hists=signal_hists, isLogY=isLogY))
    if xmin < xmax:
        dummy.GetXaxis().SetRangeUser(xmin, xmax)
    dummy.Draw("hist")

    stack.Draw("hist same")

    if data_hist is not None:
        unc_band = make_mc_unc_band(mc_hist)
        unc_band.Draw("E2 same")
        mc_hist.Draw("hist same")
    else:
        unc_band = None

    for _, sig_hist in signal_hists:
        sig_hist.Draw("hist same")

    if data_hist is not None:
        data_hist.SetStats(0)
        data_hist.Draw("PE same")

    # ---------------------------------------------------------
    # Legend: second-plot style
    # Left column = backgrounds/data
    # Right column = signal overlays
    # ---------------------------------------------------------
    bkg_entries = []

    if data_hist is not None:
        bkg_entries.append((data_hist, "Data", "PE"))

    for proc_key in LEGEND_ORDER:
        if proc_key in proc_hists:
            bkg_entries.append((proc_hists[proc_key], PROC_LABEL[proc_key], "F"))

    if data_hist is not None and unc_band is not None:
        bkg_entries.append((unc_band, "MC unc.", "F"))

    sig_entries = []
    for dset, sig_hist in signal_hists:
        sig_label = getattr(dset, "label", getattr(dset, "fileName", "Signal"))

        if signal_norm:
            pass
        elif signal_scale != 1.0:
            sig_label += f" (#times {signal_scale:g})"

        sig_entries.append((sig_hist, sig_label, "L"))

    if data_hist is not None:
        # ---------------------------------------------------------
        # Data plots only:
        # Split legend into two independent legends to avoid overlap.
        # ---------------------------------------------------------
        leg_bkg, leg_sig = make_data_legends()

        for entry in bkg_entries:
            leg_bkg.AddEntry(*entry)

        for entry in sig_entries:
            leg_sig.AddEntry(*entry)

        leg_bkg.Draw()
        leg_sig.Draw()

    else:
        # ---------------------------------------------------------
        # MC-only plots:
        # Original two-column legend untouched.
        # ---------------------------------------------------------
        leg = make_legend(with_data=False, n_signal=len(signal_hists))

        n_rows = max(len(bkg_entries), len(sig_entries))

        for i in range(n_rows):
            if i < len(bkg_entries):
                leg.AddEntry(*bkg_entries[i])
            else:
                leg.AddEntry("", "", "")

            if i < len(sig_entries):
                leg.AddEntry(*sig_entries[i])
            else:
                leg.AddEntry("", "", "")

        leg.Draw()

    
    if data_hist is None and len(signal_hists) > 0:
        draw_common_signal_text()

    draw_cms_label(
        year=year,
        extra_text="Preliminary" if data_hist is not None else "Simulation",
        show_lumi=(data_hist is not None),
    )

    ROOT.gPad.RedrawAxis()

    if do_ratio and data_hist is not None:
        pad2.cd()
        ratio = make_ratio(data_hist, mc_hist, xTitle)
        ratio.SetStats(0)
        ratio.Draw("PE")

        ratio_unc = clone_hist(mc_hist, f"ratio_unc_{histoName}")
        for ib in range(1, ratio_unc.GetNbinsX() + 1):
            mc = mc_hist.GetBinContent(ib)
            err = mc_hist.GetBinError(ib)
            ratio_unc.SetBinContent(ib, 1.0 if mc > 0 else 0.0)
            ratio_unc.SetBinError(ib, err / mc if mc > 0 else 0.0)
        ratio_unc.SetFillColor(ROOT.kBlack)
        ratio_unc.SetFillStyle(3004)
        ratio_unc.SetMarkerSize(0)
        ratio_unc.SetLineWidth(0)
        ratio_unc.Draw("E2 same")
        ratio.Draw("PE same")

        line_xmin = xmin if xmin < xmax else ratio.GetXaxis().GetXmin()
        line_xmax = xmax if xmin < xmax else ratio.GetXaxis().GetXmax()
        line = ROOT.TLine(line_xmin, 1.0, line_xmax, 1.0)
        line.SetLineStyle(ROOT.kDashed)
        line.Draw("same")
        ROOT.gPad.RedrawAxis()

    outbase = os.path.join(outputPath, histoName)
    save_canvas(canvas, outbase)
    canvas.Close()
    print(f"[OK] Wrote {outbase}.pdf/.png")


def main():
    parser = optparse.OptionParser()
    parser.add_option("-d", "--dataset", dest="dataset", default="output/Current_Model_wp90",
                      help="dataset base dir (contains YEAR[_HEM]_*.root)")
    parser.add_option("-o", "--outdir", dest="outdir", default="paper_stacks",
                      help="output directory")
    parser.add_option("--years", dest="years", default="2016,2017,2018",
                      help="comma-separated years")
    parser.add_option("--hem-period", dest="hemPeriod", default="",
                      help="optional HEM suffix, e.g. _PreHEM or _PostHEM")
    parser.add_option("--cuts", dest="cuts", default="_pre_,_pre_WNAE_",
                      help="comma-separated cuts to plot")
    parser.add_option("--ratio", action="store_true", dest="ratio", default=False,
                      help="draw data/MC ratio panel when data exists")
    parser.add_option("--nolog", action="store_true", dest="nolog", default=False,
                      help="turn off log y-axis")

    parser.add_option("--signals", dest="signals", default="",
                      help="comma-separated substrings to select signal ROOT files, "
                           "e.g. 'SVJ_mMed-200,SigX'. If empty, auto-discover all "
                           "non-bg YEAR[_HEM]_*.root files.")
    parser.add_option("--signal-scale", dest="signal_scale", type="float", default=1.0,
                      help="multiply signal histograms by this factor before drawing")
    parser.add_option("--signal-norm", action="store_true", dest="signal_norm", default=False,
                      help="normalize each signal shape to unit area before drawing")

    parser.add_option("--no-combined-years-mc", action="store_true",
                      dest="no_combined_years_mc", default=False,
                      help="do not make combined 2016+2017+2018 MC-only signal plots")

    opts, _ = parser.parse_args()

    base_dataset = opts.dataset.rstrip("/") + "/"
    os.makedirs(opts.outdir, exist_ok=True)

    years = [y.strip() for y in opts.years.split(",") if y.strip()]
    cuts = [c.strip() for c in opts.cuts.split(",") if c.strip()]
    all_vars = get_paper_vars()

    # Store all years here for the combined Run2 MC-only plots
    all_year_bgData = []
    all_year_sgData = []

    for year in years:
        print(f"\n[INFO] Processing year {year}")

        Data, sgData, bgData = getData_with_missing_data_ok(
            base_dataset,
            scale=1.0,
            year=year,
            HEMPeriod=opts.hemPeriod,
            signal_patterns=opts.signals,
        )

        ensure_dataset_labels(Data)
        ensure_dataset_labels(bgData)
        ensure_dataset_labels(sgData)

        all_year_bgData.extend(bgData)
        all_year_sgData.extend(sgData)

        for cut in cuts:
            cut_out = makeDirs(os.path.join(opts.outdir, year), cut)
            print(f"[INFO]  cut = {cut}")

            mc_sig_out = os.path.join(cut_out, "mc_signal_only")
            os.makedirs(mc_sig_out, exist_ok=True)

            for histName, details in all_vars.items():
                xlabel, nbins, xmin, xmax = details
                full_hname = "h_" + histName + cut

                # -------------------------------------------------
                # Per-year data+MC plot
                # -------------------------------------------------
                try:
                    plot_paper_stack(
                        (Data, bgData, sgData),
                        histoName=full_hname,
                        totalBin=nbins,
                        outputPath=cut_out,
                        xTitle=xlabel,
                        yTitle="Events",
                        isLogY=(not opts.nolog),
                        xmin=xmin,
                        xmax=xmax,
                        year=year,
                        do_ratio=opts.ratio,
                        signal_scale=opts.signal_scale,
                        signal_norm=opts.signal_norm,
                    )
                except Exception as e:
                    print(f"[WARN] Failed on {full_hname}: {e}")

                # -------------------------------------------------
                # Per-year MC-only + signal shape plot
                # -------------------------------------------------
                try:
                    plot_paper_stack(
                        ([], bgData, sgData),
                        histoName=full_hname,
                        totalBin=nbins,
                        outputPath=mc_sig_out,
                        xTitle=xlabel,
                        yTitle="Events",
                        isLogY=True,
                        xmin=xmin,
                        xmax=xmax,
                        year=year,
                        do_ratio=False,
                        signal_scale=opts.signal_scale,
                        signal_norm=True,
                        combine_signal_years=False,
                    )
                except Exception as e:
                    print(f"[WARN] Failed MC+signal-only on {full_hname}: {e}")

    # -------------------------------------------------------------
    # Combined-years MC-only + signal shape plots
    # -------------------------------------------------------------
    if not opts.no_combined_years_mc:
        print("\n[INFO] Making combined-years MC-only signal plots")

        ensure_dataset_labels(all_year_bgData)
        ensure_dataset_labels(all_year_sgData)

        for cut in cuts:
            run2_cut_out = makeDirs(os.path.join(opts.outdir, "Run2"), cut)
            run2_mc_sig_out = os.path.join(run2_cut_out, "mc_signal_only")
            os.makedirs(run2_mc_sig_out, exist_ok=True)

            print(f"[INFO]  combined Run2 cut = {cut}")

            for histName, details in all_vars.items():
                xlabel, nbins, xmin, xmax = details
                full_hname = "h_" + histName + cut

                print("\n[DEBUG] Combined Run2 inputs:")
                print(f"[DEBUG] Total background datasets = {len(all_year_bgData)}")
                for d in all_year_bgData:
                    print(f"   BKG: {d.fileName}")

                print(f"[DEBUG] Total signal datasets = {len(all_year_sgData)}")
                for d in all_year_sgData:
                    print(f"   SIG: {d.fileName}")
                try:
                    plot_paper_stack(
                        ([], all_year_bgData, all_year_sgData),
                        histoName=full_hname,
                        totalBin=nbins,
                        outputPath=run2_mc_sig_out,
                        xTitle=xlabel,
                        yTitle="Events",
                        isLogY=True,
                        xmin=xmin,
                        xmax=xmax,
                        year="Run2",
                        do_ratio=False,
                        signal_scale=opts.signal_scale,
                        signal_norm=True,
                        combine_signal_years=True,
                    )
                except Exception as e:
                    print(f"[WARN] Failed combined-years MC+signal-only on {full_hname}: {e}")

    print("\nAll done.")


if __name__ == "__main__":
    main()