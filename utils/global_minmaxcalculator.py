import uproot
import numpy as np
import glob
import os

# path to your directory
indir = "/uscms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis/output/DNNallbckg_posthyperparam/"
indir = '/uscms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis/output/DNN_allbckg_meteps_p001_cleps_06_nonscaled'
# list all root files
files = glob.glob(os.path.join(indir, "*.root"))
print(files)
mins = []
maxs = []

for fname in files:
    try:
        with uproot.open(fname) as f:
            if "h_dnnEventClassScore_pre" not in f:
                print(f"Skipping {fname}: histogram not found")
                continue

            hist = f["h_dnnEventClassScore_pre"]

            values = hist.values()
            edges = hist.axes[0].edges()
            centers = 0.5 * (edges[1:] + edges[:-1])

            expanded = np.repeat(centers, values.astype(int))

            if len(expanded) == 0:
                print(f"Skipping {fname}: histogram is empty")
                continue

            mins.append(expanded.min())
            maxs.append(expanded.max())

            print(f"{os.path.basename(fname)}: min={expanded.min()}, max={expanded.max()}")

    except Exception as e:
        print(f"Error reading {fname}: {e}")

if mins and maxs:
    print("\nGlobal min across files:", min(mins))
    print("Global max across files:", max(maxs))
else:
    print("No valid histograms found.")
