import uproot
import numpy as np
import glob
import os

# Path containing all .root files
base_dir = "/uscms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis/output/Potential_Model"
root_files = glob.glob(os.path.join(base_dir, "*.root"))

global_min, global_max = np.inf, -np.inf

for path in root_files:
    try:
        with uproot.open(path) as f:
            if "h_dnnEventClassScore_pre" not in f:
                print(f"⚠️ Skipping {path} (no 'h_dnnEventClassScore_pre')")
                continue

            hist = f["h_dnnEventClassScore_pre"]
            values = hist.values()          # bin counts
            edges = hist.axes[0].edges()    # bin edges
            centers = 0.5 * (edges[1:] + edges[:-1])

            # Expand histogram into a 1D array (repeated centers)
            expanded = np.repeat(centers, values.astype(int))

            # Clean up invalids
            expanded = expanded[np.isfinite(expanded)]

            if expanded.size == 0:
                print(f"⚠️ No valid entries in {path}")
                continue

            fmin, fmax = expanded.min(), expanded.max()
            global_min = min(global_min, fmin)
            global_max = max(global_max, fmax)

            print(f"✅ {os.path.basename(path)}: min={fmin:.3f}, max={fmax:.3f}")

    except Exception as e:
        print(f"❌ Failed on {path}: {e}")

print("\n===============================")
print(f"🌍 Global min across files: {global_min:.3f}")
print(f"🌍 Global max across files: {global_max:.3f}")
print("===============================")
