import uproot
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# -----------------------
# Step 1: Open file and load histogram
# -----------------------

path = "2018_allBkg.root"
file = uproot.open(path)

hist2d = file["h_METvsDNN_pre"]

# Extract bin contents and bin edges
values = hist2d.values()
MET_edges = hist2d.axes[0].edges()
DNN_edges = hist2d.axes[1].edges()

# Calculate bin centers
MET_centers = (MET_edges[:-1] + MET_edges[1:]) / 2
DNN_centers = (DNN_edges[:-1] + DNN_edges[1:]) / 2

print(f"Shape of values array: {values.shape}")
print(f"MET edges shape: {MET_edges.shape}")
print(f"DNN edges shape: {DNN_edges.shape}")
print(f"MET centers shape: {MET_centers.shape}")
print(f"DNN centers shape: {DNN_centers.shape}")

# -----------------------
# Step 2: Generate pseudo-events
# -----------------------

# Expand histogram counts into pseudo-event lists
MET_events = []
DNN_events = []

for i in range(len(MET_centers)):
    for j in range(len(DNN_centers)):
        count = int(values[i, j])
        if count > 0:
            MET_events.extend([MET_centers[i]] * count)
            DNN_events.extend([DNN_centers[j]] * count)

MET_events = np.array(MET_events)
DNN_events = np.array(DNN_events)

print(f"Total pseudo-events generated: {len(MET_events)}")

# -----------------------
# Step 3: Define copula independence test
# -----------------------

def copula_independence_test(MET_events, DNN_events, plot_hist=False, n_bins=5):
    N = len(MET_events)
    
    # Step 1: Rank and normalize
    rank_MET = stats.rankdata(MET_events, method='average')
    rank_DNN = stats.rankdata(DNN_events, method='average')
    
    MET_prime = (rank_MET - 0.5) / N
    DNN_prime = (rank_DNN - 0.5) / N
    
    # Step 2: 2D histogram with manual bins
    H, xedges, yedges = np.histogram2d(MET_prime, DNN_prime, bins=[n_bins, n_bins], range=[[0,1],[0,1]])

    if plot_hist:
        plt.figure(figsize=(7,6))
        plt.imshow(H.T, origin='lower', extent=[0,1,0,1], aspect='auto')
        plt.colorbar(label="Counts per bin")
        plt.xlabel("Ranked MET")
        plt.ylabel("Ranked DNN")
        plt.title(f"Empirical Copula Histogram ({n_bins}x{n_bins} bins)")
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.show()

    # Step 3: Chi-squared test
    b_exp = N / (n_bins * n_bins)
    chi2_stat = np.sum((H - b_exp)**2 / b_exp)
    dof = (n_bins * n_bins) - n_bins - n_bins + 1
    p_value = 1 - stats.chi2.cdf(chi2_stat, dof)

    # Step 4: Conclusion
    print(f"Chi2 Statistic: {chi2_stat:.2f}")
    print(f"Degrees of Freedom: {dof}")
    print(f"P-value: {p_value:.6f}")
    
    if p_value > 0.05:
        print("✅ MET and DNN appear statistically independent.")
    else:
        print("❌ MET and DNN appear statistically dependent.")

# -----------------------
# Step 4: Usage
# -----------------------

copula_independence_test(MET_events, DNN_events, n_bins=5)

