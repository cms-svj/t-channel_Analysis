import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import mplhep as hep
plt.style.use(hep.style.CMS)


# Define output directory
output_dir = 'Nonclosure/VRII-DNN/ControlRegion/'
os.makedirs(output_dir, exist_ok=True)  # Ensure directory exists


# Loop through each edge pair and generate a separate plot

def plot_ABCD_regions(DNN_inner_edges, DNN_outer_edges, MET_inner_edges, MET_outer_edges, Validation_Region, output_dir):
    """
    Generate and save plots for each pair of DNN inner and outer edges, visualizing ABCD regions.
    
    Parameters:
        DNN_inner_edges (list or np.array): Inner edge values.
        DNN_outer_edges (list or np.array): Outer edge values.
        MET_inner_edges (list or np.array): Inner MET edge values.
        MET_outer_edges (list or np.array): Outer MET edge values.
        Validation_Region (str): Label to show on the plot (e.g. 'VRI', 'VRII').
        output_dir (str): Directory to save the plots.
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for DNN_inner_edge, DNN_outer_edge, MET_inner_edge, MET_outer_edge in zip(DNN_inner_edges, DNN_outer_edges, MET_inner_edges, MET_outer_edges):
        regions = [
            ("A", MET_inner_edge, MET_outer_edge, DNN_inner_edge, DNN_outer_edge, "red"),
            ("B", MET_inner_edge, MET_outer_edge, 0, DNN_inner_edge, "#FFFF99"),
            ("C", 0, MET_inner_edge, DNN_inner_edge, DNN_outer_edge, "green"),
            ("D", 0, MET_inner_edge, 0, DNN_inner_edge, "#99CCFF")
        ]

        fig, ax = plt.subplots(figsize=(10, 8))

        for label, x1, x2, y1, y2, region_color in regions:
            width = x2 - x1
            height = y2 - y1
            rect = patches.Rectangle((x1, y1), width, height, color=region_color, alpha=0.7,
                                     edgecolor='black', linewidth=1)
            ax.add_patch(rect)

            ax.plot([x1, x2], [y1, y1], "k-", linewidth=1.5)
            ax.plot([x1, x2], [y2, y2], "k-", linewidth=1.5)
            ax.plot([x1, x1], [y1, y2], "k-", linewidth=1.5)
            ax.plot([x2, x2], [y1, y2], "k-", linewidth=1.5)

        hep.cms.label(rlabel="")

        # Add Validation Region text in the center of the plot
        if isinstance(Validation_Region, str):
            ax.text(0.70, 0.70, Validation_Region,fontsize=36, fontweight="bold", color="black",ha="right", va="top", transform=ax.transAxes)
                    #bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.5", alpha=0.6))

        ax.set_xlim(200, 400)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$\it{p}_{T}^{\mathrm{Miss}}$ [GeV]")
        ax.set_ylabel("Event Level DNN Score")

        xticks = list(ax.get_xticks())
        xticks = [tick for tick in xticks if tick < 400] + [400]
        xtick_labels = [str(int(tick)) if tick < 400 else "10000" for tick in xticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtick_labels)


        legend_patches = [
            patches.Patch(color="red", label="A"),
            patches.Patch(color="#FFFF99", label="B"),
            patches.Patch(color="green", label="C"),
            patches.Patch(color="#99CCFF", label="D"),
        ]
        ax.legend(handles=legend_patches, loc="upper right")

        output_path = os.path.join(output_dir, f"Regions_Inner{DNN_inner_edge:.2f}_Outer{MET_outer_edge:.2f}.pdf")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Plot saved: {output_path}")


DNN_outer_edges = np.linspace(1, 1, 3)
DNN_inner_edges = np.linspace(0.85 , 0.85,3)
MET_outer_edges = np.linspace(225,250,3)
MET_inner_edges = np.linspace(205,225,3)
plot_ABCD_regions(DNN_inner_edges,DNN_outer_edges,MET_inner_edges,MET_outer_edges,'VRI', 'SubRegions/VR1')


DNN_outer_edges = np.linspace(0.2, 0.85, 3)
DNN_inner_edges = np.linspace(0.1, 0.6, 3)
MET_outer_edges = np.linspace(1000,1000,3)
MET_inner_edges = np.linspace(250,250,3)
plot_ABCD_regions(DNN_inner_edges,DNN_outer_edges,MET_inner_edges,MET_outer_edges,'VRII', 'SubRegions/VR2')



DNN_outer_edges = np.linspace(0.2, 0.85, 3)
DNN_inner_edges = np.linspace(0.1, 0.6, 3)
MET_outer_edges = np.linspace(225,250,3)
MET_inner_edges = np.linspace(205,225,3)
plot_ABCD_regions(DNN_inner_edges,DNN_outer_edges,MET_inner_edges,MET_outer_edges,'VRIII', 'SubRegions/VR3')

DNN_outer_edges = np.linspace(1, 1, 15)
DNN_inner_edges = np.linspace(0.85, 0.85, 3)
MET_outer_edges = np.linspace(3000,3000,15)
MET_inner_edges = np.linspace(250,250,15)
plot_ABCD_regions(DNN_inner_edges,DNN_outer_edges,MET_inner_edges,MET_outer_edges,'ABCD', 'SubRegions')

'''
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.stats import gaussian_kde

def plot_ABCD_regions_with_smooth_contours(DNN_inner_edges, DNN_outer_edges, MET_inner_edges, MET_outer_edges, 
                                            Validation_Region, output_dir, x_data, y_data):
    """
    Plot ABCD regions overlaid on 2D smooth density contours of given data.

    Parameters:
        DNN_inner_edges (list or np.array): Inner edge values.
        DNN_outer_edges (list or np.array): Outer edge values.
        MET_inner_edges (list or np.array): Inner MET edge values.
        MET_outer_edges (list or np.array): Outer MET edge values.
        Validation_Region (str): Label to show on the plot (e.g. 'VRI', 'VRII').
        output_dir (str): Directory to save the plots.
        x_data (np.array): Data for x-axis (e.g. MET).
        y_data (np.array): Data for y-axis (e.g. DNN score).
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for DNN_inner_edge, DNN_outer_edge, MET_inner_edge, MET_outer_edge in zip(DNN_inner_edges, DNN_outer_edges, MET_inner_edges, MET_outer_edges):
        fig, ax = plt.subplots(figsize=(10, 8))

        # Create grid for contour
        xx, yy = np.mgrid[200:400:100j, 0:1:100j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([x_data, y_data])
        kernel = gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)

        # Plot smooth density contours
        ax.contour(xx, yy, f, levels=6, colors="black", linewidths=1)
        ax.contourf(xx, yy, f, levels=6, cmap="Reds", alpha=0.3)

        # Define ABCD regions
        regions = [
            ("A", MET_inner_edge, MET_outer_edge, DNN_inner_edge, DNN_outer_edge, "red"),
            ("B", MET_inner_edge, MET_outer_edge, 0, DNN_inner_edge, "#FFFF99"),
            ("C", 0, MET_inner_edge, DNN_inner_edge, DNN_outer_edge, "green"),
            ("D", 0, MET_inner_edge, 0, DNN_inner_edge, "#99CCFF")
        ]

        for label, x1, x2, y1, y2, region_color in regions:
            width = x2 - x1
            height = y2 - y1
            rect = patches.Rectangle((x1, y1), width, height, color=region_color, alpha=0.3,
                                     edgecolor='black', linewidth=1)
            ax.add_patch(rect)
            ax.plot([x1, x2], [y1, y1], "k-", linewidth=1.0)
            ax.plot([x1, x2], [y2, y2], "k-", linewidth=1.0)
            ax.plot([x1, x1], [y1, y2], "k-", linewidth=1.0)
            ax.plot([x2, x2], [y1, y2], "k-", linewidth=1.0)

        # Add CMS-style label manually
        ax.text(0.02, 1.02, "CMS Simulation", transform=ax.transAxes, fontsize=16, fontweight='bold', va='bottom')
        
        # Add Validation Region label
        if isinstance(Validation_Region, str):
            ax.text(0.95, 0.92, Validation_Region, fontsize=32, fontweight="bold", color="black",
                    ha="right", va="top", transform=ax.transAxes)

        ax.set_xlim(200, 400)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$\it{p}_{T}^{\mathrm{Miss}}$ [GeV]")
        ax.set_ylabel("Event Level DNN Score")

        xticks = list(ax.get_xticks())
        xticks = [tick for tick in xticks if tick < 400] + [400]
        xtick_labels = [str(int(tick)) if tick < 400 else "10000" for tick in xticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtick_labels)

        legend_patches = [
            patches.Patch(color="red", label="A"),
            patches.Patch(color="#FFFF99", label="B"),
            patches.Patch(color="green", label="C"),
            patches.Patch(color="#99CCFF", label="D"),
        ]
        ax.legend(handles=legend_patches, loc="upper right")

        output_path = os.path.join(output_dir, f"Regions_Inner{DNN_inner_edge:.2f}_Outer{MET_outer_edge:.2f}.pdf")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Plot saved: {output_path}")

# Example data
np.random.seed(0)
x_data_example = np.random.normal(370, 15, 2000)  # Fake MET
y_data_example = np.random.normal(0.9, 0.05, 2000)  # Fake DNN scores

# Example binning
DNN_outer_edges = np.linspace(1, 1, 15)
DNN_inner_edges = np.linspace(0.85, 0.85, 3)
MET_outer_edges = np.linspace(3000,3000,15)
MET_inner_edges = np.linspace(250,250,15)

# Call the function
plot_ABCD_regions_with_smooth_contours(DNN_inner_edges, DNN_outer_edges, MET_inner_edges, MET_outer_edges,
'ABCD', 'SubRegions', x_data_example, y_data_example)
'''