import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import mplhep as hep
%matplotlib inline
%config InlineBackend.figure_format='retina'

# Define output directory
output_dir = 'Nonclosure/VRII-DNN/ControlRegion/'
os.makedirs(output_dir, exist_ok=True)  # Ensure directory exists


# Loop through each edge pair and generate a separate plot

def plot_ABCD_regions(DNN_inner_edges,DNN_outer_edges,MET_inner_edges,MET_outer_edges, output_dir):
    """
    Generate and save plots for each pair of DNN inner and outer edges, visualizing ABCD regions.
    
    Parameters:
        DNN_inner_edges (list or np.array): Inner edge values.
        DNN_outer_edges (list or np.array): Outer edge values.
        met (float): MET value.
        output_dir (str): Directory to save the plots.
    """

    if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    # Loop through each edge pair and generate a separate plot
    for DNN_inner_edge,DNN_outer_edge,MET_inner_edge,MET_outer_edge in zip(DNN_inner_edges,DNN_outer_edges,MET_inner_edges,MET_outer_edges):
        regions = [ ("A", MET_inner_edge, MET_outer_edge,DNN_inner_edge ,DNN_outer_edge,"red"),
                        ("B", MET_inner_edge, MET_outer_edge, 0, DNN_inner_edge,"#FFFF99"),
                        ("C", 0, MET_inner_edge, DNN_inner_edge,DNN_outer_edge,"green"),
                        ("D", 0, MET_inner_edge, 0, DNN_inner_edge,"#99CCFF")
                    ]

        # Create a new plot for this edge pair
        fig, ax = plt.subplots(figsize=(7, 6))

        # Draw each region as a rectangle
        for label, x1, x2, y1, y2, region_color in regions:
            width = x2 - x1
            height = y2 - y1
            rect = patches.Rectangle((x1, y1), width, height, color=region_color, alpha=0.7, edgecolor='black', linewidth=1)
            ax.add_patch(rect)

            #ax.text((x1 + x2) / 2, (y1 + y2) / 2, label, fontsize=14, fontweight="bold", 
            #        color="black", ha="center", va="center", bbox=dict(facecolor="white", alpha=0.5, edgecolor="none"))


            ax.plot([x1, x2], [y1, y1], "k-", linewidth=1.5)  # Bottom border
            ax.plot([x1, x2], [y2, y2], "k-", linewidth=1.5)  # Top border
            ax.plot([x1, x1], [y1, y2], "k-", linewidth=1.5)  # Left border
            ax.plot([x2, x2], [y1, y2], "k-", linewidth=1.5)  # Right border
        # Set axis limits and labels
        hep.cms.label(rlabel="")
        #ax.set_xscale("log")
        ax.set_xlim(200, 300)  # Adjusted MET range
        ax.set_ylim(0, 1)      
        ax.set_xlabel("MET [GeV]")
        ax.set_ylabel("DNN Score")
        #ax.set_title(f'Regions for Edges {inner_edge:.2f} - {outer_edge:.2f}', fontsize=15.5)

        # Create a legend for regions
        legend_patches = [
            patches.Patch(color="red", label="A"),
            patches.Patch(color="#FFFF99", label="B"),
            patches.Patch(color="green", label="C"),
            patches.Patch(color="#99CCFF", label="D"),
        ]
        ax.legend(handles=legend_patches, loc="upper right")

        plt.grid(True, which="both", linestyle="--", linewidth=0.5)

        # Save figure with inner and outer edge values in filename
        output_path = os.path.join(output_dir, f"Regions_Inner{DNN_inner_edge:.2f}_Outer{MET_outer_edge:.2f}.jpg")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()  # Close figure to avoid display issues
        print(f'{MET_inner_edges}')
        print(f"Plot saved: {output_path}")


DNN_outer_edges = np.linspace(1, 1, 10)
DNN_inner_edges = np.linspace(0.6 , 0.6, 10)
MET_outer_edges = np.linspace(225,250,10)
MET_inner_edges = np.linspace(205,225,10)
plot_ABCD_regions(DNN_inner_edges,DNN_outer_edges,MET_inner_edges,MET_outer_edges, 'SubRegions/VR1')


DNN_outer_edges = np.linspace(0.2, 0.6, 10)
DNN_inner_edges = np.linspace(0.1, 0.3, 10)
MET_outer_edges = np.linspace(1000,1000,10)
MET_inner_edges = np.linspace(250,250,10)
plot_ABCD_regions(DNN_inner_edges,DNN_outer_edges,MET_inner_edges,MET_outer_edges, 'SubRegions/VR2')




DNN_outer_edges = np.linspace(0.2, 0.6, 10)
DNN_inner_edges = np.linspace(0.1, 0.3, 10)
MET_outer_edges = np.linspace(225,250,10)
MET_inner_edges = np.linspace(205,225,10)
plot_ABCD_regions(DNN_inner_edges,DNN_outer_edges,MET_inner_edges,MET_outer_edges, 'SubRegions/VR3')
