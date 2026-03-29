import os
import shutil
import subprocess

def hadd_and_copy(lostlepton_dir, preselection_dir, output_dir):
    # Create output folder if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    
    # Get all files from both directories
    lostlepton_files = set(os.listdir(lostlepton_dir))
    preselection_files = set(os.listdir(preselection_dir))
    
    # Merge files with matching names
    for file_name in lostlepton_files.union(preselection_files):
        lostlepton_file = os.path.join(lostlepton_dir, file_name)
        preselection_file = os.path.join(preselection_dir, file_name)
        output_file = os.path.join(output_dir, file_name)
        
        # If file exists in both directories, hadd
        if file_name in lostlepton_files and file_name in preselection_files:
            print(f"Merging {file_name} from both directories.")
            hadd_command = f"hadd -f {output_file} {lostlepton_file} {preselection_file}"
            try:
                subprocess.run(hadd_command, shell=True, check=True)
                print(f"Successfully merged into: {output_file}")
            except subprocess.CalledProcessError:
                print(f"Error merging files: {file_name}")
        elif file_name in lostlepton_files:
            # If only in lostlepton, copy to output
            print(f"Copying {file_name} from lostlepton to output directory.")
            shutil.copy(lostlepton_file, output_file)
        elif file_name in preselection_files:
            # If only in preselection, copy to output
            print(f"Copying {file_name} from preselection to output directory.")
            shutil.copy(preselection_file, output_file)

    print("Process completed.")

# User-defined input directories and output directory
lostlepton_dir = "output/t_channel_lost_lepton_control_region_PostHEM_Skims"
preselection_dir = "output/t_channel_pre_selection_PostHEM_Skims"
output_dir = "Pre_ll_output/QCD_cl0p02_net64_skims_PostHEM"  # Change this to your desired output folder name

# Run the function
hadd_and_copy(lostlepton_dir, preselection_dir, output_dir)