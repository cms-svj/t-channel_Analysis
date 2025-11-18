#Use parameters scans for plotting together ROCs for different parameters
from importlib import import_module

__config = __file__.split('/')[-2]

__config_abcd_analysis = f"configs.{__config}.abcd_analysis_options"
__opt_abcd_analysis = import_module(__config_abcd_analysis)

__config_general = f"configs.{__config}.abcdDisco_processing"
__opt_general = import_module(__config_general)
__year = __opt_general.year


param_scan_type=["mMed","rinv"]

parameters_scans = {

    "mMed" : [f"{__year}_t-channel_mMed-{x}_mDark-20_rinv-0p3_alpha-peak_yukawa-1" for x in [500, 700, 1000, 1500, 2000, 3000, 4000]],
    "rinv" : [f"{__year}_t-channel_mMed-2000_mDark-20_rinv-{x}_alpha-peak_yukawa-1" for x in ["0p1", "0p3", "0p5", "0p7", "0p9"]],
    #"rinv" : [f"{__year}_t-channel_mMed-2000_mDark-20_rinv-{x}_alpha-peak_yukawa-1" for x in ["0", "0p1", "0p2", "0p3", "0p4", "0p5", "0p6", "0p7", "0p8", "0p9", "1"]],

}

backgrounds_list = __opt_abcd_analysis.backgrounds_list
merged_backgrounds = __opt_abcd_analysis.merged_backgrounds

