#Use parameters scans for plotting together ROCs for different parameters
from importlib import import_module

__config = __file__.split('/')[-2]

__config_abcd_analysis = f"configs.{__config}.abcd_analysis_options"
__opt_abcd_analysis = import_module(__config_abcd_analysis)


param_scan_type=["mMed","rinv"]

parameters_scans = {

    "mMed" : [f"t-channel_mMed-{x}_mDark-20_rinv-0p3_alpha-peak_yukawa-1" for x in [1000, 2000, 3000, 4000]],
    "rinv" : [f"t-channel_mMed-2000_mDark-20_rinv-{x}_alpha-peak_yukawa-1" for x in ["0p1", "0p3", "0p5", "0p7"]],

}

backgrounds_list = __opt_abcd_analysis.backgrounds_list
merged_backgrounds = __opt_abcd_analysis.merged_backgrounds

