#Use parameters scans for plotting together ROCs for different parameters
from importlib import import_module

__config = __file__.split('/')[-2]

__config_abcd_analysis = f"configs.{__config}.abcd_analysis_options"
__opt_abcd_analysis = import_module(__config_abcd_analysis)


param_scan_type=["rInv-mMed"]

parameters_scans = {
    "rInv-mMed" : {
         f"t-channel_mMed-{m_med}_mDark-20_rinv-{str(r_inv).replace('.', 'p')}_alpha-peak_yukawa-1": [r_inv, m_med]
             for r_inv in [0.1, 0.3, 0.5, 0.7]
             for m_med in [1000, 2000, 3000, 4000]
    },
}

labels = {
    "rInv-mMed": ["$r_{inv}$", "$m_{\Phi}$ [GeV]"],
}


backgrounds_list = __opt_abcd_analysis.backgrounds_list
merged_backgrounds = __opt_abcd_analysis.merged_backgrounds

