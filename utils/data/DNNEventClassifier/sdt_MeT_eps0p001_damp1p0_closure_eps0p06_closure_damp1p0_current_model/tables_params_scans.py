#Use parameters scans for plotting together ROCs for different parameters
from importlib import import_module

__config = __file__.split('/')[-2]

__config_abcd_analysis = f"configs.{__config}.abcd_analysis_options"
__opt_abcd_analysis = import_module(__config_abcd_analysis)

__config_general = f"configs.{__config}.abcdDisco_processing"
__opt_general = import_module(__config_general)
__year = __opt_general.year


param_scan_type = [
    "rInv-mMed",
    "yukawa-mMed",
    "mDark-mMed",
]

parameters_scans = {
    "rInv-mMed" : {
         f"{__year}_t-channel_mMed-{m_med}_mDark-20_rinv-{str(r_inv).replace('.', 'p')}_alpha-peak_yukawa-1": [r_inv, m_med]
             for r_inv in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
             #for r_inv in [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
             for m_med in [500, 600, 700, 800, 1000, 1500, 2000, 3000, 4000]
    },
    "yukawa-mMed" : {
         f"{__year}_t-channel_mMed-{m_med}_mDark-20_rinv-0p3_alpha-peak_yukawa-{str(y).replace('.', 'p')}": [y, m_med]
             for y in [0.1, 0.5, 1, 1.5, 2, 2.5, 2, 3.5]
             for m_med in [500, 600, 700, 800, 1000, 1500, 2000, 3000, 4000]
    },
    "mDark-mMed" : {
         f"{__year}_t-channel_mMed-{m_med}_mDark-{m_dark}_rinv-0p3_alpha-peak_yukawa-1": [m_dark, m_med]
             for m_dark in [1, 5, 10, 20, 30, 40, 50, 60, 80, 100]
             for m_med in [500, 600, 700, 800, 1000, 1500, 2000, 3000, 4000]
    },
}

labels = {
    "rInv-mMed": ["$r_{inv}$", "$m_{\Phi}$ [GeV]"],
    "yukawa-mMed": ["Yukawa coupling $\lambda$", "$m_{\Phi}$ [GeV]"],
    "mDark-mMed": ["$m_{dark}$ [GeV]", "$m_{\Phi}$ [GeV]"],
}


backgrounds_list = __opt_abcd_analysis.backgrounds_list
merged_backgrounds = __opt_abcd_analysis.merged_backgrounds

