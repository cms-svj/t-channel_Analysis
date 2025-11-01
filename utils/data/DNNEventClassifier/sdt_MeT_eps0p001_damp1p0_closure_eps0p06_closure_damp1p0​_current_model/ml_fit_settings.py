from importlib import import_module
import ROOT 
#Settings for ML fit

#general settings
__config = __file__.split('/')[-2]
__config_general = f"configs.{__config}.abcdDisco_processing"
__opt_general = import_module(__config_general)
year = __opt_general.year

use_optimization_grid = False
use_optimization_grid_binned = True
use_optimization_bayesian = False
use_optimization_bayesian_binned = False

observable_name = "MET [GeV]"   #for plotting purposes
signal_region = "B"                            #signal region name

#data input settings
use_MC = 0                                 #if True, it will use MC for fit (so errors are computed based on weights)
use_data = False                               #if True, it will use data for fit (so errors are computed based on counts)
use_pseudodata = 0                          #if True, it will use toys for fit (so errors are computed based on Poisson stat)
make_asimov = 1


#actions
make_toys_study = False                        #if True, it will make toys study
test_signals_injection = 0
split_backgrounds = False                      #if True, it will split backgrounds in different categories
propagate_total_transfer_factor = False        #if True, it will propagate total transfer factor from the total background fit to the single background fits
run_variations = False                          #if True, it will run fit with variations of the weights
var_perc = 0.5                                 #percentage of variation
run_on_morphing = False                        #if True, it will run fit on morphing

#categories
categories = [
                # "0SVJ_WP80",
                # "1SVJ_WP80",
                # "2PSVJ_WP80",
                # "2SVJ_WP80",
                # "3PSVJ_WP80",
# "0SVJ_WP80_DNN85",
"1SVJ_WP80_DNN85",
"2PSVJ_WP80_DNN85",
# "2SVJ_WP80_DNN85",
# "3PSVJ_WP80_DNN85",
]


#pseudo_data options
pseudo_data_settings = {
    "n_pseudo_data": 1, 
    "pseudo_data_seed": 1234,
    "bkg_only": True,                  
}


backgrounds_to_merge = {
#    "V+jets": {
#        "backgrounds_to_merge": ["WJetsToLNu", "ZJetsToNuNu"],
#    },
}


mc_fit_settings = {
    "run_with_gaussian_approx": False,      #if True, it will run fit with gaussian approximation, else it will run the fit with the scaled poisson likelihood
                                            #with the Gaussian Poisson implementation there is no need to put a threshold on the mc events when doing the bin optimization
}

#systematics for background
add_systematics = True
per_bin_syst = True
syst = {
    "BkgNorm": ("lnN", "bkg", "SR" , 1.05),
}

#fit settings
only_norm = True                              #mode: predicts only normalization from fit
run_blind = False                             #Does not use signal region for fit
verbose = True                                 #Dump workspace
save_workspace = True                          #Save workspace

save_predictions = True                        #Save predictions ABCD method
save_norm = True                               #Save normalization
save_shapes = True                            #Save shapes
save_bare_histos = True                        #Save bare histograms



#bin optimization settings
fix_binning = False
binning_to_use = [0, 200, 250, 400, 600, 1000, 5000]   #binning to use

bins_optimization = True                              #if True, it will run bin optimization
save_binning = True                                   #if True, it will save binning in .pkl file (for example to be used in validation regions)
bins_optimization_settings = {
    "bins_count_threshold_weighted": 10,               #minimum number of events per bin 
    "bins_count_threshold_mc": 10,                     #minimum fraction of events per bin - set 10 for Gaussian approx 
    "norm_mode_opt": True,                             #if True, it will optimize the binning only for the transfering region (from where the shape is taken)
    "min_bins": 2,                                     #minimum number of bins
    "deep_search": True,                               #if True, it will try to run bin optimization lowering the bins_threshold
    "verbose": False,                                  #if True, it will print out the bin optimization results
} 


#Signals to inject
signals = [
#    "signals_t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
#    "signals_t-channel_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
#    "signals_t-channel_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
#    "signals_t-channel_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
#    "signals_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
#    "signals_t-channel_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
#    "signals_t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

#bakground for splitting
backgrounds = [
    f"{year}_QCD",
    f"{year}_TTJets",
    f"{year}_WJetsToLNu",
    f"{year}_ZJetsToNuNu",
    f"{year}_ST",
]

#stack order for plotting
stack_order = [
    f"{year}_QCD",
    f"{year}_TTJets",
    f"{year}_WJetsToLNu",
    f"{year}_ZJetsToNuNu",
    f"{year}_ST",
]

#backgrounds colors for plotting
backgrounds_colors = {
    f"{year}_QCD": ROOT.TColor.GetColor("#5790fc"),
    f"{year}_TTJets": ROOT.TColor.GetColor("#f89c20"),
    f"{year}_ZJetsToNuNu": ROOT.TColor.GetColor("#e42536"),
    f"{year}_WJetsToLNu": ROOT.TColor.GetColor("#964a8b"),
    f"{year}_ST": ROOT.TColor.GetColor("#9c9ca1"),
    #ROOT.TColor.GetColor("#7a21dd"),
}

