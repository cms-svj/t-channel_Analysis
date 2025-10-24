from importlib import import_module
import numpy as np

__config = __file__.split('/')[-2]
__config_n_svj_tag = f"configs.{__config}.n_svj_tag"
get_number_of_tagged_svjs = import_module(__config_n_svj_tag).get_number_of_pn_tagged_svjs

__config_general = f"configs.{__config}.abcdDisco_processing"
__opt_general = import_module(__config_general)
__year = __opt_general.year


make_inference_only = False
with_resonance = False
scale_lumi = 1  # Lumi already taken into account when preparing files


#plotting options for abcd analysis
save_2d_plots = False
save_corr_impacts_plots = True

save_closure_test_plots =  False
save_closure_test_plots_per_signal_hypo = False            #can be run only if save_closure_test_plots = True 

save_ks_test_maps = False

make_optimization_abcd_grid = False 
make_optimization_abcd_grid_binned = False
make_optimization_abcd_bayesian = False
make_optimization_abcd_bayesian_binned = False

produce_tables_post_optimization = 0
save_histograms_optimized_cuts = False
save_histograms_variations = False

dnnCut = 0.74
# categories = {
    # "0SVJ_WP80_DNN85"  : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.8) == 0] for branch in df.keys()},
    # "1SVJ_WP80_DNN85"  : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.8) == 1] for branch in df.keys()},
    # "2PSVJ_WP80_DNN85" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.8) >= 2] for branch in df.keys()},
    # "2SVJ_WP80_DNN85" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.8) == 2] for branch in df.keys()},
    # "3PSVJ_WP80_DNN85" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.8) >= 3] for branch in df.keys()},
    # "1SVJ_WP80" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.8) == 1] for branch in df.keys()},
    # "2PSVJ_WP80" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.8) >= 2] for branch in df.keys()},
    # "1SVJ_WP85" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.85) == 1] for branch in df.keys()},
    # "2PSVJ_WP85" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.85) >= 2] for branch in df.keys()},
    # "1SVJ_WP75" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.75) == 2] for branch in df.keys()},
    # "2PSVJ_WP75" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.75) >= 2] for branch in df.keys()},
    # "2SVJ_WP75" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.75) == 2] for branch in df.keys()},
    # "3PSVJ_WP75" : lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 0.75) >= 3] for branch in df.keys()},
# }

#optimization settings
optimization_settings = {
    "signal_region": "B",                                   #here the convention is A-B | C-D (for double disco the signal region is B)
    "maximum_non_closure": 20,
    "threshold_for_normalized_signal_contamination": 0.3,   #default value is 1.0
    "ks_test_threshold": 0.2,                               #ks test pvalue threshold, if don't want to use it set to 0.0
    "minimum_bkg_CRs": 0,
    "sigma_window": 1,
#    "min_x": 0.5,
#    "max_x": 1,
#    "min_y": 200,
#    "max_y": 800,
}

binned_optimization_settings = {
    "signal_region": "B",                                   #here the convention is A-B | C-D (for double disco the signal region is B)
    "maximum_pull": 10**6,
    "maximum_non_closure": 10**6,
    "threshold_for_normalized_signal_contamination": 10000.,   #default value is 1.0
    "chi_square_test_pvalue": 0.,                               #ks test pvalue threshold, if don't want to use it set to 0.0
    "minimum_bkg_CRs": 0,
    "set_x": dnnCut,
    "set_y": 250,
#    "min_x": 0.69,
#    "max_x": 0.71,
#    "min_y": 249,
#    "max_y": 251,
}

optimization_settings_bayes = {
    "n_workers": 12,                                            #number of parallel jobs
    "initial_points": 50,                                   #number of initial points for bayesian optimization
    "n_iter": 50,
    "signal_region": "B",                                   #here the convention is A-B | C-D (for double disco the signal region is B)
    "maximum_non_closure": 0.35,
    "threshold_for_normalized_signal_contamination": 0.5,   #default value is 1.0
    "ks_test_threshold": 0.50,                               #ks test pvalue threshold, if don't want to use it set to 0.0
    "ks_test_pvalue_n_permutations": 1000,                  #number of permuations to calculate p-value (lowered from default 9999)
    "save_bayes_plots": True,
    "penalized": True,
}

##### Settings for binned optimization with fit ##########################################################################################
#data input settings
signal_region = "B" 
use_MC = 0
use_data = False
use_pseudodata = False
make_asimov = 1

only_norm = True                              #mode: predicts only normalization from fit
run_blind = True                               #Does not use signal region for fit
verbose = False                                 #Dump workspace

mc_fit_settings = {
    "run_with_gaussian_approx": False,      #if True, it will run fit with gaussian approximation, else it will run the fit with the scaled poisson likelihood
                                            #with the Gaussian Poisson implementation there is no need to put a threshold on the mc events when doing the bin optimization
}

#systematics for background
add_systematics = True
syst = {
    "BkgNorm": ("lnN", "bkg", "SR" , 1.05),
}

#bin optimization settings
bins_optimization = True                              #if True, it will run bin optimization
bins_optimization_settings = {
    "bins_count_threshold_weighted": 10,               #minimum number of events per bin 
    "bins_count_threshold_mc": 10,                     #minimum fraction of events per bin - set 10 for Gaussian approx 
    "norm_mode_opt": True,                            #if True, it will optimize the binning only for the transfering region (from where the shape is taken)
    "min_bins": 2,                                    #minimum number of bins
    "deep_search": True,                              #if True, it will try to run bin optimization lowering the bins_threshold
    "verbose": False,                                 #if True, it will print out the bin optimization results
} 
##########################################################################################################################################


#fitting variable settings dictionary
fitting_variable_settings = {
    "fitting_variable": "MET",
    "binning": np.arange(0, 5000.1, (5000)/100),
}

backgrounds_list = [
   f"{__year}_QCD",
   f"{__year}_TTJets",
   f"{__year}_WJetsToLNu",
   f"{__year}_ZJetsToNuNu",
   f"{__year}_ST",
   f"{__year}_V+jets",
]
    


merged_backgrounds = {
    f"{__year}_V+jets": [f"{__year}_WJetsToLNu", f"{__year}_ZJetsToNuNu"],
}

#2D distributions plots settings
variables_2d = ["MET", "dnn_score1"]
colors_palettes_background_hypo = {
    "QCD": "Greens",
    "TTJets": "Oranges",
    "WJetsToLNu": "Blues",
    "ZJetsToNuNu": "Blues",
    "ST": "Purples",
    "V+jets": "Blues",
}


#1D distributions plots settings
variables_1d_corr_impacts_dict = {
    "corr_impacts_dnn_vs_MET": {
        "cut_variable": "dnn_score1",
        # "cut_values": [0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9],
        "cut_values": [0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9],
        # "cut_values": [0.75, 0.77, 0.79, 0.81, 0.83, 0.85, 0.87, 0.89],
        "cut_variable_label": "DNN score",
        "cut_variable_units": "",
        "variable_to_plot": "MET",
    },
}

#closure test plots settings
variables_closure_test_dict = {
    "x": {
        "var_name": "dnn_score1",
        "var_label": "DNN score",
        "var_units": "",
        #"var_binning": np.arange(0.65, 0.94, 0.05),
        "var_binning": [dnnCut],
    },
    "y": {
        "var_name": "MET",
        "var_label": "MET",
        "var_units": "GeV",
        "var_binning": [250],
    }
}


ks_test_and_closure_map_dict = {
    "signal_region": "B",
    "ks_test_pvalue_n_permutations": 1000, 
}

