from importlib import import_module
import numpy as np

__config = __file__.split('/')[-2]
__config_n_svj_tag = f"configs.{__config}.n_svj_tag"
get_number_of_tagged_svjs = import_module(__config_n_svj_tag).get_number_of_tagged_svjs


make_inference_only = False
with_resonance = False
scale_lumi = 1  # Lumi already taken into account when preparing files


#plotting options for abcd analysis
save_2d_plots = False
save_corr_impacts_plots = 1

save_closure_test_plots =  False
save_closure_test_plots_per_signal_hypo = False            #can be run only if save_closure_test_plots = True 

save_ks_test_maps = False

make_optimization_abcd_grid = False 
make_optimization_abcd_grid_binned = 0
make_optimization_abcd_bayesian = False
make_optimization_abcd_bayesian_binned = False

produce_tables_post_optimization = 0
save_histograms_optimized_cuts = 0
save_histograms_variations = False


# categories = {
#     "inclusive": lambda df: df,
#     "0SVJ_WNAE45": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 45) == 0] for branch in df.keys()},
#     "1SVJ_WNAE45": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 45) == 1] for branch in df.keys()},
#     "2SVJ_WNAE45": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 45) == 2] for branch in df.keys()},
#     "3PSVJ_WNAE45": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 45) >= 3] for branch in df.keys()},
    
#     "inclusive_DNN67": lambda df: df,
#     "0SVJ_WNAE45_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 45) == 0] for branch in df.keys()},
#     "1SVJ_WNAE45_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 45) == 1] for branch in df.keys()},
#     "2PSVJ_WNAE45_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 45) >= 2] for branch in df.keys()},
#     "2SVJ_WNAE45_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 45) == 2] for branch in df.keys()},
#     "3PSVJ_WNAE45_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 45) >= 3] for branch in df.keys()},
    
#     "0SVJ_WNAE30_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 30) == 0] for branch in df.keys()},
#     "1SVJ_WNAE30_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 30) == 1] for branch in df.keys()},
#     "2PSVJ_WNAE30_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 30) >= 2] for branch in df.keys()},
#     "2SVJ_WNAE30_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 30) == 2] for branch in df.keys()},
#     "3PSVJ_WNAE30_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 30) >= 3] for branch in df.keys()},
    
#     "0SVJ_WNAE25_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 25) == 0] for branch in df.keys()},
#     "1SVJ_WNAE25_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 25) == 1] for branch in df.keys()},
#     "2PSVJ_WNAE25_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 25) >= 2] for branch in df.keys()},
#     "2SVJ_WNAE25_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 25) == 2] for branch in df.keys()},
#     "3PSVJ_WNAE25_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 25) >= 3] for branch in df.keys()},
    
#     "0SVJ_WNAE20_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 20) == 0] for branch in df.keys()},
#     "1SVJ_WNAE20_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 20) == 1] for branch in df.keys()},
#     "2PSVJ_WNAE20_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 20) >= 2] for branch in df.keys()},
#     "2SVJ_WNAE20_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 20) == 2] for branch in df.keys()},
#     "3PSVJ_WNAE20_DNN67": lambda df: {branch: df[branch][get_number_of_tagged_svjs(df, 20) >= 3] for branch in df.keys()},
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
    "set_x": 0.67,
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
    "2018_QCD",
#    "TTJets",
#    "WJetsToLNu",
#    "ZJetsToNuNu",
#    "ST",
#    "V+jets",
]
merged_backgrounds = {
    "QCD": ["2016_QCD", "2017_QCD", "2018_QCD"],
    "TTJets": ["2016_TTJets", "2017_TTJets", "2018_TTJets"],
    #"V+jets": ["WJetsToLNu", "ZJetsToNuNu"],
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
        #"cut_values": [0.55, 0.58, 0.62, 0.64, 0.65, 0.66, 0.68],
        #"cut_values": [0.67],
        "cut_values": [0.6, 0.65, 0.7, 0.75, 0.8, 0.85],
        #"cut_values": [0.5, 0.6, 0.7, 0.75, 0.8],
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
        "var_binning": [0.5, 0.55, 0.6, 0.65, 0.7, 0.75],
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

