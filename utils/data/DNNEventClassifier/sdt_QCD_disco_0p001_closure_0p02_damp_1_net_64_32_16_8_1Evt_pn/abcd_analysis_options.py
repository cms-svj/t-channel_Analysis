import numpy as np

#if you don't want to run the optimization, but only apply the results of a previous optimization
make_inference_only = False
#specify paths datasets
target_datasets = {
    "test_nominal": "",
    "test_up": "",
    "test_down": ""
}
target_scores = {
    "test_nominal": "",
    "test_up": "",
    "test_down": ""
}

#path to the optimization card to use
use_optimization_card = ""
with_resonance = False

#plotting options for abcd analysis
save_2d_plots = False
save_corr_impacts_plots = False

save_closure_test_plots =  False
save_closure_test_plots_per_signal_hypo = False            #can be run only if save_closure_test_plots = True 

save_ks_test_maps = False

make_optimization_abcd_grid = True                        #this is the main switch to run the optimization with a grid search
make_optimization_abcd_bayesian = False                    #this is the main switch to run the optimization with a bayesian search
make_optimization_abcd_bayesian_binned =  False           #this is the main switch to run the optimization with a bayesian search and binned ks test
produce_tables_post_optimization = True                   #can be run only if make_optimization_abcd = True
save_histograms_optimized_cuts = True                     #can be run only if make_optimization_abcd = True
save_histograms_variations = False
optimization_to_use = "grid"
weights_variation_percentage = 0.1                                            #variation of the weights of 10% (default value)

scale_lumi = 59692.692 # 2018 full lumi

#optimization settings
optimization_settings = {
    "signal_region": "B",                                   #here the convention is A-B | C-D (for double disco the signal region is B)
    "maximum_non_closure": 0.1,
    "threshold_for_normalized_signal_contamination": 0.3,   #default value is 1.0
    "ks_test_threshold": 0.5,                               #ks test pvalue threshold, if don't want to use it set to 0.0
    "minimum_bkg_CRs": 10,
    "sigma_window": 1.5, 
}

optimization_settings_bayes = {
    "n_workers": 12,                                            #number of parallel jobs
    "initial_points": 50,                                   #number of initial points for bayesian optimization
    "n_iter": 50,
    "signal_region": "B",                                   #here the convention is A-B | C-D (for double disco the signal region is B)
    "maximum_non_closure": 0.1,
    "threshold_for_normalized_signal_contamination": 0.5,   #default value is 1.0
    "ks_test_threshold": 0.50,                               #ks test pvalue threshold, if don't want to use it set to 0.0
    "ks_test_pvalue_n_permutations": 1000,                  #number of permuations to calculate p-value (lowered from default 9999)
    "save_bayes_plots": True,
    "penalized": True,
}

#fitting variable settings dictionary
fitting_variable_settings = {
    "fitting_variable": "MET",
    #"binning":  [1500.0, 1800.0, 2100.0, 2400.0, 2700.0 , 3000.,  3900.0 ,9000.0],
    "binning":  np.arange(0, 5000 , (5000)/100),
}

backgrounds_list = [
    "QCD",
]

merged_backgrounds = {
}


#2D distributions plots settings
variables_2d = ["MET", "dnn_score1"]
colors_palettes_background_hypo = { "QCD" : 'Greens' }


#1D distributions plots settings
variables_1d_corr_impacts_dict = {

    "corr_impacts_dnn_vs_MET": {

    "cut_variable": "dnn_score1",
    "cut_values": [0.1 , 0.2 , 0.5, 0.6 , 0.8],
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
        "var_binning": np.arange(0.0, 1.1 , (1.1)/15),
    },

    "y": {
        "var_name": "MET",
        "var_label": "MET",
        "var_units": "GeV",
        "var_binning": np.arange(0.0, 5000 , (5000)/20),
    }
}


ks_test_and_closure_map_dict = {

    "signal_region": "B",
    "ks_test_pvalue_n_permutations": 1000, 

}

ks_test_and_closure_map_dict = {

    "signal_region": "B",
    "ks_test_pvalue_n_permutations": 1000, 

}
