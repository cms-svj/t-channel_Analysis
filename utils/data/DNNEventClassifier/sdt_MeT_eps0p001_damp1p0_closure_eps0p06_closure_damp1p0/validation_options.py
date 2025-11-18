#Set the ABCD region
opt_card = "bayesian_optimization"

#Validation regions to fit (here assumes ABCD plane A-B|C-D)
make_histograms_for_validation = False            #if True it will make histograms for validation regions
automatic_validation_regions = True              #if True it will build 3 validation regions automatically based on SR in card [always us this option]
validation_regions = ["C", "AC", "CD"]           #if want to ru  the fit only in specific regions (allowed ones)
steps = 10                                       #number of steps to do for each direction of the validation region
number_of_events_th = 10                         #number of events threshold to consider 

limits_abcd = {
    "x": [0.0 , 1.0],
    "y": [0.0 , 1.0],
}

abcd_regions_labels = {
    "x": "DNN1",
    "y": "DNN2",
}


bins_optimization = True                       #if True, it will run bin optimization
bins_optimization_settings = {
    "bins_count_threshold": 10,                      #minimum number of events per bin
    "min_bins": 2,                             #minimum number of bins
    "deep_search": True,                       #if True, it will try to run bin optimization lowering the bins_threshold
    "verbose": False,                            #if True, it will print out the bin edges
} 

#Settings for ML fit in the validation regions
run_ML_fits =  True                             #if True, it will run ML fits in the validation regions
analyse_results = True                         #Produce .pkl file for plotting
plot_results = True

year = "2018"
run_blind = True                                #Does not use signal region for fit
only_norm = False                                #mode: predicts only normalization from fit
verbose = False                                 #Dump workspace


