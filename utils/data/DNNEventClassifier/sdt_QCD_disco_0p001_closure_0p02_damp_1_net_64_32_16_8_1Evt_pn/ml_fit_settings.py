import ROOT 
#Settings for ML fit

#general settings
year = "2018"
use_grid_optimization = True                   #if True, it will use optimization bayesian results
use_grid_optimization_binned = False           #if True, it will use optimization bayesian results
use_bayesian_optimization = False              #if True, it will use optimization bayesian results
use_bayesian_optimization_binned = False       #if True, it will use optimization bayesian results
observable_name = "MET [GeV]"   #for plotting purposes
signal_region = "B"                            #signal region name
                           #signal region name

#data input settings
use_MC = False                                 #if True, it will use MC for fit (so errors are computed based on weights)
use_data = False                               #if True, it will use data for fit (so errors are computed based on counts)
use_pseudodata = True                          #if True, it will use toys for fit (so errors are computed based on Poisson stat)


#actions
make_toys_study = False                        #if True, it will make toys study
test_signals_injection = False                 #if True, it will test signals injection
split_backgrounds = False                      #if True, it will split backgrounds in different categories
propagate_total_transfer_factor = False        #if True, it will propagate total transfer factor from the total background fit to the single background fits
run_variations = False                          #if True, it will run fit with variations of the weights
var_perc = 0.5                                 #percentage of variation
run_on_morphing = False                        #if True, it will run fit on morphing


#pseudo_data options
pseudo_data_settings = {
    "n_pseudo_data": 1, 
    "pseudo_data_seed": 1234,
    "bkg_only": True,                  
}



backgrounds_to_merge = {

    # "background_EWK": {
    #     "backgrounds_to_merge": ["background_wjets", "background_zjets"],
    # },
}


mc_fit_settings = {
    "run_with_gaussian_approx": False,      #if True, it will run fit with gaussian approximation, else it will run the fit with the scaled poisson likelihood
                                            #with the Gaussian Poisson implementation there is no need to put a threshold on the mc events when doing the bin optimization
}

#fit settings
only_norm = True                              #mode: predicts only normalization from fit
run_blind = True                               #Does not use signal region for fit
verbose = True                                 #Dump workspace
save_workspace = True                          #Save workspace

save_predictions = True                        #Save predictions ABCD method
save_norm = True                               #Save normalization
save_shapes = False                            #Save shapes
save_bare_histos = True                        #Save bare histograms



#bin optimization settings
bins_optimization = True                              #if True, it will run bin optimization
save_binning = True                            #if True, it will save binning in .pkl file (for example to be used in validation regions)
bins_optimization_settings = {
    "bins_count_threshold_weighted": 10,              #minimum number of events per bin 
    "bins_count_threshold_mc": 10,                     #minimum fraction of events per bin - set 10 for Gaussian approx 
    "min_bins": 2,                                    #minimum number of bins
    "deep_search": True,                              #if True, it will try to run bin optimization lowering the bins_threshold
    "verbose": False,                                 #if True, it will print out the bin optimization results
} 


#Signals to inject
signals = ["signals_sChannel_leptons_mZprime-3000GeV_mDark-8GeV_rInv-0p3",
           "signals_sChannel_leptons_mZprime-3000GeV_mDark-8GeV_rInv-0p5",
           "signals_sChannel_leptons_mZprime-3000GeV_mDark-8GeV_rInv-0p7",
           "signals_sChannel_leptons_mZprime-5000GeV_mDark-8GeV_rInv-0p3",
           "signals_sChannel_leptons_mZprime-5000GeV_mDark-8GeV_rInv-0p5",
           "signals_sChannel_leptons_mZprime-5000GeV_mDark-8GeV_rInv-0p7"]

#bakground for splitting
# backgrounds = ["background_qcd", "background_ttbar", "background_wjets", "background_zjets"]
backgrounds = ["background_QCD"]

#stack order for plotting
# stack_order = ["background_EWK", "background_ttbar","background_qcd"]
stack_order = ["background_QCD"]

#backgrounds colors for plotting
backgrounds_colors = {
    "background_QCD": ROOT.kGreen+1,
    # "background_ttbar": ROOT.kOrange-3,
    # "background_EWK": ROOT.kMagenta,
}
