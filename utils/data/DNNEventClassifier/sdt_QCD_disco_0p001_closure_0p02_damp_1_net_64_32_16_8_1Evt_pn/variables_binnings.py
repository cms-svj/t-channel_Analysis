#binnings variables
import numpy as np

max_met = 2500

binning_vars_for_2d = {
    #"MET": np.arange(0., max_met, max_met/50) ,
    "MET": [0, 50, 100, 150, 250, 400, 1500],
    "dnn_score1": np.arange(0., 1.1, 1.1/20)
}


binning_vars_for_1d = {
    #"MET":  np.arange(0., max_met, max_met/20),
    "MET": [0, 50, 100, 150, 250, 400, 1500],
    #"MET": [0, 50, 100, 150, 250, 400, 600, 1000, 1500, 2500]
    "dnn_score1": np.arange(0., 1.1, 1.1/10)
}

labels = {
    "MET": "MET [GeV]",
    "dnn_score1": "DNN score",
}

