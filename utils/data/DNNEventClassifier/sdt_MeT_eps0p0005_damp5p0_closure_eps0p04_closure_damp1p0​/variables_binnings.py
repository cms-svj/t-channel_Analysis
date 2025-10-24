#binnings variables
import numpy as np

max_met = 2500

binning_vars_for_2d = {
    "MET": np.arange(0., max_met, max_met/50) ,
    #"MET": [0, 50, 100, 150, 250, 400, 1500],
    #"MET": [0, 50, 100, 150, 250, 400, 600, 1000, 1500, 2500]
    "dnn_score1": np.arange(0., 1.05, 1.05/21)
}


binning_vars_for_1d = {
    #"MET":  np.arange(0., max_met, max_met/20),
    #"MET": [0, 50, 100, 150, 250, 400, 1500],
    #"MET": [0, 50, 100, 150, 250, 400, 600, 1000, 1500, 2500],
    "MET": [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 1000.0, 1100.0, 1200.0],
    "dnn_score1": np.arange(0., 1.1, 1.1/11)
}

labels = {
    "MET": "MET [GeV]",
    "dnn_score1": "DNN score",
}

