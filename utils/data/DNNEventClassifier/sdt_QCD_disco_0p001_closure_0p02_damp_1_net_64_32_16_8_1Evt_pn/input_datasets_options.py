#Set if original files from pfnano or treemaker are used
from importlib import import_module

__config = __file__.split('/')[-2]

pfnano = False
treemaker = True

#Set lumi 
lumi = 1          #in this case the lumi is not used for correcting the weights, use scale later at optimization step (otherwise choose here lumi in pb, e.g. full Run2 138000 pb )

#Variables from input dataset

__config_files_preparation_options = f"configs.{__config}.files_preparation_options"
__opt_files_preparation_options = import_module(__config_files_preparation_options)
__features_plot = list(__opt_files_preparation_options.binning_input_features.keys())
__features_needed_for_creation = __opt_files_preparation_options.__features_needed_for_creation

__config_constraints_options = f"configs.{__config}.constraints_options"
__opt_constraints_options = import_module(__config_constraints_options)
__features_constraints = __opt_constraints_options.features_to_constrain

__n_jets_max = 4
__index_pairs = [(idx1, idx2) for idx1 in range(__n_jets_max) for idx2 in range(idx1+1, __n_jets_max)]

__features_training_in_file = (
      [f"DeltaEta{idx1}{idx2}GoodJetsAK8" for (idx1, idx2) in __index_pairs]
    + [f"DeltaPhi{idx1}{idx2}GoodJetsAK8" for (idx1, idx2) in __index_pairs]
    + [f"DeltaR{idx1}{idx2}GoodJetsAK8" for (idx1, idx2) in __index_pairs]
    + [f"LundJetPlaneZ{idx1}{idx2}GoodJetsAK8" for (idx1, idx2) in __index_pairs]
)

__other_features = (
      [f"DeltaEta{idx1}{idx2}GoodJetsAK8" for (idx1, idx2) in __index_pairs]
    + [f"DeltaPhi{idx1}{idx2}GoodJetsAK8" for (idx1, idx2) in __index_pairs]
    + [f"DeltaEtaAbs{idx1}{idx2}GoodJetsAK8" for (idx1, idx2) in __index_pairs]
    + [f"DeltaPhiAbs{idx1}{idx2}GoodJetsAK8" for (idx1, idx2) in __index_pairs]
    + [f"DijetMass{idx1}{idx2}GoodJetsAK8" for (idx1, idx2) in __index_pairs]
    + [f"LundJetPlaneZ{idx1}{idx2}GoodJetsAK8" for (idx1, idx2) in __index_pairs]
    + [
        "JetsAK8_/.fPt",
        "JetsAK8_/.fEta",
        "JetsAK8_/.fPhi",
        "JetsAK8_/.fE",
        "JetsAK8_mass",
        "JetsAK8_isGood",
        "JetsAK8_deltaPhiMET",
        "JetsAK8_LundJetPlaneZ",
        "JetsAK8_MTMET",
        "JetsAK8_pNetJetTaggerScore",
        "DeltaPhiMinGoodJetsAK8",
        "HT",
        "ST",
        "MET",
        "METPhi",
    ]
)

variables = list(set(
    __features_training_in_file
    + __features_plot
    + __features_constraints
    + __features_needed_for_creation
    + __other_features
))

