import awkward as ak
from utils.eventTaggerComputation import event_variables, jet_variables, skimmer_utils
from . import objects as ob
from pydoc import locate
import numpy as np 
import pandas as pd  
import torch

def calculate_skim_variables_from_ntuples(events, varVal):
    # Jets AK8 variables
    new_branches = {}
    obj = ob.Objects(events)
    jets_ak8 = varVal["fjets"]
    jets_ak8_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=jets_ak8.pt,
        eta=jets_ak8.eta,
        phi=jets_ak8.phi,
        mass=jets_ak8.mass,
    )
    met_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.MET,
        phi=events.METPhi,
    )

    # Kinematics
    varVal["JetsAK8_isGood"] = ak.ones_like(jets_ak8.pt,dtype=bool) # we only use good jets here
    varVal["JetsAK8_mass"] = jets_ak8_lv.mass
    varVal["JetsAK8_deltaPhiMET"] = jet_variables.calculate_delta_phi_with_met(jets_ak8_lv, met_lv)
    varVal["JetsAK8_LundJetPlaneZ"] = jet_variables.calculate_lund_jet_plane_z_with_met(jets_ak8_lv, met_lv)
    varVal["JetsAK8_MTMET"] = jet_variables.calculate_invariant_mass_with_met(jets_ak8_lv, met_lv)
    
    # Event variables
    nan_value = 0.  # Natural choice for missing values for LJP variables and delat eta / phi!
    n_jets_max = 4
    for index_0 in range(n_jets_max):
        for index_1 in range(index_0+1, n_jets_max):
            delta_eta = event_variables.calculate_delta_eta(
                physics_objects=jets_ak8_lv,
                indices=(index_0, index_1),
                absolute_value=False,
                nan_value=None,
            )
            delta_phi = event_variables.calculate_delta_phi(
                physics_objects=jets_ak8_lv,
                indices=(index_0, index_1),
                absolute_value=False,
                nan_value=None,
            )
            delta_r = event_variables.calculate_delta_r(
                physics_objects=jets_ak8_lv,
                indices=(index_0, index_1),
                nan_value=nan_value,
            )
            dijet_mass = event_variables.calculate_invariant_mass(
                physics_objects=jets_ak8_lv,
                indices=(index_0, index_1),
                nan_value=nan_value,
            )
            lund_jet_plane_z = event_variables.calculate_lund_jet_plane_z(
                physics_objects=jets_ak8_lv,
                indices=(index_0, index_1),
                nan_value=nan_value,
            )
            delta_eta_abs = abs(delta_eta)
            delta_phi_abs = abs(delta_phi)
            delta_eta = ak.fill_none(delta_eta, nan_value)
            delta_phi = ak.fill_none(delta_eta, nan_value)
            delta_eta_abs = ak.fill_none(delta_eta, nan_value)
            delta_phi_abs = ak.fill_none(delta_eta, nan_value)

            varVal[f"DeltaEta{index_0}{index_1}GoodJetsAK8"] = delta_eta
            varVal[f"DeltaPhi{index_0}{index_1}GoodJetsAK8"] = delta_phi
            varVal[f"DeltaR{index_0}{index_1}GoodJetsAK8"] = delta_r
            varVal[f"DeltaEtaAbs{index_0}{index_1}GoodJetsAK8"] = delta_eta_abs
            varVal[f"DeltaPhiAbs{index_0}{index_1}GoodJetsAK8"] = delta_phi_abs
            varVal[f"DijetMass{index_0}{index_1}GoodJetsAK8"] = dijet_mass
            varVal[f"LundJetPlaneZ{index_0}{index_1}GoodJetsAK8"] = lund_jet_plane_z
   
    varVal["DeltaPhiMinGoodJetsAK8"] = ak.min(abs(varVal["JetsAK8_deltaPhiMET"]), axis=1)
    varVal["ATLASDeltaPhiMinMax"] = event_variables.calculate_atlas_delta_phi_max_min(
        jets=jets_ak8_lv,
        met=met_lv,
        nan_value=nan_value,
    )
    varVal["ATLASPtBalance"] = event_variables.calculate_atlas_momentum_balance(
        jets=jets_ak8_lv,
        met=met_lv,
        nan_value=nan_value,
    )

def get_skim_variables_from_skims(events, varVal):
    # Kinematics
    varVal["JetsAK8_isGood"] = events.JetsAK8.isGood
    varVal["JetsAK8_mass"] = events.JetsAK8.mass
    varVal["JetsAK8_deltaPhiMET"] = events.JetsAK8.deltaPhiMET
    varVal["JetsAK8_LundJetPlaneZ"] = events.JetsAK8.LundJetPlaneZ
    varVal["JetsAK8_MTMET"] = events.JetsAK8.MTMET
    varVal["JetsAK8_pNetJetTaggerScore"] = events.JetsAK8.pNetJetTaggerScore
    # Event variables
    n_jets_max = 4
    for index_0 in range(n_jets_max):
        for index_1 in range(index_0+1, n_jets_max):
            varVal[f"DeltaEta{index_0}{index_1}GoodJetsAK8"] = events[f"DeltaEta{index_0}{index_1}GoodJetsAK8"]
            varVal[f"DeltaPhi{index_0}{index_1}GoodJetsAK8"] = events[f"DeltaPhi{index_0}{index_1}GoodJetsAK8"]
            varVal[f"DeltaR{index_0}{index_1}GoodJetsAK8"] = events[f"DeltaR{index_0}{index_1}GoodJetsAK8"]
            varVal[f"DeltaEtaAbs{index_0}{index_1}GoodJetsAK8"] = events[f"DeltaEtaAbs{index_0}{index_1}GoodJetsAK8"]
            varVal[f"DeltaPhiAbs{index_0}{index_1}GoodJetsAK8"] = events[f"DeltaPhiAbs{index_0}{index_1}GoodJetsAK8"]
            varVal[f"DijetMass{index_0}{index_1}GoodJetsAK8"] = events[f"DijetMass{index_0}{index_1}GoodJetsAK8"]
            varVal[f"LundJetPlaneZ{index_0}{index_1}GoodJetsAK8"] = events[f"LundJetPlaneZ{index_0}{index_1}GoodJetsAK8"]
    varVal["DeltaPhiMinGoodJetsAK8"] = events.DeltaPhiMinGoodJetsAK8
    varVal["ATLASDeltaPhiMinMax"] = events.ATLASDeltaPhiMinMax
    varVal["ATLASPtBalance"] = events.ATLASPtBalance

# add variables used for training the event classifier
def evtTaggerVars(events, varVal, skimSource, evtTaggerDict):
    if skimSource:
        get_skim_variables_from_skims(events, varVal)
    else:
        calculate_skim_variables_from_ntuples(events, varVal)
    # adding variables used in event tagger training
    # config_prepare_files_options = "utils.data.DNNEventClassifier.damp_1_0p001_closure_0p06_net_64_32_16_8_LJP_1EvtABCD_fixedSkimBugs.files_preparation_options"
    config_prepare_files_options = f"{evtTaggerDict['location']}.files_preparation_options"
    opt_prepare_files_options = locate(config_prepare_files_options)
    opt_prepare_files_options.add_variables(varVal)

# normalizing variables (code modeled after ABCD_Disco_framework/evaluate_DNN_ABCD_Disco.py)
def normalize(evtTaggerDict, varVal):
    scaler = evtTaggerDict["scaler"]
    branches_for_norm = scaler.feature_names_in_
    branches_to_keep = branches_for_norm
    df_target = {var: np.array(varVal[var]) for var in branches_for_norm}
    df_target = pd.DataFrame(df_target)
    df_target_std = scaler.transform(df_target)
    #convert to dataframe
    df_target_std = pd.DataFrame(df_target_std, columns=branches_to_keep)
    return df_target_std

# inferencing model
def runEventTagger(events, varVal, skimSource, evtTaggerDict):
    evtTaggerVars(events, varVal, skimSource, evtTaggerDict)
    df_target_std = normalize(evtTaggerDict, varVal)
    model = evtTaggerDict["model"]
    # config_training_options = "utils.data.DNNEventClassifier.sdt_QCD_disco_0p001_closure_0p02_damp_1_net_64_32_16_8_1Evt_pn.training_options"
    config_training_options = f"{evtTaggerDict['location']}.training_options"
    opt_training_options = locate(config_training_options)
    features = opt_training_options.features_training
    event_tagger_score = model(torch.tensor(df_target_std[features].values, dtype=torch.float32))
    varVal["dnnEventClassScore"] = event_tagger_score.detach().numpy().flatten()
