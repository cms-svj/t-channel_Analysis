import numpy as np
import awkward as ak
# Aim: to set the options for the files preparation

prepare_files = True
produce_plots_input_features = True

#general input features binning
binning_input_features = {
    "METPhi": np.linspace(-3.2, 3.2, 41),
    "MET": [50*i for i in range(20)] + [1000 + 100*i for i in range(10)] + [2000, 2500, 3000],
    "GoodJetsAK80_pt": np.linspace(0, 3000, 41),
    "GoodJetsAK81_pt": np.linspace(0, 3000, 41),
    "GoodJetsAK82_pt": np.linspace(0, 2000, 41),
    "GoodJetsAK83_pt": np.linspace(0, 2000, 41),
}


__features_needed_for_creation = [
    "JetsAK8_isGood",
    "JetsAK8_mass",
    "JetsAK8_deltaPhiMET",
    "JetsAK8_LundJetPlaneZ",
    "JetsAK8_MTMET",
    "JetsAK8_pNetJetTaggerScore",
    "DijetMass01GoodJetsAK8",
    "DijetMass02GoodJetsAK8",
    "DijetMass03GoodJetsAK8",
    "DijetMass12GoodJetsAK8",
    "DijetMass13GoodJetsAK8",
    "DijetMass23GoodJetsAK8",
    "METPhi",
    "HT",
]

__features_created = [
    "GoodJetsAK80_deltaPhiMET",
    "GoodJetsAK81_deltaPhiMET",
    "GoodJetsAK82_deltaPhiMET",
    "GoodJetsAK83_deltaPhiMET",
    "GoodJetsAK80_LundJetPlaneZ",
    "GoodJetsAK81_LundJetPlaneZ",
    "GoodJetsAK82_LundJetPlaneZ",
    "GoodJetsAK83_LundJetPlaneZ",
    "GoodJetsAK80_MTMETLog",
    "GoodJetsAK81_MTMETLog",
    "GoodJetsAK82_MTMETLog",
    "GoodJetsAK83_MTMETLog",
    "DijetMass01GoodJetsAK8Log",
    "DijetMass02GoodJetsAK8Log",
    "DijetMass03GoodJetsAK8Log",
    "DijetMass12GoodJetsAK8Log",
    "DijetMass13GoodJetsAK8Log",
    "DijetMass23GoodJetsAK8Log",
]

def add_variables(ak_array):

    import numpy as np
    import awkward as ak
    from coffea.nanoevents.methods import vector
    ak.behavior.update(vector.behavior)

    def pad_array(ak_array, idx, pad_value=0.):
        count_array = ak.count(ak_array, axis=1)
        masked_array = ak.mask(ak_array, count_array > idx)
        variable = masked_array[:, idx]
        variable = ak.fill_none(variable, pad_value)
        return variable
    
    def filter_nans(ak_array):

        filter = None
        for field in ak_array.keys():
            n_dims = str(ak.type(ak_array[field])).count("*")
            if n_dims == 1:
                filter_ = (
                    (~np.isnan(ak_array[field]))
                    & (~np.isinf(ak_array[field]))
                )
                if filter is None:
                    filter = filter_
                else:
                    filter = filter & filter_
        new_ak_array = {}
        for field in ak_array.keys():
            print("field",field)
            print("filter",filter)
            print("ak_array[field]",ak_array[field])
            print()
            print()
            new_ak_array[field] = ak_array[field][filter]
    
        if len(new_ak_array[field]) != len(ak_array[field]):
            fraction = 1 -  (len(new_ak_array) / len(ak_array))
            print(f"NaN / inf cut removed {fraction}% of the events")
    
        return new_ak_array
     

    n_jets = 8

    ak_array[f"JetsAK8_MTMETLog"] = np.log(ak_array[f"JetsAK8_MTMET"])
    variable_names = [
        "deltaPhiMET",
        "LundJetPlaneZ",
        "MTMETLog",
        "pNetJetTaggerScore",
        "/.fPt",
        "/.fEta",
        "/.fPhi",
        "mass",
    ]
    for variable_name in variable_names:
        for idx in range(n_jets):
            if variable_name == "/.fPt":
                target_variable_name = "pt"
            elif variable_name == "/.fEta":
                target_variable_name = "eta"
            elif variable_name == "/.fPhi":
                target_variable_name = "phi"
            else:
                target_variable_name = variable_name
            ak_array[f"GoodJetsAK8_{variable_name}"] = ak_array[f"JetsAK8_{variable_name}"][ak_array["JetsAK8_isGood"]]
            ak_array[f"GoodJetsAK8{idx}_{target_variable_name}"] = pad_array(ak_array[f"GoodJetsAK8_{variable_name}"], idx)
            if variable_name == "deltaPhiMET":
                ak_array[f"GoodJetsAK8{idx}_{variable_name}"] = ak.where(
                ak_array[f"GoodJetsAK8{idx}_{variable_name}"] == 0,
                ak_array["METPhi"],
                ak_array[f"GoodJetsAK8{idx}_{variable_name}"],
            )

    n_jets = 4

    for idx_1 in range(n_jets):
        for idx_2 in range(idx_1+1, n_jets):
            jet_name = "JetsAK8"
            dijet_mass = ak_array[f"DijetMass{idx_1}{idx_2}Good{jet_name}"]
            if jet_name == "JetsAK8":
                jet1_mass = ak_array[f"Good{jet_name}{idx_1}_mass"]
                jet2_mass = ak_array[f"Good{jet_name}{idx_2}_mass"]
                dijet_mass = ak.where(
                    dijet_mass == 0,
                    jet1_mass,
                    dijet_mass,
                )
                dijet_mass = ak.where(
                    dijet_mass == 0,
                    jet2_mass,
                    dijet_mass,
                )
            ak_array[f"DijetMass{idx_1}{idx_2}Good{jet_name}"] = dijet_mass
            ak_array[f"DijetMass{idx_1}{idx_2}Good{jet_name}Log"] = np.log(dijet_mass)
            ak_array[f"DijetMass{idx_1}{idx_2}Good{jet_name}Log"] = ak.nan_to_num(
                array=ak_array[f"DijetMass{idx_1}{idx_2}Good{jet_name}Log"],
                nan=0.,
                posinf=0.,
                neginf=0.,
            )

    # ak_array = filter_nans(ak_array)
    return ak_array


def select_events(ak_array):
    return ak_array

