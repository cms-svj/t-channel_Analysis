import uproot
import awkward as ak

source_path = "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims/2016/t_channel_pre_selection/nominal//HTMHT/part-0.root"
target_path = "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//HTMHT/part-0.root"
output_path = "HTMHT_2016_fixed.root"

event_keys = ["RunNum", "LumiBlockNum", "EvtNum"]
electron_branches = [
    "nElectrons",
    "Electrons_/.fPt",
    "Electrons_/.fEta",
    "Electrons_/.fPhi",
    "Electrons_/.fE",
    "Electrons_MTW",
    "Electrons_charge",
    "Electrons_iso",
    "Electrons_mediumID",
    "Electrons_passIso",
    "Electrons_pfRelIso",
    "Electrons_tightID",
    "Electrons_isVeto",
]

# --- Load source ---
print("Loading source...")
with uproot.open(source_path) as src:
    src_tree = src["Events"]
    src_events = src_tree.arrays(event_keys, library="ak")
    src_electrons = src_tree.arrays(electron_branches, library="ak")
print(f"Loaded source: {len(src_events)} events")

# --- Load target ---
print("Loading target...")
with uproot.open(target_path) as tgt:
    tgt_tree = tgt["Events"]
    tgt_events = tgt_tree.arrays(event_keys, library="ak")
    tgt_all = tgt_tree.arrays(library="ak", how=dict)  # load as dictionary of arrays
print(f"Loaded target: {len(tgt_events)} events")

# --- Build hashable tuples ---
src_tuples = [tuple(x) for x in zip(src_events["RunNum"], src_events["LumiBlockNum"], src_events["EvtNum"])]
tgt_tuples = [tuple(x) for x in zip(tgt_events["RunNum"], tgt_events["LumiBlockNum"], tgt_events["EvtNum"])]

# --- Build lookup dict ---
src_lookup = {k: i for i, k in enumerate(src_tuples)}

# --- Copy electrons by event match ---
electron_data = {}
matched_count = 0

for br in electron_branches:
    data = []
    for evt in tgt_tuples:
        idx = src_lookup.get(evt)
        if idx is not None:
            data.append(src_electrons[br][idx])
            matched_count += 1
        else:
            # Fill empty array for vector branches, 0 for scalar
            data.append(ak.Array([]) if "Electrons_/" in br else 0)
    electron_data[br] = ak.Array(data)

print(f"Matched {matched_count} electron branches to target events.")

# --- Merge branches and write output ---
new_tree = {**tgt_all, **electron_data}
with uproot.recreate(output_path) as fout:
    fout["Events"] = new_tree

print(f"✅ Done! Output written to: {output_path}")
