import uproot

# Open the ROOT file
file = uproot.open("root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//HTMHT/part-0.root")

# Access the tree (usually it's named "Events")
tree = file[file.keys()[0]]

# Print each branch name on its own line
for branch in tree.keys():
    print(branch)
