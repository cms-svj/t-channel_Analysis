#Samples used for the training and testing of the ABCD Disco algorithm
from importlib import import_module

__config = __file__.split('/')[-2]

__config_samples_nfiles = f"configs.{__config}.samples_nfiles"
__opt_samples_nfiles = import_module(__config_samples_nfiles)

__signal_names = __opt_samples_nfiles.samples["signals_samples"].keys()
__backgrounds_list = __opt_samples_nfiles.backgrounds_list

__config_general = f"configs.{__config}.abcdDisco_processing"
__opt_general = import_module(__config_general)
__output_generated_datasets = __opt_general.output_generated_datasets

input_files_backgrounds = {
    bkg: f"{__output_generated_datasets}/backgrounds/{bkg}/background_dataset.root"
    for bkg in __backgrounds_list
}

input_files_signals = {
    signal: f"{__output_generated_datasets}/signals/{signal}/{signal}_dataset.root"
    for signal in __signal_names
}

#Set labels for backgrounds and signals
labels_backgrounds_hypo = {}
for idx, background_name in enumerate(__backgrounds_list):
    labels_backgrounds_hypo[background_name] = 1+idx,

labels_signals_hypo = {}
for idx, signal_name in enumerate(__signal_names):
    labels_signals_hypo[signal_name] = -(1+idx)
