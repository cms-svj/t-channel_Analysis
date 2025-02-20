#Here need to define signal and background samples to process
from importlib import import_module

__config = __file__.split('/')[-2]

__config_samples_to_load = f"configs.{__config}.samples_to_load"
__opt_samples_to_load = import_module(__config_samples_to_load)
backgrounds_list = __opt_samples_to_load.__backgrounds_list

__years = list(__opt_samples_to_load.samples["signals_samples"].keys())
__signal_names = __opt_samples_to_load.samples["signals_samples"][_years[0]].keys()
__background_bins = __opt_samples_to_load.samples["background_samples"][_years[0]].keys()

samples = {
    year: {
        "signals_samples" : {
            signal_name: -1
            for signal_name in __signal_names
        },
        for year in __years
    },
    year: {
        "background_samples": {
            bin_name: -1
            for bin_name in __background_bins
        },
        for year in __years
    },
}
