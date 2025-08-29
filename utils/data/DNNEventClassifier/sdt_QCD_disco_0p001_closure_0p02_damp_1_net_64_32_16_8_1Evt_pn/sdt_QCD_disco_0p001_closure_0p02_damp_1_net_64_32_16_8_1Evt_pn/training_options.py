from importlib import import_module

__config = __file__.split('/')[-2]

__config_input_dataset_options = f"configs.{__config}.input_datasets_options"
__opt_input_dataset_options = import_module(__config_input_dataset_options)

__config_files_preparation_options = f"configs.{__config}.files_preparation_options"
__opt_files_preparation_options = import_module(__config_files_preparation_options)

__config_constraints_options = f"configs.{__config}.constraints_options"
__opt_constraints_options = import_module(__config_constraints_options)


# Variables for training
features_training = (
    __opt_input_dataset_options.__features_training_in_file
    + __opt_files_preparation_options.__features_created
)

all_features = list(set(
    __opt_input_dataset_options.variables
    + ["weights"]
))

features_to_add = []
features_to_remove = []


#Choose scaler
scaler = "MinMaxScaler"                    #use when adding closure_loss term
features_to_not_scale = __opt_constraints_options.features_to_constrain

#Training options
double_disco = False                        #use double disco loss
use_init_rnd_weights = True                #use random weights for initialization of the model
epochs = 200
batch_size = 8192
learning_rate_scheduler = None
learning_rate_scheduler_args = {}
learning_rate = 0.001
test_size = 0.2
validation_size = 0.2
random_state_sk = 0
number_of_workers = 0                #0 = no parallelization of processing
optimizer = "Adam"                   #if mddm is used, then optimizer is ignored

#Early stopping
use_early_stopping = False

#Proportions signal and background

#Setting the proportions of signal and background in each batch (assumes bins proportional to x-section)
weight_signal = 0.5                        
weight_background = 0.5
#Setting the proportions of background samples, if set to -1 then proportions are set according to the x-section
backgrounds_relative_proportions = {"QCD": 1.0}

#Set number of threads for training
n_threads = 4

#Set early stopping options
early_stopping_criteria = {
    "BCE_loss": {
        "patience": 100,
        "delta": 1.0,
        "min_loss": 1.1,
    },
    "closure_loss": {
        "patience": 100,
        "delta": 0.1,
        "min_loss": 0.08,
    },
}

