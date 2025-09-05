import os

#output_directory = f"/work/{os.environ['USER']}/SVJ/abcd_disco/dataset11/training_set1"
output_directory = f"/pscratch/sd/a/ashrivas/SVJ/eventLevelTagger/output/dataset11/training_set1"
year = "2018"

#general inputs/outputs
input_files_tagger = ""
output_generated_datasets = output_directory + "/input_datasets/"
output_prepared_files = output_directory + "/data_inclusive_all_years/"
output_prepared_files_norm = output_prepared_files
output_training = output_directory + "/training_results/"
output_evaluation = output_directory + "/evaluation_results/"
output_analysis = output_directory + "/analysis_results_asimov/"
output_application = output_directory + "/application_results/"

