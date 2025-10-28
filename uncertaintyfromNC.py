import os
import pandas as pd
import json
import numpy as np

# # --- Config ---
# models = [
#     # "AKSHAT_Nonclosure",
#     # "AKSHAT_Model_Nonclosure_wp_p9",
#     # "Longtrained_Florian_model_Nonclosure",
#     # "Longtrained_Florian_model_Nonclosure_wp_p9",
#     "Nonclosure",
# ]
# taggers = ["WNAE", "PNET"]
# years = ["_2016","_2017", "_2018"]
# validation_regions = ["VRI", "VRII", "VRIII"]
# categories = ["0SVJ", "1SVJ", "2SVJ", "2PSVJ", "3PSVJ"]

# # --- Main loop ---
# for model in models:
#     for tagger in taggers:
#         model_results = {}
#         excel_rows = []  # Detailed table rows
#         summary_rows = []  # Summary table rows

#         for year in years:
#             year_results = {}
#             for cat in categories:
#                 vr_results = {}
#                 max_abs_nc = 0
#                 max_abs_err = 0

#                 for vr in validation_regions:
#                     path = f"{model}/{tagger}/{vr}/{year}/ControlRegion/Ratio_Data_{cat}.txt"
#                     if not os.path.exists(path):
#                         print(f"⚠️ Missing: {path}")
#                         continue

#                     df = pd.read_csv(path, sep="\t")

#                     data_nc_cols = [c for c in df.columns if "Data" in c and "NonClosure" in c]
#                     data_err_cols = [c for c in df.columns if "Data" in c and "err" in c]
#                     boundary_cols = [c for c in df.columns if "Boundary" in c]

#                     if not (data_nc_cols and data_err_cols and boundary_cols):
#                         print(f"⚠️ Column mismatch in {path}")
#                         continue

#                     nc_col, err_col, boundary_col = data_nc_cols[0], data_err_cols[0], boundary_cols[0]

#                     df = df.dropna(subset=[nc_col, err_col, boundary_col])
#                     if df.empty:
#                         print(f"⚠️ No valid entries (NaNs only) in {path}")
#                         vr_results[vr] = {
#                             "Abs_Max_Data_NonClosure": None,
#                             "Associated_Data_Error": None
#                         }
#                         continue

#                     idx_max = df[nc_col].abs().idxmax()
#                     max_nc = float(abs(df.loc[idx_max, nc_col]))
#                     max_err = float(df.loc[idx_max, err_col])
#                     boundary = float(df.loc[idx_max, boundary_col])

#                     vr_results[vr] = {
#                         "Abs_Max_Data_NonClosure": max_nc,
#                         "Associated_Data_Error": max_err
#                     }

#                     excel_rows.append({
#                         "Year": year,
#                         "Category": cat,
#                         "ValidationRegion": vr,
#                         "Boundary": boundary,
#                         "Abs_Max_Data_NonClosure": max_nc,
#                         "Associated_Data_Error": max_err
#                     })

#                     if max_nc > max_abs_nc:
#                         max_abs_nc = max_nc
#                         max_abs_err = max_err

#                 year_results[cat] = vr_results

#                 # Add summary row (no boundary)
#                 summary_rows.append({
#                     "Year": year,
#                     "Category": cat,
#                     "Abs_Max_Data_NonClosure": max_abs_nc,
#                     "Associated_Data_Error": max_abs_err
#                 })

#             model_results[year] = year_results

#         # --- Save JSON ---
#         output_dir = f"Final-Uncertainties/{model}"
#         os.makedirs(output_dir, exist_ok=True)
#         json_path = f"{output_dir}/{tagger}_uncertainties.json"
#         with open(json_path, "w") as f:
#             json.dump(model_results, f, indent=4)
#         print(f"✅ Saved JSON: {json_path}")

#         # --- Convert to DataFrames ---
#         df_main = pd.DataFrame(excel_rows)
#         df_summary = pd.DataFrame(summary_rows)

#         # --- Derived columns for both ---
#         for df_name, df in [("Detailed", df_main), ("Summary", df_summary)]:
#             diff_sq = df["Abs_Max_Data_NonClosure"]**2 - df["Associated_Data_Error"]**2

#             diff_sq = diff_sq.replace([np.inf, -np.inf], np.nan).fillna(0)
#             num_neg = (diff_sq < 0).sum()
#             num_small = ((diff_sq >= 0) & (diff_sq < 1e-6)).sum()

#             diff_sq[diff_sq < 0] = 0
#             diff_sq[np.abs(diff_sq) < 1e-6] = 0

#             df["SqrtDiff"] = np.sqrt(diff_sq)
#             df["HalfShift"] = df["Abs_Max_Data_NonClosure"] - 0.5 * df["Associated_Data_Error"]

#             df["SqrtDiff"] = df["SqrtDiff"].fillna(0).round(4)
#             df["HalfShift"] = df["HalfShift"].fillna(0).round(4)

#             print(f"📊 {model}/{tagger}/{df_name}: {num_neg} negatives, {num_small} near-zero entries clipped.")

#         # --- Save Excel ---
#         excel_path = f"{output_dir}/{tagger}_uncertainties.xlsx"
#         with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
#             df_main.to_excel(writer, index=False, sheet_name="Detailed")
#             df_summary.to_excel(writer, index=False, sheet_name="Summary")

#         print(f"📘 Saved Excel: {excel_path}\n")
# import os
# import pandas as pd
# import json
# import numpy as np

# # --- Config ---
models = [
    "Potential_Model_wp90",
    "Current_Model_wp90",
    "AKSHAT_Nonclosure",
    "AKSHAT_Model_Nonclosure_wp_p9",
    "Longtrained_Florian_model_Nonclosure",
    "Longtrained_Florian_model_Nonclosure_wp_p9"
]
taggers = ["WNAE", "PNET"]
years = ["2017", "2018"]
validation_regions = ["VRI", "VRII", "VRIII"]
categories = ["0SVJ", "1SVJ", "2SVJ", "2PSVJ", "3PSVJ"]

# --- Main loop ---
for model in models:
    for tagger in taggers:
        model_results = {}

        for year in years:
            year_results = {}

            for cat in categories:
                max_abs_nc = 0

                for vr in validation_regions:
                    path = f"{model}/{tagger}/{vr}/{year}/ControlRegion/Ratio_Data_{cat}.txt"
                    if not os.path.exists(path):
                        print(f"⚠️ Missing: {path}")
                        continue

                    df = pd.read_csv(path, sep="\t")

                    # Identify relevant columns
                    data_nc_cols = [c for c in df.columns if "Data" in c and "NonClosure" in c]
                    data_err_cols = [c for c in df.columns if "Data" in c and "err" in c]

                    if not (data_nc_cols and data_err_cols):
                        print(f"⚠️ Column mismatch in {path}")
                        continue

                    nc_col, err_col = data_nc_cols[0], data_err_cols[0]

                    # Drop NaNs
                    df = df.dropna(subset=[nc_col, err_col])
                    if df.empty:
                        print(f"⚠️ No valid entries (NaNs only) in {path}")
                        continue

                    # Find absolute maximum NC
                    idx_max = df[nc_col].abs().idxmax()
                    max_nc = float(abs(df.loc[idx_max, nc_col]))
                    max_err = float(df.loc[idx_max, err_col])

                    # Apply capping condition
                    if max_nc > 0.4 or max_err > 0.3:
                        max_nc = 0.5

                    if max_nc > max_abs_nc:
                        max_abs_nc = max_nc

                # Store only max per category per year
                year_results[cat] = max_abs_nc

            model_results[year] = year_results

        # --- Save JSON ---
        output_dir = f"Final-Uncertainties/{model}"
        os.makedirs(output_dir, exist_ok=True)
        json_path = f"{output_dir}/{tagger}_uncertainties.json"
        with open(json_path, "w") as f:
            json.dump(model_results, f, indent=4)
        print(f"✅ Saved JSON: {json_path}")
