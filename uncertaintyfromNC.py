import os
import pandas as pd
import json

# --- Config ---
models = ["AKSHAT_Nonclosure", "AKSHAT_Model_Nonclosure_wp_p9", "Longtrained_Florian_model_Nonclosure"]
taggers = ["WNAE", "PNET"]
years = ["2017", "2018"]
validation_regions = ["VRI", "VRII", "VRIII"]
categories = ["0SVJ", "1SVJ", "2PSVJ", "3PSVJ"]

# --- Main loop ---
for model in models:
    for tagger in taggers:
        model_results = {}
        excel_rows = []  # rows for detailed table
        summary_rows = []  # rows for final summary table

        for year in years:
            year_results = {}
            for cat in categories:
                vr_results = {}
                max_abs_nc = 0
                max_abs_err = 0

                for vr in validation_regions:
                    # Path to each Ratio_Data file
                    path = f"{model}/{tagger}/{vr}/{year}/ControlRegion/Ratio_Data_{cat}.txt"
                    if not os.path.exists(path):
                        print(f"⚠️ Missing: {path}")
                        continue

                    # Read tab-separated text file
                    df = pd.read_csv(path, sep="\t")

                    # --- Make sure to use Data columns ---
                    data_nc_cols = [c for c in df.columns if "Data" in c and "NonClosure" in c]
                    data_err_cols = [c for c in df.columns if "Data" in c and "err" in c]
                    boundary_cols = [c for c in df.columns if "Boundary" in c]

                    if not (data_nc_cols and data_err_cols and boundary_cols):
                        print(f"⚠️ Column mismatch in {path}")
                        continue

                    nc_col = data_nc_cols[0]
                    err_col = data_err_cols[0]
                    boundary_col = boundary_cols[0]

                    # Find max absolute non-closure for Data
                    abs_vals = df[nc_col].abs()
                    idx_max = abs_vals.idxmax()

                    max_nc = float(abs(df[nc_col].iloc[idx_max]))
                    max_err = float(df[err_col].iloc[idx_max])
                    boundary = float(df[boundary_col].iloc[idx_max])

                    vr_results[vr] = {
                        "Abs_Max_Data_NonClosure": max_nc,
                        "Associated_Data_Error": max_err,
                        "Boundary": boundary
                    }

                    # Add row for Excel detailed table
                    excel_rows.append({
                        "Year": year,
                        "Category": cat,
                        "ValidationRegion": vr,
                        "Boundary": boundary,
                        "Abs_Max_Data_NonClosure": max_nc,
                        "Associated_Data_Error": max_err
                    })

                    # Track max across VRs for summary
                    if max_nc > max_abs_nc:
                        max_abs_nc = max_nc
                        max_abs_err = max_err

                year_results[cat] = vr_results

                # Add to summary table (one per year/category)
                summary_rows.append({
                    "Category": f"{cat} {year}",
                    "Uncertainty": round(max_abs_nc, 3),
                    "Associated_Error": round(max_abs_err, 3)
                })

            model_results[year] = year_results

        # --- Save JSON ---
        output_dir = f"Uncertainties/{model}"
        os.makedirs(output_dir, exist_ok=True)
        json_path = f"{output_dir}/{tagger}_uncertainties.json"
        with open(json_path, "w") as f:
            json.dump(model_results, f, indent=4)
        print(f"✅ Saved JSON: {json_path}")

        # --- Save Excel ---
        excel_path = f"{output_dir}/{tagger}_uncertainties.xlsx"

        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            # Main table
            df_main = pd.DataFrame(excel_rows)
            df_main.sort_values(by=["Year", "Category", "ValidationRegion"], inplace=True)
            df_main.to_excel(writer, index=False, sheet_name="Detailed")

            # Summary table at bottom
            df_summary = pd.DataFrame(summary_rows)
            df_summary.to_excel(writer, index=False, sheet_name="Summary")

        print(f"📘 Saved Excel: {excel_path}")
