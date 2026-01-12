import numpy as np
import pandas as pd
import sys
from glob import glob
import os


def replace_outliers_with_nan_cols(df, columns_to_check, sd_thresh=3):
    """
    Replaces outliers with NaN.
    Raises KeyError if a column is missing.
    Returns the modified dataframe if successful.
    """
    # 1. Validation: Check all columns exist first
    # If this fails, the error is raised and nothing is returned.
    missing_cols = [col for col in columns_to_check if col not in df.columns]
    if missing_cols:
        raise KeyError(f"Missing columns in dataframe: {missing_cols}")

    # 2. Process columns
    for column in columns_to_check:
        if not pd.api.types.is_numeric_dtype(df[column]):
            print(f"Skipping non-numeric column: {column}")
            continue

        mean = df[column].mean()
        std = df[column].std()
        
        upper = mean + (sd_thresh * std)
        lower = mean - (sd_thresh * std)

        # Vectorized replacement
        outlier_mask = (df[column] > upper) | (df[column] < lower)
        
        # Logging
        n_outliers = outlier_mask.sum()
        if n_outliers > 0:
            print(f"[Log] {column}: Replaced {n_outliers} outliers ({n_outliers/len(df):.2%} of data).")
        else:
            print(f"[Log] {column}: No outliers found.")

        df.loc[outlier_mask, column] = np.nan

    return df

session = sys.argv[1]

path = f"/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/behavior/{session}/"

unusable_soc = [3000047, 3000056, 3000180, 3000182]
exclude_id_list = [3000242, 3000254, 3000255]

matching_files = glob(f"{path}/*summary*{session}*.csv")

# Check if any files were found
if not matching_files:
    print("No matching files found.")
else:
    # Find the newest file based on modification time
    new_file_path = max(matching_files, key=os.path.getmtime)
    behavior_df = pd.read_csv(new_file_path)
    print(f"Full DF length: {behavior_df.shape[0]}")
    print(f"Removing subjects {exclude_id_list}")
    behavior_df = behavior_df[~behavior_df["sub"].isin(exclude_id_list)].reset_index(drop=True)
    print(f"New DF length: {behavior_df.shape[0]}")
    
    # print(f"Found {len(matching_files)} matching files.")
    
    print(f"The newest file is: {new_file_path}")

    behavior_df_nonsoc = behavior_df[[col for col in behavior_df.columns if ("_nonsoc" in col or "sub" in col)]]
    
    # we don't need to judge 6 errors (behavior_df_nonsoc["6_or_more_err_nonsoc"] == 1) for other analysis (e.g, EEG) based on the behavioral data; hence for general subsetting of valid participants, we don't apply 6 error rule here
    behavior_df_nonsoc = behavior_df_nonsoc[behavior_df_nonsoc["acc_nonsoc"] >= 0.6]
    # behavior_df_nonsoc = replace_outliers_with_nan(behavior_df_nonsoc) # this would have been needed for the overall outlier removal, i.e., for the behavioral analysis, but not here;
    behavior_df_nonsoc = replace_outliers_with_nan_cols(behavior_df_nonsoc, ["invalid_rt_percent_nonsoc", "skipped_percent_nonsoc"])
    behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="invalid_rt_percent_nonsoc")
    behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="skipped_percent_nonsoc")
    behavior_df_nonsoc.to_csv(f"{path}/thrive_valid_behav_nonsoc.csv", index=False)
    
    behavior_df_soc = behavior_df[[col for col in behavior_df.columns if ("_soc" in col or "sub" in col)]]
    print(f"Full soc-DF length: {behavior_df_soc.shape[0]}")
    print(f"Removing subjects {unusable_soc} from social condition data")
    behavior_df_soc = behavior_df_soc[~behavior_df_soc["sub"].isin(unusable_soc)].reset_index(drop=True)
    print(f"New soc-DF length: {behavior_df_soc.shape[0]}")
    
    behavior_df_soc = behavior_df_soc[behavior_df_soc["acc_soc"] >= 0.6]
    # behavior_df_soc = replace_outliers_with_nan(behavior_df_soc) # this would have been needed for the overall outlier removal, i.e., for the behavioral analysis, but not here;
    behavior_df_soc = replace_outliers_with_nan_cols(behavior_df_soc, ["invalid_rt_percent_soc", "skipped_percent_soc"])
    behavior_df_soc = behavior_df_soc.dropna(subset="invalid_rt_percent_soc")
    behavior_df_soc = behavior_df_soc.dropna(subset="skipped_percent_soc")
    behavior_df_soc.to_csv(f"{path}/thrive_valid_behav_soc.csv", index=False)
    
    print("Processing has been completed; 2 CSVs were created!") 
