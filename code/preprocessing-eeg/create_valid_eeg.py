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

exclude_id_list = [3000242, 3000254, 3000255] # these subs are excluded from ALL behavior + ALL EEG
unusable_soc = [3000047, 3000056, 3000180, 3000182] # these subs are excluded from SOC behavior + SOC EEG + SOC surveys
eeg_unusable_nonsoc = [3000008] # these subs are excluded from NONSOC EEG
exclude_eeg = [3000080, 3000161, 3000228] # these subs are excluded from ALL EEG

matching_files = glob(f"{path}/*summary*{session}*.csv")

# Check if any files were found
if not matching_files:
    print("No matching files found.")
else:
    # Find the newest file based on modification time
    new_file_path = max(matching_files, key=os.path.getmtime)
    behavior_df = pd.read_csv(new_file_path)
    print(f"Full DF length: {behavior_df.shape[0]}")
    print(f"Removing subjects with unusable behavior: {exclude_id_list}")
    behavior_df = behavior_df[~behavior_df["sub"].isin(exclude_id_list)].reset_index(drop=True)
    print(f"New DF length: {behavior_df.shape[0]} \n")
    
    # print(f"Found {len(matching_files)} matching files.")
    
    print(f"The newest file is: {new_file_path}")

    # the idea of removals below is that if subjects shoudl be excluded from both behav and EEG then they are remove before outlier cutoffs; however if outlier subjects excluded only from EEG, they are dropped only in the end

    behavior_df_nonsoc = behavior_df[[col for col in behavior_df.columns if ("_nonsoc" in col or "sub" in col)]]
    print(f"Full nonsoc-DF length: {behavior_df_nonsoc.shape[0]}")    

    behavior_df_nonsoc = behavior_df_nonsoc[behavior_df_nonsoc["acc_nonsoc"] >= 0.6]
    behavior_df_nonsoc = replace_outliers_with_nan_cols(behavior_df_nonsoc, ["invalid_rt_percent_nonsoc", "skipped_percent_nonsoc"])
    behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="invalid_rt_percent_nonsoc")
    behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="skipped_percent_nonsoc")
    print(f"Removing subjects with unusable EEG: {exclude_eeg}") 
    behavior_df_nonsoc = behavior_df_nonsoc[~behavior_df_nonsoc["sub"].isin(exclude_eeg)].reset_index(drop=True)
    print(f"Removing subjects with unusable NONSOC EEG data: {eeg_unusable_nonsoc}")    
    behavior_df_nonsoc = behavior_df_nonsoc[~behavior_df_nonsoc["sub"].isin(eeg_unusable_nonsoc)].reset_index(drop=True)
    print(f"New nonsoc-DF length: {behavior_df_nonsoc.shape[0]} \n")
    behavior_df_nonsoc.to_csv(f"{path}/thrive_valid_eeg_nonsoc.csv", index=False)
    
    behavior_df_soc = behavior_df[[col for col in behavior_df.columns if ("_soc" in col or "sub" in col)]]
    print(f"Full soc-DF length: {behavior_df_soc.shape[0]}")
    print(f"Removing subjects with unusable SOC data: {unusable_soc}")
    behavior_df_soc = behavior_df_soc[~behavior_df_soc["sub"].isin(unusable_soc)].reset_index(drop=True)
    print(f"New soc-DF length: {behavior_df_soc.shape[0]} \n")
    
    behavior_df_soc = behavior_df_soc[behavior_df_soc["acc_soc"] >= 0.6]
    behavior_df_soc = replace_outliers_with_nan_cols(behavior_df_soc, ["invalid_rt_percent_soc", "skipped_percent_soc"])
    behavior_df_soc = behavior_df_soc.dropna(subset="invalid_rt_percent_soc")
    behavior_df_soc = behavior_df_soc.dropna(subset="skipped_percent_soc")
    print(f"Removing subjects with unusable EEG: {exclude_eeg}")
    behavior_df_soc = behavior_df_soc[~behavior_df_soc["sub"].isin(exclude_eeg)].reset_index(drop=True)
    print(f"New soc-DF length: {behavior_df_soc.shape[0]} \n")
    behavior_df_soc.to_csv(f"{path}/thrive_valid_eeg_soc.csv", index=False)
    
    print("Processing has been completed; 2 CSVs were created!") 
