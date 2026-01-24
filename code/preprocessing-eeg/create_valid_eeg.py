import numpy as np
import pandas as pd
import sys
from glob import glob
import os

session = sys.argv[1]
ndcl_module_path = '/home/data/NDClab/analyses/thrive-theta-ddm/code/python/'
sys.path.append(ndcl_module_path)

from ndclab_py import *

path = f"/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/behavior/{session}/"

if session == "s1_r1":
    exclude_id_list = [3000242, 3000254, 3000255] # these subs are excluded from ALL behavior + ALL EEG
    unusable_soc = [3000047, 3000056, 3000180, 3000182] # these subs are excluded from SOC behavior + SOC EEG + SOC surveys
    eeg_unusable_nonsoc = [3000008] # these subs are excluded from NONSOC EEG
    exclude_eeg = [3000080, 3000161, 3000228] # these subs are excluded from ALL EEG
elif session == "s2_r1":
    exclude_id_list = [] # these subs are excluded from ALL behavior + ALL EEG
    unusable_soc = [] # these subs are excluded from SOC behavior + SOC EEG + SOC surveys
    eeg_unusable_nonsoc = [] # these subs are excluded from NONSOC EEG
    exclude_eeg = [] # these subs are excluded from ALL EEG
elif session == "s3_r1"
    exclude_id_list = [] # these subs are excluded from ALL behavior + ALL EEG
    unusable_soc = [] # these subs are excluded from SOC behavior + SOC EEG + SOC surveys
    eeg_unusable_nonsoc = [] # these subs are excluded from NONSOC EEG
    exclude_eeg = [] # these subs are excluded from ALL EEG

# Find the newest file based on modification time
new_file_path = max(matching_files, key=os.path.getmtime)
behavior_df = pd.read_csv(find_newset_file(f"{path}/*summary*{session}*.csv"))

# subset non-social valid_data
behavior_df_nonsoc = behavior_df[[col for col in behavior_df.columns if ("_nonsoc" in col or "sub" in col)]]
behavior_df_nonsoc = behavior_df_nonsoc[behavior_df_nonsoc["acc_nonsoc"] >= 0.6]

print(f"Full nonsoc-DF length: {behavior_df_nonsoc.shape[0]}")

print(f"Removing subjects {exclude_id_list} from nonsocial condition data")
behavior_df_nonsoc = behavior_df_nonsoc[~behavior_df_nonsoc["sub"].isin(exclude_id_list)].reset_index(drop=True)
print(f"New nonsoc-DF length: {behavior_df_nonsoc.shape[0]} \n")

# now we proceed to criteria-based removal
behavior_df_nonsoc = replace_outliers_with_nan_cols(behavior_df_nonsoc, ["invalid_rt_percent_nonsoc", "skipped_percent_nonsoc"])
behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="invalid_rt_percent_nonsoc")
behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="skipped_percent_nonsoc")

print(f"Removing subjects with unusable EEG: {exclude_eeg}")
behavior_df_nonsoc = behavior_df_nonsoc[~behavior_df_nonsoc["sub"].isin({exclude_eeg)].reset_index(drop=True)

print(f"Removing subjects with unusable nonsoc EEG: {eeg_unusable_nonsoc}")
behavior_df_nonsoc = behavior_df_nonsoc[~behavior_df_nonsoc["sub"].isin({eeg_unusable_nonsoc)].reset_index(drop=True)

print(f"Final non-soc-DF length: {behavior_df_nonsoc.shape[0]} \n")
behavior_df_nonsoc.to_csv(f"{path}/thrive_valid_eeg_nonsoc.csv", index=False)

# subset social valid_data
behavior_df_soc = behavior_df[[col for col in behavior_df.columns if ("_soc" in col or "sub" in col)]]
behavior_df_soc = behavior_df_soc[behavior_df_soc["acc_soc"] >= 0.6]

print(f"Full soc-DF length: {behavior_df_soc.shape[0]}")

print(f"Removing subjects {unusable_soc} from social condition data")
behavior_df_soc = behavior_df_soc[~behavior_df_soc["sub"].isin(unusable_soc)].reset_index(drop=True)

print(f"Removing subjects {exclude_id_list} from social condition data")
behavior_df_soc = behavior_df_soc[~behavior_df_soc["sub"].isin(exclude_id_list)].reset_index(drop=True)
print(f"New soc-DF length: {behavior_df_soc.shape[0]} \n")

# now we proceed to criteria-based removal
behavior_df_soc = replace_outliers_with_nan_cols(behavior_df_soc, ["invalid_rt_percent_soc", "skipped_percent_soc"])
behavior_df_soc = behavior_df_soc.dropna(subset="invalid_rt_percent_soc")
behavior_df_soc = behavior_df_soc.dropna(subset="skipped_percent_soc")

print(f"Removing subjects with unusable EEG: {exclude_eeg}")
behavior_df_soc = behavior_df_soc[~behavior_df_soc["sub"].isin({exclude_eeg)].reset_index(drop=True)

print(f"Final soc-DF length: {behavior_df_soc.shape[0]} \n")
behavior_df_soc.to_csv(f"{path}/thrive_valid_eeg_soc.csv", index=False)

print("Processing has been completed; 2 CSVs were created!") 
