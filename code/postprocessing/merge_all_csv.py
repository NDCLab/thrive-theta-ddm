import os
import pandas as pd
from glob import glob
import numpy as np
from functools import reduce
from datetime import datetime
import re
import sys

ndcl_module_path = '/home/data/NDClab/analyses/thrive-theta-ddm/code/python/'
sys.path.append(ndcl_module_path)

from ndclab_py import *

analysis_path = "/home/data/NDClab/analyses/thrive-theta-ddm/"
dataset_path = "/home/data/NDClab/datasets/thrive-dataset/"
checked_path = f"{dataset_path}/sourcedata/checked/"
derivatives_path = f"{dataset_path}/derivatives/"
iqs_p_path_s1_r1 = find_newest_file(f"{derivatives_path}/preprocessed/redcap/Thriveiqsparent*s1r1*.csv")

ndar_path = find_newest_file(f"/home/data/NDClab/tools/lab-devOps/scripts/ndar_uploads/thrive-dataset/july-2026-nda/all_sessions/ndar_subject01_combined_incomplete.csv")

ndar_data = pd.read_csv(ndar_path, header=1)
ndar_data = ndar_data.pivot(index = "src_subject_id",
    columns="timepoint_label", values=["interview_age"]
)

ndar_data.columns = ndar_data.columns.droplevel(0)
ndar_data["sub"] = ndar_data.index
ndar_data = ndar_data.reset_index(drop=True)
ndar_data = ndar_data[["sub", "s1", "s2", "s3"]]

ndar_data.index.name = None
ndar_data.columns.name = None

# session = "s1_r1"
sessions = ["s1_r1", "s2_r1", "s3_r1"]

for session in sessions:
    if session == "s1_r1":
        session_csv = "s1r1"
        exclude_id_list = [3000242, 3000254, 3000255] # these subs are excluded from ALL behavior + ALL EEG
        unusable_soc = [3000047, 3000056, 3000180, 3000182] # these subs are excluded from SOC behavior + SOC EEG + SOC surveys
        eeg_unusable_nonsoc = [3000008] # these subs are excluded from NONSOC EEG
        exclude_eeg = [3000080, 3000161, 3000228] # these subs are excluded from ALL EEG
        exclude_initstated = [3000066]
        exclude_subset_dyadb = [3000361]
    elif session == "s2_r1":
        session_csv = "s2r1"
    elif session == "s3_r1":
        session_csv = "s3r1"

    csv_output_path = f"{analysis_path}/derivatives/csv/{session}/"

    iqs_p_path = find_newest_file(f"{derivatives_path}/preprocessed/redcap/Thriveiqsparent*{session_csv}*.csv")
    iqs_ch_path = find_newest_file(f"{derivatives_path}/preprocessed/redcap/Thriveiqschild*{session_csv}*.csv")
    bbs_p_path = find_newest_file(f"{derivatives_path}/preprocessed/redcap/Thrivebbsparent*{session_csv}*.csv")
    bbs_ch_path = find_newest_file(f"{derivatives_path}/preprocessed/redcap/Thrivebbschild*{session_csv}*.csv")
    bbs_ra_path = find_newest_file(f"{derivatives_path}/preprocessed/redcap/ThrivebbsRA*{session_csv}*.csv")
    
    behavivor_summary_path = find_newest_file(f"{analysis_path}/derivatives/behavior/{session}/summary*{session}*.csv")
    behav_trial_data_path = find_newest_file(f"{analysis_path}/derivatives/behavior/{session}/full_df*.csv")
    
    ern_path = find_newest_file(f"{analysis_path}/derivatives/csv/{session}/thrive_erp_2*.csv")
    ern_laplacian_path = find_newest_file(f"{analysis_path}/derivatives/csv/{session}/thrive_erp_laplacian_2*.csv")
    tf_path = find_newest_file(f"{analysis_path}/derivatives/csv/{session}/thrive_power_itps*.csv")
    icps_path = find_newest_file(f"{analysis_path}/derivatives/csv/{session}/thrive_icps*.csv")
    
    ddm_path = find_newest_file(f"{analysis_path}/derivatives/behavior/{session}/ddm_fit_{session}*.csv")
    
    # age_data = pd.read_csv(iqs_p_path_s1_r1)
    # Identify columns explicitly to avoid index slicing errors
    # age_cols = [c for c in age_data.columns if "agemos" in c]
    
    # Select only necessary columns
    # age_data = age_data[["record_id"] + age_cols].copy()
    
    # 2. Rename and Transform ID
    # age_data = age_data.rename(columns={"record_id": "sub"})
    # age_data["sub"] = age_data["sub"] - 80000
    
    # 3. Coalesce Age Columns (The "English/Spanish" merge)
    # bfill(axis=1) fills NaNs with the next valid value in the row. 
    # We then take the first column, effectively grabbing the first non-null value found.
    # age_data["age_m"] = age_data[age_cols].bfill(axis=1).iloc[:, 0]
    
    # 4. Drop rows where we couldn't find ANY age
    # age_data = age_data.dropna(subset=["age_m"]).reset_index(drop=True)
    
    # 5. Filter final columns
    # age_data = age_data[["sub", "age_m"]]
    age_data = ndar_data[["sub", session[:2]]]
    age_data = age_data.rename({session[:2]: "age_m"}, axis=1)
    
    # 6. Apply Session Adjustments
    # if session == "s2_r1":
    #     age_data["age_m"] += 9
    # elif session == "s3_r1":
    #     age_data["age_m"] += 18
    
    # 1. Load Data
    sex_data = pd.read_csv(iqs_p_path_s1_r1)
    
    # 2. Rename ID safely
    sex_data = sex_data.rename(columns={"record_id": "sub"})
    sex_data["sub"] = sex_data["sub"] - 80000
    
    # 3. Define columns, this must be hardcoded to s1_r1
    sex_cols = [c for c in sex_data.columns if f"sexbirth_s1_r1_e1" in c]
    
    # 4. Vectorized Merge (No loops)
    # "Take values from English column. If null, fill with Spanish column."
    sex_data["sex"] = sex_data[sex_cols].bfill(axis=1).iloc[:, 0]
    
    # 5. Drop rows where sex is still NaN
    sex_data = sex_data.dropna(subset=["sex"]).reset_index(drop=True)
    
    # 6. Ensure Integer Type (Safety Step)
    # This handles the float '1.0' issue common with NaNs
    sex_data["sex"] = sex_data["sex"].astype(int)
    
    # 7. Final Selection
    sex_data = sex_data[["sub", "sex"]]
    
    # 1. Load Data
    full_behavior = pd.read_csv(behav_trial_data_path)
    
    # 2. Extract Condition from the First Trial
    first_soc = full_behavior.loc[full_behavior["trial_num"] == 1, ["sub", "condition_soc"]].copy(deep=True)
    
    # 3. Rename for clarity
    first_soc = first_soc.rename(columns={"condition_soc": "first_soc"})
    
    # Ensure every subject has exactly one entry for "trial 1". 
    # Duplicates imply data errors (e.g., merged files); Missing implies data loss.
    if first_soc["sub"].duplicated().any():
        print("Warning: Duplicate trial_num=1 found for some subjects. Check data merging.")
        
    first_soc = first_soc.reset_index(drop=True)
    
    # 1. Load Data
    bbs_ra = pd.read_csv(bbs_ra_path)
    
    # 2. Define Columns dynamically
    acid_col = f"bbsratrk_acid_{session}_e1"
    dpid_col = f"bbsratrk_dpid_{session}_e1"
    
    # 3. Filter for in-person kids (Subject ID check)
    # Safety: Ensure we are comparing numbers to numbers
    bbs_ra = bbs_ra[bbs_ra[acid_col] >= 3000000].reset_index(drop=True)
    
    # 4. Rename and Deduplicate
    bbs_ra = bbs_ra.rename(columns={acid_col: "sub"})
    bbs_ra = bbs_ra.drop_duplicates(subset=['sub'], keep='last').reset_index(drop=True)
    
    # 5. Determine Partner Mode (Vectorized & Type-Safe)
    # Convert to string first to safely check prefixes
    dpid_str = bbs_ra[dpid_col].astype(str)
    
    conditions = [
        dpid_str.str.startswith("300"), # Condition 1: In-person (starts with 300)
        dpid_str.str.startswith("10")   # Condition 2: Remote (starts with 10)
    ]
    choices = [1, 0] # 1 = In-person, 0 = Remote
    
    # np.select is much faster and cleaner than a list comprehension with nested ifs
    bbs_ra["dp_inperson"] = np.select(conditions, choices, default=np.nan)
    
    # 6. Final Selection
    dp_mode = bbs_ra[["sub", "dp_inperson"]]
    
    # here all valid IDs for EEG are created, which is based on behavioral data
    behavior_df = pd.read_csv(behavivor_summary_path)
    
    # now, to perform condition-wise outlier removal, we need to subset nonsocial and social conditions and perform removal separately
    
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
    print(f"Final non-soc-DF length: {behavior_df_nonsoc.shape[0]} \n")
    
    valid_behavior_nonsoc = behavior_df_nonsoc["sub"].to_frame()
    
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
    print(f"Final soc-DF length: {behavior_df_soc.shape[0]} \n")
    
    valid_behavior_soc = behavior_df_soc["sub"].to_frame()
    
    # all IDs from the previous step to merge with all available EEG
    id_matrix = valid_behavior_nonsoc.merge(valid_behavior_soc, on="sub", how="outer")
    
    ern_data = pd.read_csv(ern_path)
    ern_laplacian_data = pd.read_csv(ern_laplacian_path)
    tf_data = pd.read_csv(tf_path)
    icps_data = pd.read_csv(icps_path)
    
    data_frames = [
        id_matrix,
        ern_data,
        ern_laplacian_data,
        tf_data,
        icps_data,
    ]
    
    # the merge is left relative to behav_id
    # here all EEG data (ERP + power + ITPS + ICPS) is created
    eeg_data = reduce(lambda left, right: pd.merge(left, right, on="sub", how='left'), data_frames)
    
    print(f"Full EEG DF length: {eeg_data.shape[0]}")
    print(f"Removing subjects with unusable EEG: {exclude_eeg}")
    eeg_data = eeg_data[~eeg_data["sub"].isin(exclude_eeg)].reset_index(drop=True)
    print(f"Removing subjects with unusable behavior: {exclude_id_list}")
    eeg_data = eeg_data[~eeg_data["sub"].isin(exclude_id_list)].reset_index(drop=True)
    print(f"New EEG DF length: {eeg_data.shape[0]} \n")
    
    # EEG data is divided to soc and nonsoc
    eeg_data_soc = eeg_data[[i for i in eeg_data.columns if ("_soc" in i or i == "sub")]]
    eeg_data_soc = eeg_data_soc[eeg_data_soc["sub"].isin(valid_behavior_soc["sub"])].reset_index(drop=True)
    print(f"Full SOC EEG DF length: {eeg_data_soc.shape[0]}")
    print(f"Removing subjects with unusable SOC data: {unusable_soc}")
    eeg_data_soc = eeg_data_soc[~eeg_data_soc["sub"].isin(unusable_soc)].reset_index(drop=True)
    print(f"New SOC EEG DF length: {eeg_data_soc.shape[0]} \n")
    
    eeg_data_nonsoc = eeg_data[[i for i in eeg_data.columns if ("_nonsoc" in i or i == "sub")]]
    eeg_data_nonsoc = eeg_data_nonsoc[eeg_data_nonsoc["sub"].isin(valid_behavior_nonsoc["sub"])].reset_index(drop=True)
    print(f"Full NONSOC EEG DF length: {eeg_data_nonsoc.shape[0]}")
    print(f"Removing subjects with unusable NONSOC EEG data: {eeg_unusable_nonsoc}")
    eeg_data_nonsoc = eeg_data_nonsoc[~eeg_data_nonsoc["sub"].isin(eeg_unusable_nonsoc)].reset_index(drop=True)
    print(f"New NONSOC EEG DF length: {eeg_data_nonsoc.shape[0]} \n")
    
    # the code below merges valid data (based on behavior) with EEG, separately soc and non based on valid IDs
    valid_eeg_soc = valid_behavior_soc.merge(eeg_data_soc, on="sub", how="inner")
    
    print(f"New SOC EEG DF length: {valid_eeg_soc.shape[0]} \n")
    
    valid_eeg_nonsoc = valid_behavior_nonsoc.merge(eeg_data_nonsoc, on="sub", how="inner")
    print(f"New NONSOC EEG DF length: {valid_eeg_nonsoc.shape[0]} \n")
    
    # this code merges SOC and NONSOC EEG so that resulting DF has EEG data for both conditions and is ready for list-wise outlier removal
    merged_valid_eeg_data = valid_eeg_soc.merge(valid_eeg_nonsoc, on="sub", how="outer")
    
    # the code below renames columns to put err and corr in the end of the column name and then creates pairs out of these columns for subsequent difference score computation
    
    tf_columns_original = [i for i in merged_valid_eeg_data.columns if not (
        # exclude ERN/CRN
        "RN_" in i or \
        # exclude sub ID column
        i=="sub"
    )]
    tf_columns_renamed = ["_".join(c.split("_err_")) + "_err" if "_err_" in c else "_".join(c.split("_corr_")) + "_corr" if "_corr_" in c else np.nan for c in tf_columns_original]
    
    for i, orig_c in enumerate(tf_columns_original):
        merged_valid_eeg_data.rename({orig_c: tf_columns_renamed[i]}, axis=1, inplace=True)
    
    columns = tf_columns_renamed
    
    # Function to isolate pairs
    def isolate_pairs(columns):
        pairs = []
        seen = set()
    
        for col in columns:
            parts = col.split('_')
            measure, condition, window, accuracy = parts[0], parts[1], parts[2], parts[-1]
    
            # Create a base identifier without the accuracy part
            base_id = '_'.join(parts[:-1])
    
            if base_id in seen:
                continue
    
            # Find the corresponding pair
            if accuracy == 'err':
                corr_col = f"{base_id}_corr"
            else:
                corr_col = f"{base_id}_err"
    
            if corr_col in columns:
                pairs.append((col, corr_col))
                seen.add(base_id)
    
        return pairs
    
    # Get the pairs
    pairs = isolate_pairs(columns)
    
    # Print the pairs
    # for pair in pairs:
    #     print(pair)
    
    # computation of difference score for all TF measures
    # before computing difference score, outlier removal must be done on the EEG measures (this dataset should only contain EEG)
    merged_valid_eeg_data = replace_outliers_with_nan(merged_valid_eeg_data, sd_thresh=3, exclude_cols=["sub"])
    # Iterate only through the column pairs
    for p in pairs:
        # 1. Define the new column name dynamically
        # Example: 'tf_soc' -> 'tf_soc_diff'
        diff_col_name = "_".join(p[0].split("_")[:-1]) + "_diff"
        
        # 2. Vectorized Subtraction
        # Pandas automatically handles NaNs: 
        # If p[0] is NaN OR p[1] is NaN, the result is automatically NaN.
        merged_valid_eeg_data[diff_col_name] = (
            merged_valid_eeg_data[p[0]] - merged_valid_eeg_data[p[1]]
        )
    
    # computation of difference score for all ERP
    suffixes = ['soc', 'nonsoc', 'soc_laplacian', 'nonsoc_laplacian']
    
    for suffix in suffixes:
        # specific column names
        ern_col = f'ERN_{suffix}'
        crn_col = f'CRN_{suffix}'
        diff_col = f'ERN_min_CRN_{suffix}'
        
        # Vectorized subtraction
        merged_valid_eeg_data[diff_col] = (
            merged_valid_eeg_data[ern_col] - merged_valid_eeg_data[crn_col]
        )
    
    # outlier removal for all difference scores
    merged_valid_eeg_data = replace_outliers_with_nan_cols(merged_valid_eeg_data,
        columns_to_check = [c for c in merged_valid_eeg_data if ("diff" in c or "ERN_min_CRN" in c)])
    
    # computation of collapsed (L+R)/2 scores for ICPS measures
    # note that there is no outlier removal for those measures
    # Identify all 'Left' columns to drive the loop
    left_cols = [c for c in merged_valid_eeg_data.columns if '_L_' in c]
    
    for l_col in left_cols:
        # Construct the corresponding 'Right' column name
        r_col = l_col.replace('_L_', '_R_')
        
        # Check if the Right pair exists
        if r_col in merged_valid_eeg_data.columns:
            # Create new name: 
            # 1. Replace '_L_' with '_' (ICPS_soc_early_DLPFC_L_err -> ICPS_soc_early_DLPFC_err)
            # 2. Append '_collapsed'
            new_col = l_col.replace('_L_', '_') + '_collapsed'
            
            # Calculate mean (NaN + Value = NaN)
            merged_valid_eeg_data[new_col] = (merged_valid_eeg_data[l_col] + merged_valid_eeg_data[r_col]) / 2
    
    # this code creates behavioral data to analyze strictly flanker behavioral measures (hence it doesn't care about unusable EEG-only subjects)
    behavior_df = pd.read_csv(behavivor_summary_path)
    
    # subset nonsocial valid_data
    behavior_df_nonsoc = behavior_df[[col for col in behavior_df.columns if ("_nonsoc" in col or "sub" in col)]]
    
    behavior_df_nonsoc = behavior_df_nonsoc[behavior_df_nonsoc["acc_nonsoc"] >= 0.6]
    behavior_df_nonsoc = behavior_df_nonsoc[behavior_df_nonsoc["6_or_more_err_nonsoc"] == 1]
    
    print(f"Full nonsoc-DF length: {behavior_df_nonsoc.shape[0]}")
    print(f"Removing subjects {exclude_id_list} from nonsocial condition data")
    behavior_df_nonsoc = behavior_df_nonsoc[~behavior_df_nonsoc["sub"].isin(exclude_id_list)].reset_index(drop=True)
    print(f"New nonsoc-DF length: {behavior_df_nonsoc.shape[0]} \n")
    
    # criteria-based removals
    behavior_df_nonsoc = replace_outliers_with_nan_cols(behavior_df_nonsoc, ["invalid_rt_percent_nonsoc", "skipped_percent_nonsoc"])
    behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="invalid_rt_percent_nonsoc")
    behavior_df_nonsoc = behavior_df_nonsoc.dropna(subset="skipped_percent_nonsoc")
    print(f"Final non-soc-DF length: {behavior_df_nonsoc.shape[0]} \n")
    
    # subset only data-related columns
    behavior_df_nonsoc[[i for i in behavior_df_nonsoc if ("sub" in i or\
                                                        ("acc" in i and "con" in i) or\
                                                         ("peri" in i or "pea" in i or "pes" in i))]]
    
    # subset nonsocial valid_data
    behavior_df_soc = behavior_df[[col for col in behavior_df.columns if ("_soc" in col or "sub" in col)]]
    
    behavior_df_soc = behavior_df_soc[behavior_df_soc["acc_soc"] >= 0.6]
    behavior_df_soc = behavior_df_soc[behavior_df_soc["6_or_more_err_soc"] == 1]
    
    print(f"Full soc-DF length: {behavior_df_soc.shape[0]}")
    print(f"Removing subjects {unusable_soc} from social condition data")
    behavior_df_soc = behavior_df_soc[~behavior_df_soc["sub"].isin(unusable_soc)].reset_index(drop=True)
    print(f"Removing subjects {exclude_id_list} from social condition data")
    behavior_df_soc = behavior_df_soc[~behavior_df_soc["sub"].isin(exclude_id_list)].reset_index(drop=True)
    print(f"New soc-DF length: {behavior_df_soc.shape[0]} \n")
    
    # criteria-based removals
    behavior_df_soc = replace_outliers_with_nan_cols(behavior_df_soc, ["invalid_rt_percent_soc", "skipped_percent_soc"])
    behavior_df_soc = behavior_df_soc.dropna(subset="invalid_rt_percent_soc")
    behavior_df_soc = behavior_df_soc.dropna(subset="skipped_percent_soc")
    print(f"Final soc-DF length: {behavior_df_soc.shape[0]} \n")
    
    # subset only data-related columns
    behavior_df_soc[[i for i in behavior_df_soc if ("sub" in i or\
                                                    ("acc" in i and "con" in i) or\
                                                    ("peri" in i or "pea" in i or "pes" in i))]]
    
    merged_valid_behav_data = behavior_df_soc.merge(behavior_df_nonsoc, on="sub", how="outer")
    
    merged_valid_behav_data = replace_outliers_with_nan(merged_valid_behav_data, exclude_cols=["sub"])
    
    id_matrix = pd.Series(list(range(3000000, 3000401)), name = "sub")
    
    iqs_p_data = pd.read_csv(iqs_p_path)
    iqs_p_data = iqs_p_data.rename({"record_id" : "sub"}, axis=1)
    iqs_p_data["sub"] = iqs_p_data["sub"] - 80000
    iqs_p_data = iqs_p_data[
        [i for i in iqs_p_data.columns if (i == "sub" or "scrd" in i)]
    ]
    
    iqs_ch_data = pd.read_csv(iqs_ch_path)
    iqs_ch_data = iqs_ch_data.rename({"record_id" : "sub"}, axis=1)
    iqs_ch_data = iqs_ch_data[
        [i for i in iqs_ch_data.columns if (i == "sub" or "scrd" in i)]
    ]
    
    bbs_p_data = pd.read_csv(bbs_p_path)
    bbs_p_data = bbs_p_data.rename({"record_id" : "sub"}, axis=1)
    bbs_p_data["sub"] = bbs_p_data["sub"] - 80000
    bbs_p_data = bbs_p_data[
        [i for i in bbs_p_data.columns if (i == "sub" or "scrd" in i)]
    ]
    
    bbs_ch_data = pd.read_csv(bbs_ch_path)
    bbs_ch_data = bbs_ch_data.rename({"record_id" : "sub"}, axis=1)
    bbs_ch_data = bbs_ch_data[
        [i for i in bbs_ch_data.columns if (i == "sub" or "scrd" in i)]
    ]
    
    data_frames = [
        id_matrix,
        iqs_p_data,
        iqs_ch_data,
        bbs_p_data,
        bbs_ch_data
    ]
    
    # the merge is left relative to behav_id
    redcap_data = reduce(lambda left, right: pd.merge(left, right, on="sub", how='left'), data_frames)
    # remove IDs that were not present
    redcap_data = redcap_data.dropna(how="all", subset = [c for c in redcap_data if c != "sub"]).reset_index(drop=True)
    redcap_data.shape[0]
    
    # merge spanish columns with corresponding english columns
    
    # 1. Identify Spanish columns
    spanish_columns = [col for col in redcap_data.columns if "es_" in col]
    
    for sp_col in spanish_columns:
        # 2. Derive English column name
        # Your logic: replaces "es_" with "_"
        # Example: "demoes_gender" -> "demo_gender"
        eng_col = "_".join(sp_col.split("es_"))
        
        # Safety Check: Ensure the target English column actually exists
        if eng_col not in redcap_data.columns:
            print(f"[Warning] Could not find English match '{eng_col}' for '{sp_col}'. Skipping.")
            continue
    
        # 3. Vectorized Merge (No loops)
        redcap_data[eng_col] = redcap_data[eng_col].combine_first(redcap_data[sp_col])
    
    # 4. Cleanup
    # Only drop the columns we actually processed
    redcap_data = redcap_data.drop(columns=spanish_columns, errors='ignore')
    
    # transform state survey data to reflect the order of conditions (alone/observed) within a session
    state_surveys = pd.read_csv(bbs_ch_path)
    state_surveys = state_surveys.rename({"record_id": "sub"}, axis=1)
    
    state_surveys = state_surveys[
    [i for i in state_surveys.columns if (i == "sub" or "selfnowa" in i or "initstatec" in i or "posttaske" in i or "dyada" in i or "initstated" in i or "posttaskf" in i or "dyadb" in i)\
     and ("timestamp" not in i) and ("_complete" not in i)]
    ]
    
    state_surveys = state_surveys.merge(first_soc, on="sub", how="left")
    
    new_state_survey_df = pd.DataFrame()
    # state_surveys[[i for i in state_surveys.columns if ("initstatec" in i or "first_soc" in i)]]
    for c, num_items in zip(["initstatec", "posttaske"], [5, 10]):
        for i in range(state_surveys.shape[0]):
            new_state_survey_df.loc[i, "sub"] = state_surveys.loc[i, "sub"]
            for item in range(1, num_items + 1):
                if state_surveys.loc[i, "first_soc"] == 1:
                    new_state_survey_df.loc[i, f"{c}_i{item}_{session}_soc"] = state_surveys.loc[i, f"{c}_i{item}_{session}_e1"]
                    new_state_survey_df.loc[i, f"{c}_i{item}_{session}_nonsoc"] = state_surveys.loc[i, f"{c}_i{item}_{session}_e2"]
                elif state_surveys.loc[i, "first_soc"] == 0:
                    new_state_survey_df.loc[i, f"{c}_i{item}_{session}_soc"] = state_surveys.loc[i, f"{c}_i{item}_{session}_e2"]
                    new_state_survey_df.loc[i, f"{c}_i{item}_{session}_nonsoc"] = state_surveys.loc[i, f"{c}_i{item}_{session}_e1"]
                else:
                    new_state_survey_df.loc[i, f"{c}_i{item}_{session}_soc"] = np.nan
                    new_state_survey_df.loc[i, f"{c}_i{item}_{session}_nonsoc"] = np.nan
    
    new_state_survey_df = new_state_survey_df.dropna(how="all", subset = new_state_survey_df.columns[1:]).reset_index(drop=True)
    state_surveys = new_state_survey_df.merge(state_surveys[[i for i in state_surveys.columns if not ("initstatec" in i or "posttaske" in i or "first_soc" in i)]], on="sub", how="left")
    redcap_data = redcap_data.merge(state_surveys, on="sub", how="outer")
    
    # --- PART 1: Broad Exclusion (unusable_soc) ---
    # Columns to wipe for the "Social" exclusion group
    exclude_patterns = ["initstatec", "posttaske", "dyada", "dyadb", "initstated"]
    cols_broad = [c for c in redcap_data.columns 
                  if np.any([pat in c for pat in exclude_patterns]) 
                  and "_nonsoc" not in c]
    
    # Apply NaN for 'unusable_soc' subjects
    rows_broad = redcap_data["sub"].isin(unusable_soc)
    if cols_broad:
        redcap_data.loc[rows_broad, cols_broad] = np.nan
        print(f"Broad clean: Wiped {len(cols_broad)} columns for subjects {unusable_soc}.")
    
    # --- PART 2: Specific Exclusion (exclude_initstated) ---
    # Identify ONLY columns containing 'initstated' (and not nonsocial, if applicable)
    cols_init = [c for c in redcap_data.columns 
                 if "initstated" in c 
                 and "_nonsoc" not in c]
    
    # Apply NaN for 'exclude_initstated' subjects
    rows_init = redcap_data["sub"].isin(exclude_initstated)
    
    if cols_init:
        redcap_data.loc[rows_init, cols_init] = np.nan
        print(f"Specific clean: Wiped 'initstated' columns for subjects {exclude_initstated}.")
    else:
        print("Warning: No 'initstated' columns found.")
    
    
    cols_init = [c for c in redcap_data.columns 
                 if "initstated" in c 
                 and "_nonsoc" not in c]
    
    # Apply NaN for 'exclude_subset_dyadb' subjects
    rows_dyadb = redcap_data["sub"].isin(exclude_subset_dyadb)
    items_to_exclude = [2, 3, 4, 6, 7, 8]
    
    cols_dyadb = [f"dyadb_i{item}_{session}_e1" for item in items_to_exclude]
    
    if cols_init:
        redcap_data.loc[rows_dyadb, cols_dyadb] = np.nan
        print(f"Specific clean: Wiped 'dyadb' columns for subjects {exclude_subset_dyadb}.")
    else:
        print("Warning: No 'dyadb' columns found.")
    
    redcap_data = replace_outliers_with_nan(redcap_data, exclude_cols=["sub"])
    
    # here all valid IDs for DDM are created, which is based on behavioral data
    behavior_df = pd.read_csv(behavivor_summary_path)
    
    # now, to perform condition-wise outlier removal, we need to subset nonsocial and social conditions and perform removal separately
    
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
    print(f"Final non-soc-DF length: {behavior_df_nonsoc.shape[0]} \n")
    
    valid_behavior_nonsoc = behavior_df_nonsoc["sub"].to_frame()
    
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
    print(f"Final soc-DF length: {behavior_df_soc.shape[0]} \n")
    
    valid_behavior_soc = behavior_df_soc["sub"].to_frame()
    
    ddm_df = pd.read_csv(ddm_path)
    ddm_df = ddm_df.drop("seed", axis=1)
    
    # subset subjects based on valid behavioral flanker data condition-wise
    ddm_df = ddm_df[
        (ddm_df['soc'] == 0) | 
        ((ddm_df['soc'] == 1) & (ddm_df['sub'].isin(valid_behavior_soc["sub"])))
    ].reset_index(drop=True)
    print(f"Full DF length: {ddm_df.shape[0]}")
    
    ddm_df = ddm_df[
        (ddm_df['soc'] == 1) | 
        ((ddm_df['soc'] == 0) & (ddm_df['sub'].isin(valid_behavior_nonsoc["sub"])))
    ].reset_index(drop=True)
    print(f"Full DF length: {ddm_df.shape[0]}")
    
    # subset subjects who have enough post-error trials
    people_passed_cutoff = pd.DataFrame()
    cutoff = 16
    
    full_behavior = pd.read_csv(behav_trial_data_path)
    
    counter = 0
    for i, sub in enumerate(full_behavior["sub"].unique()):
        # print(sub)
        sub_data = full_behavior[full_behavior["sub"] == sub]
        
        for cond in [0, 1]:
            n_posterr_trials = []
            data = sub_data[sub_data["condition_soc"] == cond]
    
            n_trials = data[(data["sub"] == sub) & (data["pre_valid_rt"] == 1) & (data["pre_extra_resp"] == 0)\
            & (data["pre_no_resp"] == 0) & (data["pre_congruent"] == 0) & (data["valid_rt"] == 1) & (data["no_resp"] == 0)\
            & (data["pre_accuracy"] == 0)].shape[0]
    
            if n_trials >= cutoff:
                people_passed_cutoff.loc[counter, "sub"] = sub
                people_passed_cutoff.loc[counter, "soc"] = cond
                people_passed_cutoff.loc[counter, "n_posterr_trial"] = n_trials
                counter += 1
    
    ddm_df = ddm_df.merge(people_passed_cutoff[['sub', 'soc']], on=['sub', 'soc'], how='inner')
    # remove subjects with bad fit
    ddm_df = ddm_df[ddm_df["fitStat"] <= 200].reset_index(drop=True)
    print(f"Full DF length: {ddm_df.shape[0]}")
    
    # convert data to wide format
    
    # ddm_df_soc = ddm_df[ddm_df["soc"] == 1]
    # ddm_df_nonsoc = ddm_df[ddm_df["soc"] == 0]
    
    def get_condition_suffix(row):
        acc_str = 'postcorr' if row['pre_accuracy'] == 1 else 'posterr'
        soc_str = 'soc' if row['soc'] == 1 else 'nonsoc'
        return f"{acc_str}_{soc_str}"
    
    ddm_df['condition'] = ddm_df.apply(get_condition_suffix, axis=1)
    
    df_pivoted = ddm_df.pivot(index='sub', columns='condition')
    
    df_pivoted.columns = [f"{col[0]}_{col[1]}" for col in df_pivoted.columns]
    
    df_pivoted.reset_index(inplace=True)
    
    cols_to_drop = [c for c in df_pivoted.columns if c.startswith('pre_accuracy') or c.startswith('soc_')]
    ddm_df = df_pivoted.drop(columns=cols_to_drop, errors='ignore')
    
    ddm_df = replace_outliers_with_nan(ddm_df, exclude_cols=["sub"])
    
    # create attentional ratio columns
    # Identify all 'Left' columns to drive the loop
    rd_cols = [c for c in ddm_df.columns if 'rd_' in c]
    
    for rd_col in rd_cols:
        # Construct the corresponding 'Right' column name
        sda_col = rd_col.replace('rd_', 'sda_')
        
        # Check if the Right pair exists
        if sda_col in ddm_df.columns:
            # Create new name: 
            # 1. Replace '_L_' with '_' (ICPS_soc_early_DLPFC_L_err -> ICPS_soc_early_DLPFC_err)
            # 2. Append '_collapsed'
            new_col = rd_col.replace('rd_', 'reversed_ratio_')
            
            # Calculate mean (NaN + Value = NaN)
            ddm_df[new_col] = (ddm_df[sda_col] / ddm_df[rd_col]) * (-1)
    
    # compute difference scores
    # Parameters and conditions to iterate over
    params = ['a', 'ter', 'p', 'sda', 'rd', 'reversed_ratio']
    conditions = ['soc', 'nonsoc']
    
    for param in params:
        for cond in conditions:
            # Construct existing column names
            posterr_col = f'{param}_posterr_{cond}'
            postcorr_col = f'{param}_postcorr_{cond}'
            
            # Define new difference column name
            # e.g., a_posterr_min_postcorr_soc
            diff_col = f'{param}_diff_{cond}'
            
            # Vectorized subtraction
            # (posterr - postcorr)
            ddm_df[diff_col] = ddm_df[posterr_col] - ddm_df[postcorr_col]
    
    ddm_df = replace_outliers_with_nan_cols(ddm_df, columns_to_check = [c for c in ddm_df if "_diff_" in c])
    
    id_matrix = pd.Series(list(range(3000000, 3000401)), name = "sub")
    
    data_frames = [
        id_matrix,
        age_data,
        sex_data,
        first_soc,
        dp_mode,
        merged_valid_eeg_data,
        merged_valid_behav_data,
        ddm_df,
        redcap_data
    ]
    
    merged_df = reduce(lambda left, right: pd.merge(left, right, on="sub", how='left'), data_frames)
    print(merged_df.shape)
    
    # remove IDs that were not present
    merged_df = merged_df.dropna(how="all", subset = [c for c in merged_df if c != "sub"]).reset_index(drop=True)
    print(merged_df.shape)
    
    date_time = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
    merged_df.to_csv(f"{csv_output_path}/thrive_wide_{session}_{date_time}.csv", index=False)
    
    import pandas as pd
    
    # Load data
    df = pd.read_csv(find_newest_file(f"{csv_output_path}/thrive_wide_{session}_*.csv"))
    
    # 1. Define Column Categories
    id_cols = ['sub', 'age_m', 'sex', 'first_soc', 'dp_inperson']
    redcap_cols = [c for c in df if session in c]
    
    # Social-Only Variables (Repeated for acc=0 and acc=1)
    # These vary by social condition but are NOT split by accuracy in the wide format
    social_only_bases_input = [
        'n_trials', 'invalid_rt_percent', 'skipped_percent', 'acc', 'acc_con', 
        'acc_incon', 'rt_con', 'rt_incon', 'rt_con_log', 'rt_incon_log', 
        'rt_corr_log', 'rt_err_log', 'pes', 'pea', 'peri_acc', 'peri_rt', 
        '6_or_more_err'
    ]
    
    col_map = {}
    
    for c in df.columns:
        # ID Columns
        if c in id_cols:
            col_map[c] = {'type': 'id'}
            continue
            
        mapped = False
        
        # --- 1. Check Explicit Social-Only List ---
        for base in social_only_bases_input:
            if c == f"{base}_soc":
                # Rename 'acc' to 'accuracy_score' to avoid conflict with index 'acc'
                b = 'accuracy_score' if base == 'acc' else base
                col_map[c] = {'type': 'social', 'base': b, 'soc': 1}
                mapped = True
                break
            elif c == f"{base}_nonsoc":
                b = 'accuracy_score' if base == 'acc' else base
                col_map[c] = {'type': 'social', 'base': b, 'soc': 0}
                mapped = True
                break
        if mapped: continue
        
        # --- 2. Check Social-Only Patterns (diff, ERN_min_CRN) ---
        # Matches 'a_diff_soc', 'power_soc_early_diff', etc.
        if 'diff' in c or 'ERN_min_CRN' in c:
            if '_soc' in c:
                base = c.replace('_soc', '') 
                col_map[c] = {'type': 'social', 'base': base, 'soc': 1}
            elif '_nonsoc' in c:
                base = c.replace('_nonsoc', '')
                col_map[c] = {'type': 'social', 'base': base, 'soc': 0}
            mapped = True
            continue
        
        # --- 3. Check Post-Error / Post-Correct (NEW) ---
        # Handled before generic _err/_corr to ensure correct base extraction
        if '_posterr' in c or '_postcorr' in c:
            is_err = '_posterr' in c
            
            if '_soc' in c:
                soc = 1
                # Remove suffixes to get base (e.g., 'a_posterr_soc' -> 'a')
                base = c.replace('_soc', '').replace('_posterr', '').replace('_postcorr', '')
            elif '_nonsoc' in c:
                soc = 0
                base = c.replace('_nonsoc', '').replace('_posterr', '').replace('_postcorr', '')
            else:
                continue
                
            col_map[c] = {'type': 'crossed', 'base': base, 'soc': soc, 'acc': 0 if is_err else 1}
            continue
    
        # --- 4. Check ERN / CRN (Crossed) ---
        if c.startswith('ERN') or c.startswith('CRN'):
            is_ern = c.startswith('ERN')
            
            # Determine social condition
            if '_soc' in c:
                soc = 1
                suffix = c.replace('ERN', '').replace('CRN', '').replace('_soc', '')
            elif '_nonsoc' in c:
                soc = 0
                suffix = c.replace('ERN', '').replace('CRN', '').replace('_nonsoc', '')
            else:
                continue
    
            # Clean suffix to get base (empty suffix -> 'amplitude')
            base = suffix.strip('_')
            if base == '': base = 'amplitude'
            if base == 'laplacian': base = 'laplacian'
            
            col_map[c] = {'type': 'crossed', 'base': base, 'soc': soc, 'acc': 0 if is_ern else 1}
            continue
    
        # --- 5. Check Generic _err / _corr Suffix (Crossed) ---
        if '_err' in c or '_corr' in c:
            is_err = '_err' in c
            if '_soc' in c:
                soc = 1
                base = c.replace('_soc', '').replace('_err', '').replace('_corr', '')
            elif '_nonsoc' in c:
                soc = 0
                base = c.replace('_nonsoc', '').replace('_err', '').replace('_corr', '')
            else:
                continue
            
            col_map[c] = {'type': 'crossed', 'base': base, 'soc': soc, 'acc': 0 if is_err else 1}
            continue
    
    # --- Transformation ---
    
    melted = df.melt(id_vars=id_cols, var_name='original_col', value_name='value')
    meta = melted['original_col'].map(col_map)
    melted = melted[meta.notna()]
    meta = meta[meta.notna()]
    
    melted['base'] = meta.apply(lambda x: x.get('base'))
    melted['soc'] = meta.apply(lambda x: x.get('soc'))
    melted['acc_spec'] = meta.apply(lambda x: x.get('acc'))
    melted['type'] = meta.apply(lambda x: x.get('type'))
    
    # Split rows: Duplicate Social-only rows
    social_rows = melted[melted['type'] == 'social'].copy()
    crossed_rows = melted[melted['type'] == 'crossed'].copy()
    
    social_0 = social_rows.copy()
    social_0['acc'] = 0
    social_1 = social_rows.copy()
    social_1['acc'] = 1
    
    # Crossed rows use their specific accuracy
    crossed_rows['acc'] = crossed_rows['acc_spec']
    
    final_long = pd.concat([social_0, social_1, crossed_rows], ignore_index=True)
    
    # Pivot
    pivot_index = id_cols + ['soc', 'acc']
    pivot_cols = 'base'
    pivot_values = 'value'
    
    df_long = final_long.pivot_table(index=pivot_index, columns=pivot_cols, values=pivot_values, aggfunc='first')
    
    # Cleanup Index
    df_long.index.names = [n + '_idx' for n in df_long.index.names] # Rename to avoid collision
    df_long = df_long.reset_index()
    df_long.rename(columns={'soc_idx': 'soc', 'acc_idx': 'acc'}, inplace=True)
    df_long.rename(columns={c: c.replace('_idx', '') for c in df_long.columns if c.endswith('_idx')}, inplace=True)
    df_long = df_long.merge(df[['sub'] + redcap_cols].drop_duplicates('sub'), on='sub', how='left')
    # Ensure Integer Types
    df_long['soc'] = df_long['soc'].astype(int)
    df_long['acc'] = df_long['acc'].astype(int)
    
    # Save
    df_long = df_long.rename({"amplitude": "ERN", "laplacian": "ERN_laplacian"}, axis=1)
    df_long.to_csv(f"{csv_output_path}/thrive_long_{session}_{date_time}.csv", index=False)
    # print("Transformation complete. Columns:", df_long.columns.tolist())

list_of_wide_dfs = []
# list_of_long_dfs = []
for session in sessions:
    csv_output_path = f"{analysis_path}/derivatives/csv/{session}/"

    wide_data = pd.read_csv(find_newest_file(f"{csv_output_path}/thrive_wide_{session}*.csv"))
    wide_data.columns = [c + f"_{session}" if (not(session in c) and c!="sub") else c for c in wide_data]
    list_of_wide_dfs.append(wide_data)
    
    # long_data = pd.read_csv(find_newest_file(f"{csv_output_path}/thrive_long_{session}*.csv"))
    # long_data.columns = [c + f"_{session}" if (not(session in c) and c!="sub") else c for c in long_data]
    # list_of_long_dfs.append(long_data)

mega_wide_data = reduce(lambda left, right: pd.merge(left, right, on="sub", how='left'), list_of_wide_dfs)
mega_wide_data.to_csv(f"{analysis_path}/derivatives/csv/thrive_mega_wide_{date_time}.csv", index=False)

# mega_long_data = reduce(lambda left, right: pd.merge(left, right, on="sub", how='left'), list_of_long_dfs)
# mega_long_data.to_csv(f"{analysis_path}/derivatives/csv/thrive_mega_long_{date_time}.csv", index=False)
