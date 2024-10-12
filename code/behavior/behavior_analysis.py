import sys
import os
import glob
import pandas as pd
import numpy as np
import re
import time
import datetime
pd.options.mode.chained_assignment = None

def convert_to_list_rt(series):
    float_list = []
    for value in series:
        if isinstance(value, str):
            if "," in value.strip("[]"):
                float_list.append([float(v) for v in value.strip("[]").split(",")][0]) # Check if the value is a string
            else:
                float_list.append(float(value.strip("[]"))) # Convert string to float and remove brackets
        elif isinstance(value, list): # Check if the value is a list
            float_list.extend([float(v) for v in value]) # Convert each element of the list to float
        else: # Handle NaN values
            float_list.append(np.nan) # Append NaN if value is NaN
    return float_list

def convert_to_list_resp(series):
    resp_list = []
    for value in series:
        if isinstance(value, str):
            converted_row = list(map(int, re.findall(r'\d+', value)))
            resp_list.append(converted_row)
        else: # Handle NaN values
            resp_list.append(np.nan) # Append NaN if value is NaN
    return resp_list

start = time.time()

input_dataset_path = "/home/data/NDClab/datasets/thrive-dataset/"
output_dataset_path = "/home/data/NDClab/analyses/thrive-theta-ddm/"
data_path = "sourcedata/checked/"
sub_path = "s1_r1/psychopy/"
output_path = "derivatives/behavior/"

date_time = datetime.datetime.now().strftime("%d-%m-%Y_%H_%M_%S")
sys.stdout = open(f"{output_dataset_path}{output_path}{date_time}_log.txt","wt")

n_blocks = 20
n_trials = 40
valid_rt_thresh = 0.150

sub_folders = [i for i in os.listdir(input_dataset_path + data_path) if i.startswith("sub-")]
subjects = sorted([re.findall(r'\d+', item)[0] for item in sub_folders])

processing_log = dict()
summary_columns = [
            "n_trials", "invalid_rt_percent", "skipped_percent",
            "acc", "acc_con", "acc_incon", "rt_con", "rt_incon", "rt_corr", "rt_err",
            "pes", "pea", "peri_acc", "peri_rt", "6_or_more_err",
            ]
processing_log["sub"] = []
processing_log["success"] = []
processing_log["n_trials"] = []

for condition in [0, 1]:
    if condition == 0:
        prefix = "_nonsoc"
    elif condition == 1:
        prefix = "_soc"
    for colname in summary_columns:
        processing_log[colname + prefix] = []

for sub in subjects:
    processing_log["sub"].append(sub)    
    subject_folder = (input_dataset_path + data_path + "sub-" + sub + os.sep + sub_path)
    num_files = len(os.listdir(subject_folder))
    
    if (num_files != 3) and (sub not in ["3000124", "3000008", "3000014"]):
        processing_log["success"].append(0)
        print("sub-{} has a deviation in psychopy data ({} files), skipping ...".format(sub, num_files))
        [processing_log[i].append(np.nan) for i in list(processing_log.keys())[2:]]
        pass
    else:
        print("Processing sub-{}...".format(sub))
        processing_log["success"].append(1)
        
        pattern = "{}sub-{}_arrow-alert-v1-*_psychopy_s1_r1_e1.csv".format(subject_folder, sub)
        filename = glob.glob(pattern)
        data = pd.read_csv(filename[0])
        start_index = data["task_blockText.started"].first_valid_index()
        data = data.iloc[start_index:, :].dropna(subset = "middleStim")
        data = data[data["conditionText"].isin(["Observed", "Alone"])].reset_index(drop=True)
        assert (len(data) == n_blocks * n_trials), "Check your data length!"
        processing_log["n_trials"].append(len(data))
        
        trial_data = data[[
            "target",
            "congruent",
            "accuracy",
            "task_stim_keyResp.rt",
            "task_stim_keyResp.stopped",
            "task_stim_keyResp.keys",
            "conditionText",
        ]]
        trial_data["rt"] = convert_to_list_rt(trial_data["task_stim_keyResp.rt"])
        trial_data.drop("task_stim_keyResp.rt", axis = 1, inplace = True)
        assert (np.sum([type(i) != float for i in trial_data["rt"]]) == 0), "Check your RT!"

        trial_data["resp_direction_R"] = convert_to_list_resp(trial_data["task_stim_keyResp.keys"])
        trial_data.drop("task_stim_keyResp.keys", axis = 1, inplace = True)
        
        trial_data.columns = [
            "target",
            "congruent",
            "accuracy",
            "task_stim_keyResp.stopped",
            "condition",
            "rt",
            "resp_direction_R",
        ]
        trial_data["condition_soc"] = [1 if i == "Observed" else 0 for i in trial_data["condition"]]
        trial_data.drop("condition", axis = 1, inplace = True)
        
        trial_data["target_R"] = [0 if i == "left" else 1 for i in trial_data["target"]]
        trial_data.drop("target", axis = 1, inplace = True)
        
        trial_data["fl_direction_R"] = [
                                        0 if 
                                        (
                                            (trial_data.loc[i, 'target_R'] == 0 and trial_data.loc[i, 'congruent'] == 1) or
                                            (trial_data.loc[i, 'target_R'] == 1 and trial_data.loc[i, 'congruent'] == 0)
                                        )
                                        else 1 if 
                                        (
                                            (trial_data.loc[i, 'target_R'] == 0 and trial_data.loc[i, 'congruent'] == 0) or
                                            (trial_data.loc[i, 'target_R'] == 1 and trial_data.loc[i, 'congruent'] == 1)
                                        )
                                        else None
                                        for i in range(len(trial_data))
                                    ]
        trial_data["valid_rt"] = [0 if i < valid_rt_thresh else 1 for i in trial_data["rt"]]
        trial_data["no_resp"] = [1 if np.isnan(i) else 0 for i in trial_data["rt"]]
        
        trial_data["block_num"] = sum([[i] * n_trials for i in range(1, n_blocks+1)], [])
        trial_data["trial_num"] = [i for i in range(1, len(trial_data)+1)]
        trial_data["first_trial"] = [1 if i == 0 else 0 for i in range(len(trial_data))]
        trial_data["last_trial"] = [1 if i == (len(trial_data)-1) else 0 for i in range(len(trial_data))]

        extra_resp = []
        resp_direction = []
        for i in range(len(trial_data)):
            row = trial_data.loc[i, "resp_direction_R"]
            if type(row) == list:
                if row[0] == 1:
                    resp_direction.append(0)
                elif row[0] == 8:
                    resp_direction.append(1)
                if len(row) > 1:
                    extra_resp.append(1)
                else:
                    extra_resp.append(0)
            elif np.isnan(row):
                extra_resp.append(np.nan)
                resp_direction.append(np.nan)

        trial_data["resp_direction_R"] = resp_direction
        trial_data["extra_resp"] = extra_resp
        
        assert (len(trial_data) == n_blocks * n_trials), "Check your data length!"

        current_cols = trial_data.columns
        for col_name in current_cols:
            trial_data["pre_" + col_name] = "None"
            trial_data["next_" + col_name] = "None"

        # Iterate through each row of the dataframe
        for i in range(len(trial_data)):
            # Check for previous trial (n-1) if it exists and is in the same block
            if i > 0 and (trial_data.loc[i, 'task_stim_keyResp.stopped'] - trial_data.loc[i-1, 'task_stim_keyResp.stopped']) <= 3\
            and trial_data.loc[i, 'valid_rt'] == 1 and trial_data.loc[i, 'no_resp'] == 0:
                for col_name in current_cols:
                    trial_data.loc[i, 'pre_' + col_name] = trial_data.loc[i-1, col_name]
            else:
                for col_name in current_cols:
                    trial_data.loc[i, 'pre_' + col_name] = np.nan
        for i in range(len(trial_data)):
            # Check for next trial (n+1) if it exists and is in the same block
            if i < len(trial_data)-1 and (trial_data.loc[i+1, 'task_stim_keyResp.stopped'] - trial_data.loc[i, 'task_stim_keyResp.stopped']) <= 3\
            and trial_data.loc[i, 'valid_rt'] == 1 and trial_data.loc[i, 'no_resp'] == 0:
                for col_name in current_cols:
                    trial_data.loc[i, 'next_' + col_name] = trial_data.loc[i+1, col_name]
            else:
                for col_name in current_cols:
                    trial_data.loc[i, 'next_' + col_name] = np.nan

        # Check if the string "None" exists anywhere in the DataFrame to make sure all cells were properly populated in the above step
        assert not ((trial_data == "None").any().any()), "Check your data!"

        trial_data.drop(['pre_task_stim_keyResp.stopped', 'next_task_stim_keyResp.stopped'], axis = 1, inplace = True)
        
        trial_data["sub"] = sub
        all_cols = list(trial_data.columns)[:-1]
        all_cols.insert(0, "sub")
        trial_data = trial_data[all_cols]

        trial_data.to_csv(f"{output_dataset_path}{output_path}sub-{sub}_trial_data.csv", index=False)
        condition = []
        for condition in [0, 1]:
            if condition == 0:
                prefix = "_nonsoc"
            elif condition == 1:
                prefix = "_soc"
            condition_data = trial_data[trial_data["condition_soc"] == condition]
            processing_log["n_trials"+prefix].append(len(condition_data))
            processing_log["skipped_percent"+prefix].append(np.round(condition_data["no_resp"].sum() / len(condition_data) * 100, 3))
            processing_log["invalid_rt_percent"+prefix].append(np.round((1 - (sum(condition_data["valid_rt"]) / len(condition_data))) * 100, 3))
            condition_data = condition_data[(condition_data["valid_rt"] == 1)]
            processing_log["6_or_more_err"+prefix].append(1 if len(condition_data[(condition_data["no_resp"] == 0) & (condition_data["accuracy"] == 0)]) >= 6 else 0)
            processing_log["acc"+prefix].append(np.round(condition_data.accuracy.mean(), 3))
            processing_log["acc_con"+prefix].append(np.round(condition_data[condition_data["congruent"] == 1].accuracy.mean(), 3))
            processing_log["acc_incon"+prefix].append(np.round(condition_data[condition_data["congruent"] == 0].accuracy.mean(), 3))
            processing_log["rt_con"+prefix].append(np.round(condition_data[(condition_data["congruent"] == 1) & (condition_data["accuracy"] == 1)]["rt"].mean() * 1000, 3))
            processing_log["rt_incon"+prefix].append(np.round(condition_data[(condition_data["congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"].mean() * 1000, 3))
            processing_log["rt_corr"+prefix].append(np.round(condition_data[(condition_data["congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"].mean() * 1000, 3))
            processing_log["rt_err"+prefix].append(np.round(condition_data[(condition_data["congruent"] == 0) & (condition_data["accuracy"] == 0)]["rt"].mean() * 1000, 3))
            condition_data = condition_data[(condition_data["pre_valid_rt"] == 1) & (condition_data["pre_extra_resp"] == 0)]
            processing_log["pes"+prefix].append(np.round(
                np.log(
                    condition_data[(condition_data["accuracy"] == 1) & (condition_data["pre_accuracy"] == 0) &\
                    (condition_data["pre_congruent"] == 0)]["rt"]
                ).mean()\
                - np.log(
                    condition_data[(condition_data["accuracy"] == 1) & (condition_data["pre_accuracy"] == 1) &\
                    (condition_data["pre_congruent"] == 0)]["rt"]
                ).mean(), 5
            ))
            processing_log["pea"+prefix].append(np.round(
                condition_data[(condition_data["pre_accuracy"] == 0) & (condition_data["pre_congruent"] == 0)]["accuracy"].mean()\
                - condition_data[(condition_data["pre_accuracy"] == 1) & (condition_data["pre_congruent"] == 0)]["accuracy"].mean(), 5
            ))
    
            processing_log["peri_acc"+prefix].append(np.round(
                (
                    condition_data[(condition_data["pre_accuracy"] == 0) & (condition_data["congruent"] == 0) &\
                    (condition_data["pre_congruent"] == 0)]["accuracy"].mean()\
                 - condition_data[(condition_data["pre_accuracy"] == 0) & (condition_data["congruent"] == 1) &\
                    (condition_data["pre_congruent"] == 0)]["accuracy"].mean()
                )\
                - (
                    condition_data[(condition_data["pre_accuracy"] == 1) & (condition_data["congruent"] == 0) &\
                    (condition_data["pre_congruent"] == 0)]["accuracy"].mean()\
                 - condition_data[(condition_data["pre_accuracy"] == 1) & (condition_data["congruent"] == 1) &\
                    (condition_data["pre_congruent"] == 0)]["accuracy"].mean()
                ), 5
            ))
    
            processing_log["peri_rt"+prefix].append(np.round(
                (
                    np.log(
                    condition_data[(condition_data["pre_accuracy"] == 0) & (condition_data["congruent"] == 0) &\
                        (condition_data["pre_congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"]
                    ).mean()\
                 - np.log(
                     condition_data[(condition_data["pre_accuracy"] == 0) & (condition_data["congruent"] == 1) &\
                     (condition_data["pre_congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"]
                 ).mean()
                )\
                - (
                    np.log(
                    condition_data[(condition_data["pre_accuracy"] == 1) & (condition_data["congruent"] == 0) &\
                    (condition_data["pre_congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"]
                    ).mean()\
                 - np.log(
                     condition_data[(condition_data["pre_accuracy"] == 1) & (condition_data["congruent"] == 1) &\
                     (condition_data["pre_congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"]
                 ).mean()
                  ), 5
            ))
        
        print(f"sub-{sub} has been processed")

pd.DataFrame(processing_log).to_csv(f"{output_dataset_path}{output_path}summary.csv", index=False)

[i for i in os.listdir(f"{output_dataset_path}{output_path}") if "sub-" in i]
list_of_ind_csv = []
for df in [i for i in os.listdir(f"{output_dataset_path}{output_path}") if "sub-" in i]:
    list_of_ind_csv.append(pd.read_csv(f"{output_dataset_path}{output_path}{df}"))
full_df = pd.concat(list_of_ind_csv)
# full_df = full_df[(full_df["pre_accuracy"] == 1) | (full_df["pre_accuracy"] == 0)]
full_df.to_csv(f"{output_dataset_path}{output_path}full_df.csv", index = False)

end = time.time()
print(f"Executed time {np.round(end - start, 2)} s")
