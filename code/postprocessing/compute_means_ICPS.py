import io
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
from glob import glob
import datetime
import time
import re
import h5py
import os
import sys

"""
Script to compute mean Inter-Channel Phase Synchrony (ICPS) values.

This script loads precomputed ICPS arrays (from MAT/HDF5 files), extracts mean values
for specific frequency bands (theta, delta), time windows (early, late), and
electrode clusters (DLPFC, OCC, MOTOR). It merges these values with behavioral
summary data and saves the result to a CSV file.

Usage:
    python compute_means_ICPS.py <session_id>

Arguments:
    session_id (str): The session identifier.
"""

session = sys.argv[1]
dataset_path = "/home/data/NDClab/analyses/thrive-theta-ddm/"

arr_path = f"{dataset_path}/derivatives/preprocessed/TF_arrays/{session}/"
helper_data = h5py.File(
    glob(f"{dataset_path}/derivatives/preprocessed/TF_outputs/{session}/resp/seed_1/TF/sub-*.mat")[0]
)

freqs = helper_data['frequency'][:]
times = helper_data['ds_time'][:]
assert np.max(np.abs(times)) > 50, (
    f"Time unit warning: Max time is {np.max(np.abs(times)):.2f}. "
    "This looks like SECONDS. Convert 'times' to MS (times * 1000)."
)
ch_locs = [str(i) for i in range(1, 65)]

# NOT DIFFERENCE ICPS

matching_files = glob(f"{dataset_path}/derivatives/behavior/{session}/*summary*{session}*.csv")

# Check if any files were found
if not matching_files:
    print("No matching files found.")
else:
    # Find the newest file based on modification time
    new_file_path = max(matching_files, key=os.path.getmtime)
    thrive_data = pd.read_csv(new_file_path)["sub"].to_frame()
    print(f"The newest file is: {new_file_path}")

measures = [
    "ICPS",
]
conditions = [
    "resp_s_i_0",
    "resp_s_i_1",
    # "resp_s_c_1",
    "resp_ns_i_0",
    "resp_ns_i_1",
    # "resp_ns_c_1",
            ]

for m in measures:
    for c in conditions:
        tf_arr = scipy.io.loadmat(f"{arr_path}/{m}_{c}.mat")
        tf_data = tf_arr[f"{m}_{c}"]       
        assert tf_data.shape[1:] == (64, 375, 59), f"Check your {m} data!"
        sub_ids = tf_arr['subjects']
        pattern = r'sub-(\d+)'
        sub_ids = [int(re.search(pattern, i).group(1)) for i in sub_ids]

        for band in [
            "theta",
             # "delta"
                    ]:
            for window in [
                "early",
                "late"
            ]:
                for cluster in [
                    "DLPFC_L",
                    "DLPFC_R",
                    "OCC_L",
                    "OCC_R",
                    "MOTOR_L",
                    "MOTOR_R",
                    #"CENTRAL"
                ]:
                    ch = None

                    if cluster == "DLPFC_L":
                        ch = ['6', '9']
                    elif cluster == "DLPFC_R":
                        ch = ['39', '42']
                    elif cluster == "OCC_L":
                        ch = ['22', '24']
                    elif cluster == "OCC_R":
                        ch = ['53', '55']
                    elif cluster == "MOTOR_L":
                        ch = ['3', '7']
                    elif cluster == "MOTOR_R":
                        ch = ['35', '40']
                    elif cluster == "CENTRAL":
                        ch = ['19', '50']

                    assert ch is not None, f"Cluster {cluster} is not defined in channel map!"

                    if band == "theta":
                        fmin = 4
                        fmax = 7
                    elif band == "delta":
                        fmin = 1
                        fmax = 3
    
                    if window == "early":
                        tmin = 0
                        tmax = 250
                    elif window == "late":
                        tmin = 256
                        tmax = 504
                        
                    fmin_idx = np.argmin(np.abs(freqs-fmin))
                    assert freqs[fmin_idx] == fmin, "Check your freqs!"
                    fmax_idx = np.argmin(np.abs(freqs-fmax))
                    assert freqs[fmax_idx] == fmax, "Check your freqs!"
                    
                    tmin_idx = np.argmin(np.abs(times-tmin))
                    # assert times[tmin_idx] == tmin, "Check your times!"
                    tmax_idx = np.argmin(np.abs(times-tmax))
                    # assert times[tmax_idx] == tmax, "Check your times!"
                    
                    ch_idx = []
                    for channel in ch:
                        if channel in ch_locs:
                            ch_idx.append(ch_locs.index(channel))
                    
                    # sub_idx = scipy.io.loadmat(f"{arr_path}/idx_{c}.mat")["sub_idx"][0]-1 # make it 0-based again
                    tf_df = pd.DataFrame(columns = ["sub", f"{m}_{c}_{band}_{window}_{cluster}"])
                    
                    # for sub_id in sub_idx:
                    for sub_id in range(tf_data.shape[0]):
                        # sub_avg = np.mean(tf_data[sub_id, :, :, :], 0)
                        sub_avg = tf_data[sub_id, :, :, :]
                        assert sub_avg.shape == (64, 375, 59), f"Check your {m} data!"
                        
                        ch_avg = np.mean(sub_avg[ch_idx, :, :], 0)
                        assert ch_avg.shape == (375, 59), f"Check your {m} data!"
                        
                        time_avg = np.mean(ch_avg[tmin_idx:tmax_idx+1, :], 0)
                        assert len(time_avg) == 59 and time_avg.ndim == 1, f"Check your {m} data!"
                        freq_avg = np.mean(time_avg[fmin_idx:fmax_idx+1], 0)
                    
                        tf_df.loc[sub_id, "sub"] = sub_ids[sub_id]
                        tf_df.loc[sub_id, f"{m}_{c}_{band}_{window}_{cluster}"] = freq_avg
                    
                    thrive_data = thrive_data.merge(tf_df, on="sub", how="left")

thrive_data = thrive_data[
[i for i in thrive_data.columns if ("delta" not in i or i == "sub")]
]

# note that this renaming logic would not work correctly if congruent conditions are requested above
colnames = list(thrive_data.columns)
for i, c in enumerate(colnames[1:]):
    i+=1
    splitted_list = c.split("_")
    if splitted_list[2] == "s":
        splitted_list[2] = "soc"
    elif splitted_list[2] == "ns":
        splitted_list[2] = "nonsoc"
    if splitted_list[4] == "0":
        splitted_list[4] = "err"
    elif splitted_list[4] == "1":
        splitted_list[4] = "corr"
    splitted_list[1] = ""
    splitted_list[5] = ""
    splitted_list[3] = ""
    splitted_list = [i for i in splitted_list if i!=""]
    colnames[i] = "_".join(splitted_list)

thrive_data.columns = colnames

thrive_data.to_csv(f"{dataset_path}/derivatives/csv/{session}/thrive_icps_{datetime.datetime.now()}.csv", index=False)
