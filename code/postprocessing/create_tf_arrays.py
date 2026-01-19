import pandas as pd
import numpy as np
from glob import glob
import scipy.io
import h5py
import os
from tqdm import tqdm
import re
import sys

"""
Script to aggregate individual Time-Frequency (TF) data into group-level arrays.

This script iterates through all subjects and conditions, loads their individual
TF output files (HDF5/MAT), and concatenates them into large matrices (Subject x Channel x Time x Freq).
These aggregate arrays are saved for faster group-level analysis.

Usage:
    python create_tf_arrays.py <session_id>

Arguments:
    session_id (str): The session identifier.
"""

# This code creates big arrays of n_sub * chan * times * freqs data from individual time-frequecy data arrays by concatenating individual participants' data.

#session = "s2_r1"
session = sys.argv[1]

# csv_path = f"derivatives/behavior/{session}/"
output_path = f"/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/preprocessed/TF_arrays/{session}/"

# Define the regex pattern to match 'sub-' followed by digits
pattern = r'sub-\d+'

# tf_files = sorted(glob(f"{data_path}/sub-*{condition}*.mat"))
for measure in [
    "TF",
    "ITPS",
    "ICPS",
    "wPLI"
]:
    print(f"Working on {measure} ... ")
    if measure == "ITPS" or measure == "ICPS":
        key_idx = 1
    else:
        key_idx = -1
    data_path = f"/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/preprocessed/TF_outputs/{session}/resp/seed_1/{measure}/"
    for condition in tqdm(["resp_ns_c_1", "resp_ns_i_0", "resp_ns_i_1",
                           "resp_s_i_1", "resp_s_c_1", "resp_s_i_0"]):
        
        arr_list = []
        subjects_with_data = []
        # Extract subject ids that have TF data
        matched_parts = [
            re.search(pattern, s).group(0) if re.search(pattern, s) else None for s in glob(
                f"{data_path}/sub-*all_eeg_processed_data*{measure}*{condition}*.mat"
            )
        ]
        for sub_id in sorted(matched_parts):
            # check if there is data for that subject for that condition
            try:
                # sort all TF files by sub_id first
                tf_files = sorted(glob(f"{data_path}/{sub_id}*all_eeg_processed_data*{measure}*{condition}*.mat"))
                assert len(tf_files) == 1, "Check your tf_files length!"

                # read TF array
                data_file = h5py.File(tf_files[0])
                key_list = list(data_file.keys())
                data = data_file[key_list[key_idx]]
                # take only actual data with channels * times * freqs
                assert data.shape == (64, 375, 59), "Check your data!"
                arr_list.append(data)
                subjects_with_data.append(sub_id) # that way, the actual data and participant ids will go in the same order
            except: continue
        
        # concatenate all valid subject data for a given condition 
        full_data = np.stack(arr_list, axis=0)
        # make sure the number of arrays is the same as number of subject ids saved previously
        assert full_data.shape[0] == len(subjects_with_data), "Check your data!"
        print(f"# of subjects to have condition {condition}: {len(arr_list)}")
        # save resulting subs * channels * times * freqs
        scipy.io.savemat(f"{output_path}/{measure}_{condition}.mat",
                         {
                             f"{measure}_{condition}": full_data,
                             f"subjects": subjects_with_data,
                         })
