import mne
import io
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
from glob import glob
import datetime
import time
import h5py
import sys
import os

ndcl_module_path = '/home/data/NDClab/analyses/thrive-theta-ddm/code/python/'
sys.path.append(ndcl_module_path)

from ndclab_py import *

session = sys.argv[1]
laplacian = int(sys.argv[2])

analysis_path = "/home/data/NDClab/analyses/thrive-theta-ddm/"
fig_path = f"{analysis_path}/derivatives/figs/{session}/"

if laplacian == 1:
    path_to_mat = find_newest_file(f"{analysis_path}/derivatives/preprocessed/erp_check/{session}/thrive_Resp_erps_csd_min_6t_*.mat")
    plot_name = "_csd"
elif laplacian == 0:
    path_to_mat = find_newest_file(f"{analysis_path}/derivatives/preprocessed/erp_check/{session}/thrive_Resp_erps_min_6t_*.mat")
    plot_name = ""

path_to_eeg = glob(f"{analysis_path}/derivatives/preprocessed/csd_data/{session}/sub-3000001_all_eeg_processed_data_{session}_e1.set")[0]
csv_path = f"{analysis_path}/derivatives/csv/{session}/"

component = "ERN"

if component == "ERN":
    chan = [1, 2, 5, 37, 34] # define channel cluster as per layout
    xtitle = "Time Relative to Response (ms)"
    xlim = [-400, 500]
    if "csd" in path_to_mat:
        ylim = [-40, 5]
    else:
        ylim = [-8, 8]
elif component == "Pe":
    chan = [17, 49, 50, 19, 18] # define channel cluster as per layout
    xtitle = "Time Relative to Response (ms)"
    xlim = [-400, 800]
    ylim = [-20, 20]

mat = scipy.io.loadmat(path_to_mat)
allData = mat['erpDat_data']

# take IDs from EEG (all people > 6 trials)
sub_from_eeg = [int(mat["erpDat_subIds"][i].item()[0]) for i in range(len(mat["erpDat_subIds"]))] 

# take IDs from fully processed behavioral data (checked for accuracy, validRT, missed responses) separately for each condition
valid_subject_path = f"/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/behavior/{session}/"
thrive_data_soc = pd.read_csv(find_newest_file(f"{valid_subject_path}/thrive_valid_eeg_soc.csv"))
thrive_data_nonsoc = pd.read_csv(find_newest_file(f"{valid_subject_path}/thrive_valid_eeg_nonsoc.csv"))

sub_soc = list(thrive_data_soc["sub"])
sub_nonsoc = list(thrive_data_nonsoc["sub"])

# find overlapping IDs, i.e., people with valid behavioral data who are also present in the EEG
sub_ind_nonsoc = []
for sub_id in sub_nonsoc:
    if int(sub_id) in sub_from_eeg:
        sub_ind_nonsoc.append(sub_from_eeg.index(int(sub_id)))

sub_ind_soc = []
for sub_id in sub_soc:
    if int(sub_id) in sub_from_eeg:
        sub_ind_soc.append(sub_from_eeg.index(int(sub_id)))

EEG = mne.io.read_epochs_eeglab(path_to_eeg, verbose=False)

EEG_times = EEG.times * 1000
startTime = -400
endTime = -200

startIdx = np.argmin(np.abs(EEG_times-startTime)) # get start index for baseline
endIdx = np.argmin(np.abs(EEG_times-endTime)) # get end index for baseline

allBase = np.squeeze(np.nanmean(allData[:, :, :, startIdx:endIdx+1], 3))
allBase = np.nanmean(allData[:, :, :, startIdx:endIdx+1], 3)
newData = np.zeros_like(allData)

for i in range(allData.shape[3]):
    newData[:, :, :, i] = allData[:, :, :, i] - allBase # baseline correction

chan = [i - 1 for i in chan] # make it work for python
chan = newData[:, :, chan, :] # take data from that cluster
chan = np.nanmean(chan, 2, keepdims=True) # average data across that cluster

# 1. Filter Social Condition (Requires valid data for BOTH Cond 0 and Cond 1)
final_ind_soc = []
for idx in sub_ind_soc:
    # Check if all values in the epoch are NaN for Error (0) or Correct (1)
    # If .all() is True, it means the entire epoch is NaN (missing data)
    is_missing_err = np.isnan(chan[idx, 0, :, :]).all()
    is_missing_corr = np.isnan(chan[idx, 1, :, :]).all()
    
    # Only keep subject if BOTH conditions have data
    if not is_missing_err and not is_missing_corr:
        final_ind_soc.append(idx)

# 2. Filter Nonsocial Condition (Requires valid data for BOTH Cond 2 and Cond 3)
final_ind_nonsoc = []
for idx in sub_ind_nonsoc:
    is_missing_err = np.isnan(chan[idx, 2, :, :]).all()
    is_missing_corr = np.isnan(chan[idx, 3, :, :]).all()
    
    if not is_missing_err and not is_missing_corr:
        final_ind_nonsoc.append(idx)

# --- DATA SLICING & AVERAGING ---

# Use the new filtered lists (final_ind_*) to slice the data
s_resp_incon_error = chan[final_ind_soc, 0:1, :, :] 
s_resp_incon_corr = chan[final_ind_soc, 1:2, :, :]
ns_resp_incon_error = chan[final_ind_nonsoc, 2:3, :, :]
ns_resp_incon_corr = chan[final_ind_nonsoc, 3:4, :, :]

# Average data for each condition across subjects
# Note: Since we filtered out the NaNs above, np.mean and np.nanmean will now give the same result,
# but keeping nanmean is safer.
s_resp_incon_error_Mean = np.squeeze(np.nanmean(s_resp_incon_error, 0))
s_resp_incon_corr_Mean = np.squeeze(np.nanmean(s_resp_incon_corr, 0))
ns_resp_incon_error_Mean = np.squeeze(np.nanmean(ns_resp_incon_error, 0))
ns_resp_incon_corr_Mean = np.squeeze(np.nanmean(ns_resp_incon_corr, 0))

# --- PLOTTING ---

plt.rcParams['figure.dpi'] = 900
plt.figure(figsize=(14,8))
# Updated labels to use len(final_ind_*) so the legend shows the correct N
plt.plot(EEG_times, s_resp_incon_error_Mean, 'red', label = f"Social-Error (n={len(final_ind_soc)})")
plt.plot(EEG_times, s_resp_incon_corr_Mean, 'blue', label = f"Social-Correct (n={len(final_ind_soc)})")
plt.plot(EEG_times, ns_resp_incon_error_Mean, 'red', linestyle="dotted", label = f"Alone-Error (n={len(final_ind_nonsoc)})")
plt.plot(EEG_times, ns_resp_incon_corr_Mean, 'blue', linestyle="dotted", label = f"Alone-Correct (n={len(final_ind_nonsoc)})")
plt.legend()
plt.grid(True)
plt.xlim(xlim[0], xlim[1])
plt.ylim(ylim[0], ylim[1])
plt.ylabel('Amplitude in µV', fontsize=14)
plt.xlabel(xtitle, fontsize=14)
plt.savefig(f"{fig_path}{component}_avg{plot_name}_{session}.png")
plt.show()
