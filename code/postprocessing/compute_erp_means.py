import mne
import io
import numpy as np
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

dataset_path = "/home/data/NDClab/datasets/thrive-dataset/"
analysis_path = "/home/data/NDClab/analyses/thrive-theta-ddm/"

outputHeader = [
    'id',
    'ERN_soc', 'CRN_soc', 'ERN_nonsoc', 'CRN_nonsoc',
    'ERN_min_CRN_diff_soc', 'ERN_min_CRN_diff_nonsoc',
    # 'PE_error_soc', 'PE_corr_soc', 'PE_error_nonsoc', 'PE_corr_nonsoc',
    # 'PE_err_min_corr_diff_soc', 'PE_err_min_corr_diff_nonsoc'
]

output_data = pd.DataFrame()

clustCell= [
    [i-1 for i in [1, 2, 33, 34]],
    # [i-1 for i in [17, 49, 50, 19, 18]],
]

timeCell = [
    [0, 100], # ERN cluster
    # [300, 500], # PE cluster
]

if laplacian == 1:
    path_to_mat = find_newest_file(f"{analysis_path}/derivatives/preprocessed/erp_check/{session}/thrive_Resp_erps_csd_min_6t_*.mat")
elif laplacian == 0:
    path_to_mat = find_newest_file(f"{analysis_path}/derivatives/preprocessed/erp_check/{session}/thrive_Resp_erps_min_6t_*.mat")

path_to_eeg = glob(f"{dataset_path}/derivatives/preprocessed/sub-3000001/{session}/eeg/sub-3000001_all_eeg_processed_data_{session}_e1.set")[0]

mat = scipy.io.loadmat(path_to_mat)
allData = mat['erpDat_data']

# take IDs from EEG (all people > 6 trials)
sub_from_eeg = [int(mat["erpDat_subIds"][i].item()[0]) for i in range(len(mat["erpDat_subIds"]))] 

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

# %round EEG.times to nearest whole ms to make easier to work with
# EEG.times = round(EEG.times);

output_data[outputHeader[0]] = sub_from_eeg

# initialize index var at 1 because i=0 is the column for subject ids
i = 1
for comp in range(len(clustCell)):

    cluster= clustCell[comp]
    times = timeCell[comp]

    compStartTime = times[0] # in ms
    compEndTime = times[1] # in ms

    compStartIdx = np.argmin(np.abs(EEG_times-compStartTime))
    compEndIdx = np.argmin(np.abs(EEG_times-compEndTime))

    s_resp_incon_error_avgTime = np.nanmean(newData[:, 0:1, :, compStartIdx:compEndIdx+1], 3)
    s_resp_incon_corr_avgTime = np.nanmean(newData[:, 1:2, :, compStartIdx:compEndIdx+1], 3)
    ns_resp_incon_error_avgTime = np.nanmean(newData[:, 2:3, :, compStartIdx:compEndIdx+1], 3)
    ns_resp_incon_corr_avgTime = np.nanmean(newData[:, 3:4, :, compStartIdx:compEndIdx+1], 3)

    # average cluster of interest
    s_resp_incon_error_avgTimeClust = np.nanmean(s_resp_incon_error_avgTime[:, :, cluster], 2)
    s_resp_incon_corr_avgTimeClust = np.nanmean(s_resp_incon_corr_avgTime[:, :, cluster], 2)
    ns_resp_incon_error_avgTimeClust = np.nanmean(ns_resp_incon_error_avgTime[:, :, cluster], 2)
    ns_resp_incon_corr_avgTimeClust = np.nanmean(ns_resp_incon_corr_avgTime[:, :, cluster], 2)

    # compute difference scores (not included to the final csv)
    s_resp_incon_error_avgTimeClust_diff = s_resp_incon_error_avgTimeClust - s_resp_incon_corr_avgTimeClust
    ns_resp_incon_error_avgTimeClust_diff = ns_resp_incon_error_avgTimeClust - ns_resp_incon_corr_avgTimeClust

    output_data[outputHeader[i]] = s_resp_incon_error_avgTimeClust
    output_data[outputHeader[i+1]] = s_resp_incon_corr_avgTimeClust
    output_data[outputHeader[i+2]] = ns_resp_incon_error_avgTimeClust
    output_data[outputHeader[i+3]] = ns_resp_incon_corr_avgTimeClust
    output_data[outputHeader[i+4]] = s_resp_incon_error_avgTimeClust_diff
    output_data[outputHeader[i+5]] = ns_resp_incon_error_avgTimeClust_diff
    i+=6

output_data
output_data = output_data.iloc[:, :5]

if laplacian == 1:
    output_data.columns = [i + "_laplacian" if i != "id" else i for i in output_data.columns]
    output_data = output_data.rename({"id": "sub"}, axis=1)
    output_data.to_csv(f"{analysis_path}/derivatives/csv/{session}/thrive_erp_laplacian_{datetime.now()}.csv", index=False)
elif laplacian == 0:
    output_data = output_data.rename({"id": "sub"}, axis=1)
    output_data.to_csv(f"{analysis_path}/derivatives/csv/{session}/thrive_erp_{datetime.now()}.csv", index=False)
