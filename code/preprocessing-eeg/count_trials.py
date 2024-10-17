import mne
import io
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
from glob import glob
import datetime
import time

def eeg_point2lat(lat_array, epoch_array, srate, timewin=None, timeunit=1):
    """
    Convert latency in data points to latency in ms relative to the time locking.

    Parameters:
    lat_array (array-like): Latency array in data points assuming concatenated data epochs.
    epoch_array (array-like or None): Epoch number corresponding to each latency value.
    srate (float): Data sampling rate in Hz.
    timewin (list or None): [min, max] time limits in 'timeunit' units. Default is None.
    timeunit (float): Time unit in seconds. Default is 1 (seconds).

    Returns:
    numpy.ndarray: Converted latency values (in 'timeunit' units) for each epoch.
    """
    
    lat_array = np.array(lat_array)
    
    if epoch_array is None:
        epoch_array = np.ones_like(lat_array)
    else:
        epoch_array = np.array(epoch_array)

    if timewin is None:
        timewin = [0, 0]
    
    if len(lat_array) != len(epoch_array):
        if len(epoch_array) != 1:
            raise ValueError("Latency and epoch arrays must have the same length")
        else:
            epoch_array = np.ones_like(lat_array) * epoch_array[0]
    
    if len(timewin) != 2:
        raise ValueError("Timelimits array must have length 2")
    
    timewin = np.array(timewin) * timeunit
    
    if len(timewin) == 2:
        pnts = int((timewin[1] - timewin[0]) * srate + 1)
        pnts = (timewin[1] - timewin[0]) * srate + 1
    else:
        pnts = 0
    
    newlat = ((lat_array - (epoch_array - 1) * pnts - 1) / srate + timewin[0]) / timeunit
    
    return np.round(newlat * 1e9) / 1e9

def disp_diff_arr(arr1, arr2, round = 5):
    mask = np.round(arr1, round) != np.round(arr2, round)
    
    diff_elements_1 = arr1[mask]
    
    diff_elements_2 = arr2[mask]
    
    # Get the indices where differences occur
    diff_indices = np.argwhere(mask)
    print(len(diff_indices))
    for idx, (val1, val2) in enumerate(zip(diff_elements_1, diff_elements_2)):
        print(f"Difference at {diff_indices[idx]}: {np.round(val1, round)} vs {np.round(val2, round)}")

sys.stdout = open(f"{analysis_path}derivatives/preprocessed/erp_check/{datetime.datetime.now()_log.txt","wt")

trial_data = dict({
        "sub": [],
        "s_resp_incon_error": [],
        "s_resp_incon_corr": [],
        "ns_resp_incon_error": [],
        "ns_resp_incon_corr": [],
        "s_stim_incon_corr": [],
        "s_stim_con_corr": [],
        "ns_stim_incon_corr": [],
        "ns_stim_con_corr": [],
})

dataset_path = "/home/data/NDClab/datasets/thrive-dataset/"
analysis_path = "/home/data/NDClab/analyses/thrive-theta-ddm/"

sub_ids = sorted([i.split("/")[-1] for i in glob(
        f"{dataset_path}derivatives/preprocessed/sub-*")])

list_of_eeg_file = sorted(
    glob(
        f"{dataset_path}derivatives/preprocessed/*/s1_r1/eeg/*all_eeg_processed_data*.set")
)

start = time.time()

for file_idx, filename in enumerate(list_of_eeg_file):
    sub_id = sub_ids[file_idx].split("-")[-1]
    trial_data["sub"].append(sub_id)
    EEG = scipy.io.loadmat(filename, squeeze_me=True, struct_as_record=False)["EEG"]
    EEG_mne = mne.io.read_epochs_eeglab(filename, verbose = 'ERROR',)
    
    events = EEG.event
    n_times = EEG.pnts
    sr = EEG.srate
    num_ch = EEG.nbchan

    drop_idx = []
    for i in range(len(events)):
        latency = eeg_point2lat(
            [events[i].latency],
            [events[i].epoch],
            sr,
            timewin = [EEG.xmin*1000, EEG.xmax*1000],
            timeunit = 1e-3,
             )
        if latency >= -.1 and latency <= .1:
            drop_idx.append(i)
    
    events = [ev for ev in events if list(events).index(ev) in drop_idx]
    print(f"sub-{}: {len(events)} good events were found!")
    
    trial_data["s_resp_incon_error"].append(len(
        [ev for ev in events if\
        (ev.observation == "s") & (ev.eventType == "resp") & (ev.congruency == "i")\
        & (ev.accuracy == 0) & (ev.responded == 1) & (ev.validRt == 1)
    ]
    ))
    
    trial_data["s_resp_incon_corr"].append(len(
        [ev for ev in events if\
        (ev.observation == "s") & (ev.eventType == "resp") & (ev.congruency == "i")\
        & (ev.accuracy == 1) & (ev.responded == 1) & (ev.validRt == 1)
    ]
    ))
    
    trial_data["ns_resp_incon_error"].append(len(
        [ev for ev in events if\
        (ev.observation == "ns") & (ev.eventType == "resp") & (ev.congruency == "i")\
        & (ev.accuracy == 0) & (ev.responded == 1) & (ev.validRt == 1)
    ]
    ))
    
    trial_data["ns_resp_incon_corr"].append(len(
        [ev for ev in events if\
        (ev.observation == "ns") & (ev.eventType == "resp") & (ev.congruency == "i")\
        & (ev.accuracy == 1) & (ev.responded == 1) & (ev.validRt == 1)
    ]
    ))
    
    trial_data["s_stim_incon_corr"].append(len(
        [ev for ev in events if\
        (ev.observation == "s") & (ev.eventType == "stim") & (ev.congruency == "i")\
        & (ev.accuracy == 1) & (ev.responded == 1) & (ev.validRt == 1)
    ]
    ))
    
    trial_data["s_stim_con_corr"].append(len(
        [ev for ev in events if\
        (ev.observation == "s") & (ev.eventType == "stim") & (ev.congruency == "c")\
        & (ev.accuracy == 1) & (ev.responded == 1) & (ev.validRt == 1)
    ]
    ))
    
    trial_data["ns_stim_incon_corr"].append(len(
        [ev for ev in events if\
        (ev.observation == "ns") & (ev.eventType == "stim") & (ev.congruency == "i")\
        & (ev.accuracy == 1) & (ev.responded == 1) & (ev.validRt == 1)
    ]
    ))
    
    trial_data["ns_stim_con_corr"].append(len(
        [ev for ev in events if\
        (ev.observation == "ns") & (ev.eventType == "stim") & (ev.congruency == "c")\
        & (ev.accuracy == 1) & (ev.responded == 1) & (ev.validRt == 1)
    ]
    ))

end = time.time()
print(f"Executed time {np.round(end - start, 2)} s")

trial_data.to_csv(f"{analysis_path}derivatives/preprocessed/erp_check/thrive_trialCounts_RespAndStim_{datetime.datetime.now()}_py.csv", index = False)