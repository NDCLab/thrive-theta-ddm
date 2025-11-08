import os
import re
from glob import glob
from datetime import datetime as dt
from pathlib import Path
import sys

session = sys.argv[1]
dataset_path = "/home/data/NDClab/datasets/thrive-dataset/"
data_path = "derivatives/preprocessed/"
start_date = dt.strptime("18/09/34", "%d/%m/%y")

dir_list = sorted(
    [
        p for p in glob(
    f"{dataset_path}{data_path}sub-*/{session}/eeg/MADE_preprocessing_report*.csv")
    ], # if "ERROR" not in os.path.basename(p)],
                  key=os.path.getmtime, reverse=True
                 )

subs_to_process = []
for path in dir_list:
    # print(path)
    sub_id = re.search(r"sub-(\d+)", path)
    # print(sub_id[1])
    timestamp = os.path.getmtime(path)
    if dt.fromtimestamp(timestamp) < start_date:
        subs_to_process.append(sub_id[1])
subs_to_process = sorted(list(set(subs_to_process)))
# print("/".join(subs_to_process))
print("")
print(f"Subjects already preprocessed in {dataset_path+data_path}: {len(subs_to_process)}")
print("")
print("/".join(subs_to_process))

derivatives = subs_to_process.copy()

dataset_path = "/home/data/NDClab/datasets/thrive-dataset/"
data_path = "sourcedata/checked/"

dir_list = sorted(
    [
        p for p in glob(
    f"{dataset_path}{data_path}sub-*/{session}/eeg/*all_eeg_{session}_e1*.eeg")
    ], # if "ERROR" not in os.path.basename(p)],
                  key=os.path.getmtime, reverse=True
                 )

subs_to_process = []
for path in dir_list:
    # print(path)
    sub_id = re.search(r"sub-(\d+)", path)
    # print(sub_id[1])
    timestamp = os.path.getmtime(path)
    if dt.fromtimestamp(timestamp) < start_date:
        subs_to_process.append(sub_id[1])
subs_to_process = sorted(list(set(subs_to_process)))
# print("/".join(subs_to_process))
print("")
print(f"Subjects total in {dataset_path+data_path}: {len(subs_to_process)}")
sourcedata = subs_to_process.copy()

difference = sorted(list(set(sourcedata) - set(derivatives)))
print("")
print("/".join(difference))
print("")
print(f"Subjects left to preprocess: {len(difference)}")
