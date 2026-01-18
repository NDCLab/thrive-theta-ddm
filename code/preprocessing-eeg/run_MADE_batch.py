
import subprocess
import sys
import time
import os
import re
from glob import glob
from pathlib import Path

# --- Configuration ---
if len(sys.argv) < 2:
    print("Usage: python script.py <session>")
    sys.exit(1)

session = sys.argv[1]
dataset = "thrive-dataset"
# Using absolute paths as base strings to ensure exact match with original glob logic
base_path = "/home/data/NDClab/datasets/" + dataset
processed_dir = f"{base_path}/derivatives/preprocessed/"
source_dir = f"{base_path}/sourcedata/checked/"

sub_id_pattern = re.compile(r"sub-(\d+)")

def get_subjects(glob_pattern):
    """
    Globs files and returns a set of unique subject IDs.
    """
    subjects = set()
    files = glob(glob_pattern)
    
    for path in files:
        # Extract ID
        match = sub_id_pattern.search(path)
        if match:
            subjects.add(match.group(1))
    return subjects

# --- 1. Get Lists ---

# A. Processed (Derivatives)
# Original pattern: sub-*/{session}/eeg/MADE_preprocessing_report_all_eeg_{session}_e1.csv
pat_processed = f"{processed_dir}sub-*/{session}/eeg/MADE_preprocessing_report_all_eeg_{session}_e1.csv"

# you may argue that it's better to use actual fully processed .set files but those won't include cases when no usable data was found during preproessing, therefore, we will use processing reports
#pat_processed = f"{processed_dir}sub-*/{session}/eeg/sub-*all_eeg_processed_data_{session}_e1*.set"

subs_processed = get_subjects(pat_processed)

# B. Total Source (Source Data)
# Original pattern: sub-*/{session}/eeg/*all_eeg_{session}_e1*.eeg
pat_source = f"{source_dir}sub-*/{session}/eeg/*all_eeg_{session}_e1*.eeg"
subs_source = get_subjects(pat_source)

# C. Deviations
# Original pattern: sub-*/{session}/eeg/*deviation*.txt
pat_deviations = f"{source_dir}sub-*/{session}/eeg/*deviation*.txt"
subs_deviations = get_subjects(pat_deviations)

# --- 2. Calculate Subsets ---

# Subjects left to process = (All Source) - (Already Processed)
subs_left = subs_source - subs_processed

# Left AND have deviations = (Left) INTERSECTION (Deviations)
subs_left_with_dev = subs_left.intersection(subs_deviations)

# Left AND NO deviations = (Left) DIFFERENCE (Deviations)
subs_left_no_dev = subs_left - subs_deviations

# --- 3. Output ---

def print_group(title, data_set, show_ids=True):
    sorted_list = sorted(list(data_set))
    print("")
    print(f"{title}: {len(sorted_list)}")
    if show_ids and len(sorted_list) > 0:
        print("/".join(sorted_list))

# 1. How many subjects were processed (without listing IDs)
print_group(f"Subjects already preprocessed in {processed_dir}", subs_processed, show_ids=False)

# 2. How many subjects left to process with their ids
print_group("Subjects left to preprocess", subs_left, show_ids=True)

# 3. How many subjects left to process with deviations with their ids
print_group("Subjects left to preprocess & have deviations", subs_left_with_dev, show_ids=True)

# 4. How many subjects left to process without deviations with their ids
print_group("Subjects left to preprocess WITHOUT deviations", subs_left_no_dev, show_ids=True)

print("")
# Define your parameter space
print("Input subjects to process in the form XXXXXXX/XXXXXXX/XXXXXXX")

subjects_to_process = str(input())

# Construct the sbatch command with --export
# ALL inherits the default environment, then we add our custom vars
cmd = [
    "sbatch",
    f"--export=ALL,DATASET={dataset},SUBS={subjects_to_process},SESS={session}",
    "eeg_processing_batch.sub"
    ]

# Execute the command
result = subprocess.run(cmd, check=True, capture_output=True, text=True)

if result.returncode == 0:
    print(f"Submitted: {session} ({len(subjects_to_process.split('/'))} subjects): {subjects_to_process} -> {result.stdout.strip()}")
else:
    print(f"Error submitting {session}: {result.stderr}")
