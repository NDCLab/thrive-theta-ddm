
import os
import glob
import pandas as pd
import sys

"""
Script to aggregate deviation reports from subject directories.

This script scans for 'deviation.txt' files in the EEG directories of subjects,
reads their content, and compiles a CSV report listing deviations and the number
of EEG files present for each subject.

Usage:
    python deviation_output.py <session_id>

Arguments:
    session_id (str): The session identifier.
"""

session = sys.argv[1]
base_dir = "/home/data/NDClab/datasets/thrive-dataset/sourcedata/checked/"
subjects = sorted(glob.glob(os.path.join(base_dir, "sub-3*")))

data = []

for sub_path in subjects:
    sub_id = os.path.basename(sub_path)
    eeg_path = os.path.join(sub_path, "s1_r1", "eeg")
    
    # Skip if EEG folder doesn't exist
    if not os.path.isdir(eeg_path):
        continue

    # Find deviation file(s)
    deviation_files = glob.glob(os.path.join(eeg_path, "*deviation.txt"))
    if not deviation_files:
        continue  # Skip subjects without deviation.txt

    # Read deviation text (if multiple files, join)
    deviations = []
    for dev_file in deviation_files:
        with open(dev_file, "r") as f:
            deviations.append(f.read().strip())
    deviation_text = "\n".join(deviations)

    # Collect EEG file info (excluding .txt)
    eeg_files = [f for f in os.listdir(eeg_path)
                 if os.path.isfile(os.path.join(eeg_path, f)) and not f.endswith(".txt")]
    
    data.append({
        "sub": sub_id.split("-")[1],
        "deviation": deviation_text,
        "num_of_eeg_files": len(eeg_files),
        "names_of_eeg_files": ", ".join(sorted(eeg_files))
    })

df = pd.DataFrame(data)

# Optional: display or save to CSV
print(df)
df.to_csv(f"subjects_with_deviations_{session}.csv", index=False)

