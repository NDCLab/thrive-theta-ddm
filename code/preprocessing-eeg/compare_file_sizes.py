
import os
from pathlib import Path
import pandas as pd

"""
Script to compare file sizes between two directory trees.

Useful for verifying data integrity after copying or moving files, specifically
checking if processed files in the analysis directory match those in the dataset directory.
"""

def compare_subject_files(root_path_a, root_path_b, specific_subpath="**/*"):
    """
    Iterates through subject folders (sub-*) in root_path_a, finds the corresponding
    subject folder in root_path_b, and compares file sizes.

    Args:
        root_path_a (str): The base path containing the sub-XXXX folders (Source).
        root_path_b (str): The base path containing the sub-XXXX folders (Target).
        specific_subpath (str): Glob pattern to limit search inside subject folders.
                                Use "**/*" for everything.
                                Use "**/eeg/*" to only look inside eeg folders.

    Returns:
        pd.DataFrame: A DataFrame containing details of any size mismatches or missing files.
    """
    base_a = Path(root_path_a)
    base_b = Path(root_path_b)
    
    data = []
    
    print(f"--- Scanning ---")
    print(f"Source: {base_a}")
    print(f"Target: {base_b}\n")

    # 1. Find all subject directories in the source root
    # We assume folders start with "sub-" based on your description
    subjects = [d for d in base_a.iterdir() if d.is_dir() and d.name.startswith("sub-")]

    if not subjects:
        print("No 'sub-*' directories found in source path.")
        return pd.DataFrame()

    for sub_dir_a in subjects:
        sub_name = sub_dir_a.name  # e.g., "sub-3000001"
        sub_dir_b = base_b / sub_name
        
        # Check if this subject exists in the target location
        if not sub_dir_b.exists():
            print(f"[MISSING SUBJECT] {sub_name} not found in target root.")
            continue

        # 2. Walk through files inside this subject's folder
        # specific_subpath allows you to filter for 's1_r1/eeg' if needed
        for file_a in sub_dir_a.glob(specific_subpath):
            if file_a.is_file():
                # Get the path relative to the subject folder
                # e.g., if file is .../sub-001/s1_r1/eeg/data.mat
                # relative path is "s1_r1/eeg/data.mat"
                rel_path = file_a.relative_to(sub_dir_a)
                
                # Construct the expected path in the second directory
                file_b = sub_dir_b / rel_path
                
                if file_b.exists():
                    size_a = file_a.stat().st_size
                    size_b = file_b.stat().st_size
                    
                    if size_a != size_b:
                        data.append({
                            "Subject": sub_name,
                            "Relative_Path": str(rel_path),
                            "Status": "SIZE_MISMATCH",
                            "Size_A": size_a,
                            "Size_B": size_b,
                            "Diff": size_a - size_b
                        })
                else:
                    data.append({
                        "Subject": sub_name,
                        "Relative_Path": str(rel_path),
                        "Status": "MISSING_IN_TARGET",
                        "Size_A": file_a.stat().st_size,
                        "Size_B": None,
                        "Diff": None
                    })

    return pd.DataFrame(data)

# --- Configuration ---
if __name__ == "__main__":
    # PATH 1: The 'analyses' folder
    DIR_A = "/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/preprocessed"
    
    # PATH 2: The 'datasets' folder
    DIR_B = "/home/data/NDClab/datasets/thrive-dataset/derivatives/preprocessed"

    # search_pattern:
    # Use "**/*" to check every single file under the subject.
    # Use "*/eeg/*" if you specifically want to check the EEG folders inside any session.
    # Use "s1_r1/eeg/*" if you only want that specific session/run.
    PATTERN = "*/eeg/*" 

    df = compare_subject_files(DIR_A, DIR_B, PATTERN)

    if not df.empty:
        print("\n--- ISSUES FOUND ---")
        # Set pandas to display full path width
        pd.set_option('display.max_colwidth', None) 
        print(df.sort_values(by="Subject").to_string(index=False))
        
        # Optional: Save to CSV
        # df.to_csv("file_mismatches.csv", index=False)
    else:
        print("\nSuccess: All files matching the pattern have identical sizes.")
