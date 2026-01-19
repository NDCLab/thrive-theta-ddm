import argparse
from pathlib import Path
import sys

"""
Script to verify the existence and count of preprocessed EEG files.

This script scans the derivatives directory for preprocessed EEG data (.set and .fdt files)
and alerts if a subject has more than the expected number of files, which might indicate
duplicates or processing issues.

Usage:
    python check_preprocessed_files.py [root_dir]

Arguments:
    root_dir (str, optional): Root directory of the dataset. Defaults to standard path.
"""

def check_eeg_files(root_dir):
    """
    Scans the directory structure for EEG files and reports counts per subject.

    Args:
        root_dir (str): The root directory to start the scan from.
    """
    root_path = Path(root_dir)
    
    if not root_path.exists():
        print(f"Error: Directory '{root_dir}' does not exist.")
        sys.exit(1)

    # Walk through the directory structure
    # We look for patterns matching: derivatives/preprocessed/sub-*/s1_r1/eeg
    # Adjust glob pattern if directory depth varies
    target_pattern = "derivatives/preprocessed/sub-*/s1_r1/eeg"
    
    # Use rglob if the 'derivatives' folder isn't the immediate child of root_dir
    # or simple glob if root_dir points directly to the dataset root
    found_any = False
    
    print(f"Scanning {root_path} for pattern: {target_pattern}...")
    
    for eeg_dir in sorted(root_path.glob(target_pattern)):
        found_any = True
        
        # Extract subject ID from the parent path (sub-XXXXXXX)
        # Structure: .../sub-XXXXXXX/s1_r1/eeg
        # eeg_dir.parts[-3] should be sub-XXXXXXX
        try:
            sub_id = eeg_dir.parts[-3] 
        except IndexError:
            print(f"Skipping malformed path: {eeg_dir}")
            continue
            
        if not sub_id.startswith("sub-"):
            continue

        # Find .set and .fdt files containing the subject ID
        set_files = list(eeg_dir.glob(f"*{sub_id}*.set"))
        fdt_files = list(eeg_dir.glob(f"*{sub_id}*.fdt"))
        
        count_set = len(set_files)
        count_fdt = len(fdt_files)
        
        # Check condition: More than 3 sets (where a "set" implies a .set file)
        # Note: You can adjust logic if you strictly need pairs. 
        # Here we flag if .set files > 3 OR .fdt files > 3
        if count_set > 3 or count_fdt > 3:
            print(f"\n[ALERT] {sub_id}")
            print(f"Path: {eeg_dir}")
            print(f"Count: {count_set} .set files, {count_fdt} .fdt files")
            print("Files found:")
            for f in set_files + fdt_files:
                print(f" - {f.name}")

    if not found_any:
        print("No directories matching the hierarchy found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check EEG file counts per subject.")
    parser.add_argument("root_dir", nargs="?", default="/home/data/NDClab/datasets/thrive-dataset", 
                        help="Root directory of the dataset (default: /home/data/NDClab/datasets/thrive-dataset)")
    
    args = parser.parse_args()
    check_eeg_files(args.root_dir)
