
import argparse
import re
from pathlib import Path

def get_args():
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(description="Check behavior processing status.")
    parser.add_argument("session", type=str, help="The session ID (e.g., ses-01)")
    return parser.parse_args()

def extract_sub_id(path_obj):
    """Extracts subject ID (digits) from a Path object."""
    match = re.search(r"sub-(\d+)", str(path_obj))
    #match = re.search(r"sub-(\d+)", path_obj.name)
    return match.group(1) if match else None

def main():
    args = get_args()
    session = args.session

    # --- Configuration ---
    base_path = Path("/home/data/NDClab")
    source_dir = base_path / "datasets/thrive-dataset/sourcedata/checked"
    deriv_dir = base_path / "analyses/thrive-theta-ddm/derivatives/behavior"

    print(f"--- Processing Check for Session: {session} ---\n")

    # --- 1. Find Source Data ---
    # Look for raw PsychoPy files: sub-*/{session}/eeg/*all_eeg_{session}_e1*.eeg
    # Using rglob or glob with specific pattern
    source_pattern = f"sub-*/{session}/psychopy/*{session}_e1*.csv"
    source_files = sorted(source_dir.glob(source_pattern))
    
    source_subs = {extract_sub_id(p) for p in source_files}
    source_subs.discard(None) # Safety cleanup
    
    print(f"Subjects with Source Data: {len(source_subs)}")

    # --- 2. Find Deviations ---
    # Look for deviation txt files
    dev_pattern = f"sub-*/{session}/psychopy/*deviation*.txt"
    dev_files = source_dir.glob(dev_pattern)
    
    dev_subs = {extract_sub_id(p) for p in dev_files}
    dev_subs.discard(None)

    if dev_subs:
        print(f"Subjects with Deviations:   {len(dev_subs)}")
        print(f"IDs with Deviations:       {'/'.join(sorted(dev_subs))}\n")

    # --- 3. Find Successfully Processed Data ---
    # Look for the specific MADE report csv
    report_pattern = f"_trial_data.csv"
    # We look into sub-*/{session}/eeg/ inside derivatives
    processed_files = deriv_dir.glob(f"{session}/sub-*{report_pattern}")
    
    processed_subs = {extract_sub_id(p) for p in processed_files}
    processed_subs.discard(None)
    
    print(f"Subjects Fully Processed: {len(processed_subs)}")

    # --- 4. Calculate Difference (Pending) ---
    # Source - Processed = Pending
    pending_subs = sorted(list(source_subs - processed_subs))
    
    # Identify which of the pending subjects have deviations
    pending_with_devs = sorted(list(set(pending_subs) & dev_subs))

    # --- 5. Output Results ---
    print("-" * 30)
    print(f"Subjects LEFT to process:  {len(pending_subs)}")
    print("-" * 30)
    
    if pending_subs:
        print("/".join(pending_subs))
    else:
        print("No pending subjects. All caught up!")

    if pending_with_devs:
        print("\n" + "!"*30)
        print(f"WARNING: {len(pending_with_devs)} of the pending subjects have deviation files:")
        print("/".join(pending_with_devs))
        print("!"*30)

if __name__ == "__main__":
    main()
