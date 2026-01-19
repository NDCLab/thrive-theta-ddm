import sys
import subprocess
import time
import pandas as pd
from glob import glob
from pathlib import Path
import pandas as pd

"""
Script to submit DDM fitting jobs to the SLURM scheduler.

This script identifies subjects that have not yet been fitted for the DDM model
and submits batch jobs for them. It allows for submitting jobs in batches.

Usage:
    python run_ddm_batch.py <session_id>

Arguments:
    session_id (str): The session identifier.
"""

session = sys.argv[1]
#session = "s2_r1"
data_dir = "/home/data/NDClab/analyses/thrive-theta-ddm/"
thrive_id_soc = pd.read_csv(f"{data_dir}/derivatives/behavior/{session}/thrive_data_soc.csv")
thrive_id_nonsoc = pd.read_csv(f"{data_dir}/derivatives/behavior/{session}/thrive_data_nonsoc.csv")

#already_fitted = pd.read_csv(f"{data_dir}/code/ddm/fitted_subjects_{session}_2025_12_27_23_54_17.csv")
#already_fitted = pd.read_csv(f"{data_dir}/code/ddm/fitted_subjects_{session}_2025_12_27_23_54_17.csv")

# Define path and pattern
search_dir = Path(f"{data_dir}/code/ddm/")
pattern = f"fitted_subjects_{session}_*.csv"

batch_size = 1
idDat = pd.concat([thrive_id_soc, thrive_id_nonsoc])

# Get all matching files, sort them (works because date format is YYYY_MM_DD...), and take the last one
files = sorted(search_dir.glob(pattern))

if files:
    latest_file = files[-1]
    print(f"Loading: {latest_file.name}")
    already_fitted = pd.read_csv(latest_file)
    print("Continue with this file: ")
    decision = int(input())
    if decision == 1:
        idDat = idDat[~idDat["sub"].isin(already_fitted["sub"])]
        idDat = idDat.sort_values(by="sub")
        total_n = len(idDat["sub"].unique())

        # Define the argument pairs for start_idx and end_idx
        first_idx = list(range(1, total_n+1, batch_size))
        second_idx = list(range(batch_size, total_n+batch_size, batch_size))

        arg_list = [(i, j) for i, j in zip(first_idx, second_idx)]

        # Path to your SLURM script
        slurm_script = "fit_ddm_batch.sub"

        # Iterate over the argument pairs and submit jobs
        for start_idx, end_idx in arg_list:
            # Construct the sbatch command with arguments
            sbatch_command = [
                "sbatch",
                slurm_script,
                str(start_idx),
                str(end_idx),
                session
            ]

            # Submit the job using subprocess.run
            try:
                result = subprocess.run(sbatch_command, check=True, text=True, capture_output=True)
                print(f"Job submitted successfully: {result.stdout.strip()}")
            except subprocess.CalledProcessError as e:
                print(f"Error submitting job: {e.stderr.strip()}")
            time.sleep(4)

    else:
        print("Exit")

else:
    print("No matching files found.")
