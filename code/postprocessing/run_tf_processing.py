import subprocess
import time

"""
Script to submit batch jobs for Time-Frequency mean computation.

This script submits SLURM jobs to run `compute_means_TF.py` and `compute_means_ICPS.py`
for defined sessions.

Usage:
    python run_tf_processing.py
"""

# Define your iteration parameters
sessions = [
"s1_r1",
#"s2_r1",
#"s3_r1"
]

slurm_scripts = [
"compute_means_TF.sub",
"compute_means_ICPS.sub",
]

for script in slurm_scripts:
    for ses in sessions:
        print(f"Submitting job for Subject: {ses}")
    
        # The command: sbatch script_name arg1 arg2
        cmd = ["sbatch", script, ses]
    
        try:
	    # Run the sbatch command
            result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
            )

            # Print the Slurm job ID (standard output from sbatch)
            print(f"Success: {result.stdout.strip()}")
	
        except subprocess.CalledProcessError as e:
            print(f"Error submitting job: {e.stderr}")

        # Optional: slight delay to be gentle on the scheduler
        time.sleep(0.5)

print("\nAll jobs submitted.")
