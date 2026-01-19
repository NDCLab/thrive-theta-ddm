import subprocess
import time

"""
Script to submit batch jobs for ERP processing and plotting.

This script iterates over defined sessions and conditions (Laplacian vs standard)
and submits SLURM jobs for computing ERP means and generating ERP plots.

Usage:
    python run_erp_processing.py
"""

# Define your iteration parameters
sessions = [
"s1_r1",
#"s2_r1",
#"s3_r1"
]
laplacian = ["0", "1"]

slurm_scripts = [
"compute_erp_means.sub",
"/home/data/NDClab/analyses/thrive-theta-ddm/code/figures/plot_erp.sub",
]

for script in slurm_scripts:
    for ses in sessions:
        for cond in laplacian:
            print(f"Submitting job for Subject: {ses}, Condition: {cond}")
        
            # The command: sbatch script_name arg1 arg2
            cmd = ["sbatch", script, ses, cond]
        
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
