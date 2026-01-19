import subprocess

"""
Script to bulk cancel SLURM jobs within a specified range of job IDs.

This script generates a list of job IDs from a start ID to an end ID and
executes the 'scancel' command for all of them.

Usage:
    python bulk_scancel.py
"""

start_id = 2648371
end_id = 2648451

# Generate a list of all job IDs as strings (range is inclusive of end_id)
# range() in Python is exclusive at the end, so we add +1
job_ids = [str(i) for i in range(start_id, end_id + 1)]

print(f"Cancelling {len(job_ids)} jobs...")

# Run a single scancel command with all IDs as arguments
# This executes: scancel 2648371 2648372 ... 2648451
try:
    subprocess.run(["scancel"] + job_ids, check=True)
    print("Success.")
except subprocess.CalledProcessError as e:
    print(f"Error occurred: {e}")
