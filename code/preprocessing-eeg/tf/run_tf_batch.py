import subprocess
import sys
import time

session = sys.argv[1]

# Define your parameter space
conditions = [
'resp_s_i_0',
'resp_s_i_1',
'resp_s_c_1',
'resp_ns_i_0',
'resp_ns_i_1',
'resp_ns_c_1',
]

for cond in conditions:
    # Construct the sbatch command with --export
    # ALL inherits the default environment, then we add our custom vars
    cmd = [
	"sbatch",
	f"--export=ALL,SESS={session},COND={cond}",
	"run_tf_P.sub"
    ]
    
    # Execute the command
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"Submitted: Session {session}, Condition {cond} -> {result.stdout.strip()}")
    else:
        print(f"Error submitting {sess}: {result.stderr}")
    time.sleep(4)
