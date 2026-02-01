from glob import glob
import pandas as pd
from datetime import datetime
import sys

session = sys.argv[1]

current_datetime = datetime.now()
formatted_date = current_datetime.strftime("%Y_%m_%d_%H_%M_%S")

analysis_path = "/home/data/NDClab/analyses/thrive-theta-ddm/"
data_path = f"derivatives/behavior/{session}/fitting_2018/"
output_path = f"{analysis_path}/derivatives/behavior/{session}/"

file_list = sorted(glob(f"{analysis_path}/{data_path}/*ddm_output_data*"))

csv_list = []
for f in file_list:
    csv_list.append(pd.read_csv(f))

ddm_df = pd.concat(csv_list)
ddm_df.columns = ['sub', 'a', 'ter', 'p', 'rd', 'sda', 'fitStat', 'iterNum', 'pre_accuracy', 'soc', 'seed']
# drop garbage output files which were created due to an error when submitting incorrect lists of subjects (empty subject id)
ddm_df = ddm_df.dropna(subset="sub").sort_values(by="sub").reset_index(drop=True)
ddm_df.to_csv(f"{output_path}/ddm_fit_{session}_{formatted_date}.csv", index = False)
