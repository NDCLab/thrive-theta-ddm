from glob import glob
import pandas as pd
from datetime import datetime
import sys

session = sys.argv[1]
current_datetime = datetime.now()
formatted_date = current_datetime.strftime("%Y_%m_%d_%H_%M_%S")

dataset_path = "/home/data/NDClab/analyses/thrive-theta-ddm/"
data_path = f"derivatives/behavior/{session}/fitting_2018/"
output_path = f"{dataset_path}/code/ddm/"

file_list = sorted(glob(f"{dataset_path}/{data_path}/*ddm_output_data*"))

csv_list = []
if len(file_list) != 0:
    for f in file_list:
        csv_list.append(pd.read_csv(f))

    ddm_df = pd.concat(csv_list)
    ddm_df.columns = ['sub', 'a', 'ter', 'p', 'rd', 'sda', 'fitStat', 'iterNum', 'pre_accuracy', 'soc', 'seed']
    ddm_df = ddm_df.dropna(subset="sub").sort_values(by="sub").reset_index(drop=True)
    pd.Series(ddm_df["sub"].unique(), name = "sub").to_csv(f"{output_path}/fitted_subjects_{session}_{formatted_date}.csv", index = False)

else:
    pd.Series(name = "sub").to_csv(f"{output_path}/fitted_subjects_{session}_{formatted_date}.csv", index = False)
    print("No subjects with DDM fit were found")
