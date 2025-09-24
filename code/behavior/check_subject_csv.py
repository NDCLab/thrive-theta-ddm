import os
import re
import glob

input_dataset_path = "/home/data/NDClab/datasets/thrive-dataset/"
data_path = "sourcedata/checked/"
sub_path = "s1_r1/psychopy/"

sub_folders = [i for i in os.listdir(input_dataset_path + data_path) if i.startswith("sub-")]
subjects = sorted([re.findall(r'\d+', item)[0] for item in sub_folders])
for sub in subjects:
    subject_folder = (input_dataset_path + data_path + "sub-" + sub + os.sep + sub_path)
    num_files = len(os.listdir(subject_folder))
    if (num_files != 3) and (sub not in ["3000124", "3000008", "3000014"]):
        print("sub-{} has unresolved deviation in psychopy data ({} files), skipping ...".format(sub, num_files))
        pass
    else:
        print("sub-{} checked".format(sub))    
        pattern = "{}sub-{}_arrow-alert-v1-*_psychopy_s1_r1_e1.csv".format(subject_folder, sub)
        assert len(glob.glob(pattern)) != 0, f"sub-{sub} main .csv has deviation in filename"