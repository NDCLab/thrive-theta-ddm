%
% Modified on 2024/04/11 to process thrive dataset flanker data
%
% This script was created by George Buzzell for the NDC Lab EEG Training
% Workshop on 02/22. This script uses parts of the "set up" structure from
% the MADE preprocessing pipeline (Debnath, Buzzell, et. al., 2020)

clear % clear matlab workspace
clc % clear matlab command window

%% Setting up other things

%Location of MADE and ADJUSTED-ADJUST scripts
% addpath(genpath([main_dir filesep 'MADE-EEG-preprocessing-pipeline']));% enter the path of the EEGLAB folder in this line
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeg_preprocessing'));% enter the path of the folder in this line

addpath(genpath('/home/data/NDClab/analyses/thrive-theta-ddm/code/matlab/'))

%Location of "EEG
% addpath(genpath([main_dir filesep 'eeglab13_4_4b']));% enter the path of the EEGLAB folder in this line
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b'));% enter the path of the EEGLAB folder in this line

%remove path to octave functions inside matlab to prevent errors when
% rmpath([main_dir filesep 'eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])
rmpath(['/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])

%% setup; run this section before any other section below

%location of analysis folder
analysis_dir = '/home/data/NDClab/analyses/thrive-theta-ddm';

%location of dataset folder
dataset_dir = '/home/data/NDClab/datasets/thrive-dataset';

% Setting up other things

% 1. Enter the path of the folder that has the data to be analyzed
data_location = [dataset_dir filesep 'derivatives' filesep 'preprocessed'];

% 2. Enter the path of the folder where you want to save the postprocessing outputs
output_location = [analysis_dir filesep 'derivatives' filesep 'preprocessed/erp_check'];
save_location = [analysis_dir filesep 'derivatives' filesep 'preprocessed/csd_data'];
% 3. this is the correct channel location file BUT INCORRECT PATH!

%modifying above, to account for files named differently
%specify parameters of data to process
task = 'all';
procStage = 'processed_data';
visitDirName = 's1_r1'; %visit folder does not list "e1"
visitFileName = 's1_r1_e1'; %file names include "e1" designation

% Read files to analyses
datafile_info=dir([data_location filesep 'sub-*' filesep visitDirName filesep 'eeg' filesep 'sub-*_' task '_eeg_*' procStage '_' visitFileName '.set']);
datafile_info=datafile_info(~ismember({datafile_info.name},{'.', '..', '.DS_Store'}));
datafile_names={datafile_info.name};
datafile_paths={datafile_info.folder};
[filepath,name,ext] = fileparts(char(datafile_names{1}));

% Check whether EEGLAB and all necessary plugins are in Matlab path.
if exist('eeglab','file')==0
    error(['Please make sure EEGLAB is on your Matlab path. Please see EEGLAB' ...
        'wiki page for download and instalation instructions']);
end

% Create output folders to save data
if exist(save_location, 'dir') == 0
    mkdir(save_location);
end

for site = 1:64;
    trodes{site} = num2str(site);
end
Montage_64=ExtractMontage('/home/data/NDClab/analyses/thrive-theta-ddm/code/preprocessing-eeg/64ch_bv_montage_csd.csd', trodes');
% MapMontage(Montage_64);
[G, H] = GetGH(Montage_64);

diary(sprintf('csd_log_%s.log', datestr(now, 'mm-dd-yyyy_HH_MM_SS')))

% loop through each participant in the study
for subject = 1:length(datafile_names)

    % extract participant number
    subNumText = datafile_names{subject}(5:11);

    %load the original data set
    EEG = pop_loadset('filename', datafile_names{subject}, 'filepath', datafile_paths{subject});
    EEG = eeg_checkset(EEG);
    disp(datafile_names{subject})
    disp(save_location)
    %remove all the non-stim-locking markers (should have done already...)
    EEG = pop_selectevent(EEG, 'latency','-.1 <= .1','deleteevents','on');
    EEG = eeg_checkset(EEG);

    for ne = 1:length(EEG.epoch)
        myEEG = single(EEG.data(:, :, ne));
        MyResults = CSD(myEEG, G, H);            % compute CSD for <channels-by-samples> 2-D epoch
        data(:, :, ne) = MyResults;
    end
    EEG.data = data;

    data(:,:,:) = NaN;

    EEG = pop_editset(EEG, 'setname', datafile_names{subject});
    EEG = pop_saveset(EEG, 'filename', datafile_names{subject}, 'filepath', save_location);
end
