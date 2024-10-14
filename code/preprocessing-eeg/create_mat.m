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

%Location of "EEG
% addpath(genpath([main_dir filesep 'eeglab13_4_4b']));% enter the path of the EEGLAB folder in this line
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b'));% enter the path of the EEGLAB folder in this line

%remove path to octave functions inside matlab to prevent errors when
% rmpath([main_dir filesep 'eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])
rmpath(['/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])

%% setup; run this section before any other section below

%location of analysis folder
% analysis_dir = '/Users/fzaki001/thrive-theta-ddm';
analysis_dir = '/home/data/NDClab/analyses/thrive-theta-ddm';

%location of dataset folder
% dataset_dir = '/Users/fzaki001/thrive-theta-ddm';
dataset_dir = '/home/data/NDClab/datasets/thrive-dataset/derivatives/preprocessed/';
% summary_csv_path = '/Users/fzaki001/thrive-theta-ddm/derivatives/behavior/summary.csv';
summary_csv_path = '/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/behavior/summary-eeg.csv';

% Setting up other things

% 1. Enter the path of the folder that has the data to be analyzed
data_location = [dataset_dir filesep 'derivatives' filesep 'preprocessed'];

% 2. Enter the path of the folder where you want to save the postprocessing outputs
output_location = [analysis_dir filesep 'derivatives' filesep 'preprocessed/erp_check'];

% 3. this is the correct channel location file BUT INCORRECT PATH!
% channel_locations = loadbvef('/Users/fzaki001/Downloads/thrive/code/CACS-128-X7-FIXED-64only.bvef');
channel_locations = loadbvef('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeg_preprocessing/chan_locs_files/electrode_locs_files/CACS-128-X7-FIXED-64only.bvef');

% %specify parameters of data to process

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
if exist(output_location, 'dir') == 0
    mkdir(output_location)
end

%% Count trials
subjects = [50, 53, 60, 62, 64, 66] % these are just ids to put to the resulting csv, they are not used for computations
% n_subjects = length(subjects) % this number will be used for loops because the initial subjects list may be changed based on num of trials
% switch to output directory
cd(output_location);

% %create variable names for count trials output and write to disk
% outputHeader = {'id, s_resp_incon_error, s_resp_incon_corr, ns_resp_incon_error, ns_resp_incon_corr'};
% dlmwrite(strcat('thrive_trialCounts_respOnly', date, '.csv'), outputHeader, 'delimiter', '', '-append');

for subject=1:length(datafile_names)
    EEG=[];

    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names{subject});

    %load in raw data that is alread in eeglab (.set) format)
    EEG = pop_loadset( 'filename', datafile_names{subject}, 'filepath', datafile_paths{subject});
    EEG = eeg_checkset(EEG);

    %convert subject name to number
    %note:would be nice to modify line below to not be hard-coded for finding
    %location of subject id. eg, use some combination of strtok instead
    subIdNum = str2double(datafile_names{subject}(5:11));

    %remove all the non-stim-locking markers (should have done already...)
    EEG = pop_selectevent( EEG, 'latency','-.1 <= .1','deleteevents','on');
    EEG = eeg_checkset( EEG );

    % latency
    % duration
    % channel
    % bvtime
    % bvmknum
    % visible
    % type
    % code
    % urevent
    % observation
    % eventType
    % targetDir
    % congruency
    % responded
    % accuracy
    % rt
    % validRt
    % extraResponse
    % validTrial
    % prevTargetDir
    % prevCongruency
    % prevResponded
    % prevAccuracy
    % prevRt
    % prevValidRt
    % prevExtraResponse
    % prevValidTrial
    % nextTargetDir
    % nextCongruency
    % nextResponded
    % nextAccuracy
    % nextRt
    % nextValidRt
    % nextExtraResponse
    % nextValidTrial
    % epoch

    %count how many of each event type (combination of event types) of
    %interest are present
    s_resp_incon_error = length(find( (strcmp({EEG.event.observation}, "s")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 0) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    s_resp_incon_corr = length(find( (strcmp({EEG.event.observation}, "s")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    ns_resp_incon_error = length(find( (strcmp({EEG.event.observation}, "ns")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 0) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    ns_resp_incon_corr = length(find( (strcmp({EEG.event.observation}, "ns")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    s_stim_incon_corr = length(find( (strcmp({EEG.event.observation}, "s")) & (strcmp({EEG.event.eventType}, "stim")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    s_stim_con_corr = length(find( (strcmp({EEG.event.observation}, "s")) & (strcmp({EEG.event.eventType}, "stim")) & (strcmp({EEG.event.congruency}, "c")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    ns_stim_incon_corr = length(find( (strcmp({EEG.event.observation}, "ns")) & (strcmp({EEG.event.eventType}, "stim")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    ns_stim_con_corr = length(find( (strcmp({EEG.event.observation}, "ns")) & (strcmp({EEG.event.eventType}, "stim")) & (strcmp({EEG.event.congruency}, "c")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));

    %Create the trial counts table for trial counts
    counts_table=table({datafile_names{subject}}, {s_resp_incon_error}, {s_resp_incon_corr}, {ns_resp_incon_error}, {ns_resp_incon_corr}, {s_stim_incon_corr}, {s_stim_con_corr}, {ns_stim_incon_corr}, {ns_stim_con_corr});

    %create variable names for count trials output and write to disk
    counts_table.Properties.VariableNames = {'fileName', 's_resp_incon_error', 's_resp_incon_corr', 'ns_resp_incon_error', 'ns_resp_incon_corr', 's_stim_incon_corr', 's_stim_con_corr', 'ns_stim_incon_corr', 'ns_stim_con_corr'};

    %write/append table to disk
    writetable(counts_table, [output_location filesep 'thrive_trialCounts_RespAndStim_', date, '.csv'], "WriteMode", "append");

end