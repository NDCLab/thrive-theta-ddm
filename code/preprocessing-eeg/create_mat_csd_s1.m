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
session = 's1_r1';
analysis_dir = '/home/data/NDClab/analyses/thrive-theta-ddm';

%location of dataset folder
dataset_dir = '/home/data/NDClab/analyses/thrive-theta-ddm';
%dataset_dir = '/home/data/NDClab/datasets/thrive-dataset';
summary_csv_path = [analysis_dir filesep 'derivatives/behavior' filesep session filesep 'summary.csv'];

% Setting up other things

% 1. Enter the path of the folder that has the data to be analyzed
%data_location = [dataset_dir filesep 'derivatives' filesep 'preprocessed'];
data_location = [dataset_dir filesep 'derivatives' filesep 'preprocessed' filesep 'csd_data' session;
% 2. Enter the path of the folder where you want to save the postprocessing outputs
output_location = [analysis_dir filesep 'derivatives' filesep 'preprocessed/erp_check' filesep session];

%modifying above, to account for files named differently
%specify parameters of data to process
task = 'all';
procStage = 'processed_data';
visitDirName = 's1_r1'; %visit folder does not list "e1"
visitFileName = 's1_r1_e1'; %file names include "e1" designation

% Read files to analyses
%datafile_info=dir([data_location filesep 'sub-*' filesep visitDirName filesep 'eeg' filesep 'sub-*_' task '_eeg_*' procStage '_' visitFileName '.set']);
datafile_info=dir([data_location filesep 'sub-*_' task '_eeg_*' procStage '_' visitFileName '.set']);
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
    mkdir(output_location);
end

%% Count trials
% switch to output directory
cd(output_location);

% %create variable names for count trials output and write to disk
% outputHeader = {'id, s_resp_incon_error, s_resp_incon_corr, ns_resp_incon_error, ns_resp_incon_corr'};
% dlmwrite(strcat('thrive_trialCounts_respOnly', date, '.csv'), outputHeader, 'delimiter', '', '-append');

diary(sprintf('erp_log_%s.log', datestr(now, 'mm_dd_yyyy_HH_MM_SS')))

%% pull resp-locked erp mat file

%read in behavioral data for participants
behavior_info = readtable(summary_csv_path);

%specify min number of trials per condition (if file contains less than
%this number for ANY condition, then they will be skipped for ALL conditions
minTrials = 6;

%specify min accuracy per condition (if file contains less than
%this number for ANY condition, then they will be skipped for ALL conditions
acc_cutoff = .6;

%initialize participant counter variable (used for indexing into large mat
%file that data is saved into)
pIdx = 1;

%initialize matrices to hold erp data and corresponding sub ids
erpDat_data = [];
erpDat_subIds = [];

% loop through each participant in the study
for subject = 1:length(datafile_names)

    %initialize numTrials for this participant/file
    numTrials = [];

    % extract participant number
    subNumText = datafile_names{subject}(5:11);

    %find row in behavior file corresponding to this participant
    behavior_id_match_idxs = find(behavior_info{:,'sub'} == str2num(subNumText));

    EEG = pop_loadset( 'filename', datafile_names{subject}, 'filepath', datafile_paths{subject});
    EEG = eeg_checkset( EEG );

    %remove all the non-stim-locking markers (should have done already...)
    EEG = pop_selectevent( EEG, 'latency','-.1 <= .1','deleteevents','on');
    EEG = eeg_checkset( EEG );


    %count trials for each condition of interest and store in numTrials vector
    numTrials(1) = length(find( (strcmp({EEG.event.observation}, "s")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 0) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1) & ([EEG.event.extraResponse] == 0) ));
    numTrials(2) = length(find( (strcmp({EEG.event.observation}, "s")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1) & ([EEG.event.extraResponse] == 0) ));
    numTrials(3) = length(find( (strcmp({EEG.event.observation}, "ns")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 0) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1) & ([EEG.event.extraResponse] == 0) ));
    numTrials(4) = length(find( (strcmp({EEG.event.observation}, "ns")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1) & ([EEG.event.extraResponse] == 0) ));

    %logical test if the number of trials for each condition (numTrials vector)
    %are NOTE all >= minTrials. If statement is true, then participant/file
    %is skipped and for loop over files continues to next file
    if ~(sum(numTrials >= minTrials) == length(numTrials))
        continue
    end

    conditionNums = 1:4;
    for c = conditionNums

        if (c==1) % social error
            observation = 's';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 0;
            responded = 1;
            validRt = 1;
            extraResponse = 0;
        elseif (c==2) % social correct
            observation = 's';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 1;
            responded = 1;
            validRt = 1;
            extraResponse = 0;
        elseif (c==3) % nonsocial error
            observation = 'ns';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 0;
            responded = 1;
            validRt = 1;
            extraResponse = 0;
        elseif (c==4) % nonsocial correct
            observation = 'ns';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 1;
            responded = 1;
            validRt = 1;
            extraResponse = 0;
        end

        %select combination of event types of interest based on vars above
        EEG1 = pop_selectevent( EEG, 'latency','-1<=1','observation', observation, 'eventType', eventType, 'congruency', congruency, 'accuracy', accuracy, 'responded', responded,'validRt', validRt, 'extraResponse', extraResponse, 'deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG1 = eeg_checkset( EEG1 );

        % Average across epoch dimension
        % this all Channel ERP only needs to be computed once
        % per condition
        meanEpochs = mean(EEG1.data, 3);

        %store data for this condition in array
        erpDat_data(pIdx,c,:,:)= meanEpochs;

        %store participant number for corresponding row in erpdat
        erpDat_subIds{pIdx,1} = datafile_names{subject}(5:11);

        %iterate idx counter IMPORTANT: ONLY ITERATE COUNTER WHEN
        %ON LAST CONDITION
        if c == length(conditionNums)%if this is the last condition of condition loop
            pIdx = pIdx + 1;
        end
        %end loop through conditions
    end
    %end loop through participants
end

%save the erps and subject list
save(sprintf('thrive_Resp_erps_csd_min_6t_%s.mat', datestr(now, 'mm_dd_yyyy_HH_MM_SS')), 'erpDat_data', 'erpDat_subIds')
