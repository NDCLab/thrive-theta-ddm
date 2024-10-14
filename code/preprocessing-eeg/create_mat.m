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
summary_csv_path = '/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/behavior/summary.csv';

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

%% pull resp-locked erp mat file

%read in behavioral data for participants
behavior_info = readtable(summary_csv_path);

%specify min number of trials per condition (if file contains less than
%this number for ANY condition, then they will be skipped for ALL conditions
minTrials = 8;

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

    %if participant has low accuracy in either condition, skip that
    %participant for ALL conditions
    % if (behavior_info{behavior_id_match_idxs,'acc_nonsoc'} < acc_cutoff || behavior_info{behavior_id_match_idxs,'acc_soc'} < acc_cutoff)
    %     continue
    % end
    % 
    % if (behavior_info{behavior_id_match_idxs,'x6_or_more_err_nonsoc'} < 6 || behavior_info{behavior_id_match_idxs,'x6_or_more_err_soc'} < 6)
    %     continue
    % end
    %load the original data set
    EEG = pop_loadset( 'filename', datafile_names{subject}, 'filepath', datafile_paths{subject});
    EEG = eeg_checkset( EEG );

    %remove all the non-stim-locking markers (should have done already...)
    EEG = pop_selectevent( EEG, 'latency','-.1 <= .1','deleteevents','on');
    EEG = eeg_checkset( EEG );

    % NOTE %
    %the logic of checking conditions and then looping over conditions
    %below is fairly hard-coded and could be be substantially improved to
    %allow for easier reuse when number of conditions or number of
    %variables per condition changes.
    %
    %before pulling trials of interest, for any conditions, check to make
    %sure this file/participant has more than minTrials for EACH condition.
    %If the file/participant is below minTrials for even one of the
    %conditions that will be pulled, then the file/participant is skipped
    %entirely and no condition data at all will be pulled for this specific
    %file (but the participant can still have data pulled for another one
    %of their files from another visit).
    %
    %count trials for each condition of interest and store in numTrials vector
    numTrials(1) = length(find( (strcmp({EEG.event.observation}, "s")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 0) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    numTrials(2) = length(find( (strcmp({EEG.event.observation}, "s")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    numTrials(3) = length(find( (strcmp({EEG.event.observation}, "ns")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 0) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    numTrials(4) = length(find( (strcmp({EEG.event.observation}, "ns")) & (strcmp({EEG.event.eventType}, "resp")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));

    %logical test if the number of trials for each condition (numTrials vector)
    %are NOTE all >= minTrials. If statement is true, then participant/file
    %is skipped and for loop over files continues to next file
    % if ~(sum(numTrials >= minTrials) == length(numTrials))
    %     continue
    % end

    % loop through conditions of interest for this file (combo of event types)
    %
    % specify number of conditions using a seperate conditionNums var, so
    % that it can be referenced below when iterating idx counters (to only
    %iterate when c == length(conditionNums);
    conditionNums = 1:4;
    %
    for c = conditionNums

        if (c==1) % social error
            observation = 's';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 0;
            responded = 1;
            validRt = 1;
        elseif (c==2) % social correct
            observation = 's';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 1;
            responded = 1;
            validRt = 1;
        elseif (c==3) % nonsocial error
            observation = 'ns';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 0;
            responded = 1;
            validRt = 1;
        elseif (c==4) % nonsocial correct
            observation = 'ns';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 1;
            responded = 1;
            validRt = 1;
        end

        %select combintion of event types of interest based on vars above
        EEG1 = pop_selectevent( EEG, 'latency','-1<=1','observation',observation,'eventType',eventType,'congruency',congruency,'accuracy',accuracy,'responded',responded,'validRt',validRt,'deleteevents','on','deleteepochs','on','invertepochs','off');
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
save('thrive_Resp_erps_min_8t_60acc.mat','erpDat_data', 'erpDat_subIds')