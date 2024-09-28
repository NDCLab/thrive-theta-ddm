%
% Modified on 2024/04/11 to process thrive dataset flanker data
%
% This script was created by George Buzzell for the NDC Lab EEG Training
% Workshop on 02/22. This script uses parts of the "set up" structure from
% the MADE preprocessing pipeline (Debnath, Buzzell, et. al., 2020)

clear % clear matlab workspace
clc % clear matlab command window

%running in "EEG_training" folder on your computer
main_dir = '/Users/fzaki001/Downloads/EEG_training';

%% Setting up other things

%Location of MADE and ADJUSTED-ADJUST scripts
addpath(genpath([main_dir filesep 'MADE-EEG-preprocessing-pipeline']));% enter the path of the EEGLAB folder in this line

%Location of "EEG
addpath(genpath([main_dir filesep 'eeglab13_4_4b']));% enter the path of the EEGLAB folder in this line

%remove path to octave functions inside matlab to prevent errors when
rmpath([main_dir filesep 'eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])


%% setup; run this section before any other section below

%location of analysis folder
% analysis_dir = '/home/NDClab/analyses/thrive-theta-ddm/';
analysis_dir = '/Users/fzaki001/thrive-theta-ddm';

%location of dataset folder
% dataset_dir = '/home/NDClab/analyses/thrive-theta-ddm/';
dataset_dir = '/Users/fzaki001/thrive-theta-ddm';
summary_csv_path = '/Users/fzaki001/thrive-theta-ddm/derivatives/behavior/summary.csv';

% Setting up other things

% %Location of "EEG
% addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b'));% enter the path of the EEGLAB folder in this line
%
% %remove path to octave functions inside matlab to prevent errors when
% rmpath(['/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])

% 1. Enter the path of the folder that has the data to be analyzed
data_location = [dataset_dir filesep 'derivatives' filesep 'preprocessed'];

% 2. Enter the path of the folder where you want to save the postprocessing outputs
output_location = [analysis_dir filesep 'derivatives' filesep 'preprocessed/erp_check'];

% 3. this is the correct channel location file BUT INCORRECT PATH!
% channel_locations = loadbvef(strcat(dataset_dir, '/code/eeg_preprocessing/chan_locs_files/electrode_locs_files/CACS-128-X7-FIXED-64only.bvef'));
% channel_locations = loadbvef('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeg_preprocessing/chan_locs_files/electrode_locs_files/CACS-128-X7-FIXED-64only.bvef');
channel_locations = loadbvef('/Users/fzaki001/Downloads/thrive/code/CACS-128-X7-FIXED-64only.bvef');

% % 4. THESE ARE INCORRECT Markers
% stimulus_markers = {'11', '12', '21', '22'};
% respose_markers = {'111', '112', '121', '122','211', '212', '221', '222'};

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
    writetable(counts_table, [output_location filesep 'thrive_trialCounts_RespAndStim', date, '.csv'], "WriteMode", "append");

    % %format trial counts into output vector
    % output = [subIdNum, s_resp_incon_error, s_resp_incon_corr, ns_resp_incon_error, ns_resp_incon_corr];
    % %write to disk the trial counts for this participant
    % dlmwrite(strcat('thrive_trialCounts_respOnly', date, '.csv'), output, 'delimiter', ',', '-append');

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
    if (behavior_info{behavior_id_match_idxs,'acc_nonsoc'} < acc_cutoff || behavior_info{behavior_id_match_idxs,'acc_soc'} < acc_cutoff)
        continue
    end

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
    if ~(sum(numTrials >= minTrials) == length(numTrials))
        continue
    end

    % loop through conditions of interest for this file (combo of event types)
    %
    % specify number of conditions using a seperate conditionNums var, so
    % that it can be referenced below when iterating idx counters (to only
    %iterate when c == length(conditionNums);
    conditionNums = 1:4;
    %
    for c = conditionNums

        if (c==1)
            observation = 's';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 0;
            responded = 1;
            validRt = 1;
        elseif (c==2)
            observation = 's';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 1;
            responded = 1;
            validRt = 1;
        elseif (c==3)
            observation = 'ns';
            eventType = 'resp';
            congruency = 'i';
            accuracy = 0;
            responded = 1;
            validRt = 1;
        elseif (c==4)
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

%% pull STIM-locked erp mat file

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
    if (behavior_info{behavior_id_match_idxs,'acc_nonsoc'} < acc_cutoff || behavior_info{behavior_id_match_idxs,'acc_soc'} < acc_cutoff)
        continue
    end

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
    numTrials(1) = length(find( (strcmp({EEG.event.observation}, "s")) & (strcmp({EEG.event.eventType}, "stim")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    numTrials(2) = length(find( (strcmp({EEG.event.observation}, "s")) & (strcmp({EEG.event.eventType}, "stim")) & (strcmp({EEG.event.congruency}, "c")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    numTrials(3) = length(find( (strcmp({EEG.event.observation}, "ns")) & (strcmp({EEG.event.eventType}, "stim")) & (strcmp({EEG.event.congruency}, "i")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));
    numTrials(4) = length(find( (strcmp({EEG.event.observation}, "ns")) & (strcmp({EEG.event.eventType}, "stim")) & (strcmp({EEG.event.congruency}, "c")) & ([EEG.event.accuracy] == 1) & ([EEG.event.responded] == 1) & ([EEG.event.validRt] == 1)   ));


    %logical test if the number of trials for each condition (numTrials vector)
    %are NOTE all >= minTrials. If statement is true, then participant/file
    %is skipped and for loop over files continues to next file
    if ~(sum(numTrials >= minTrials) == length(numTrials))
        continue
    end

    % loop through conditions of interest for this file (combo of event types)
    %
    % specify number of conditions using a seperate conditionNums var, so
    % that it can be referenced below when iterating idx counters (to only
    %iterate when c == length(conditionNums);
    conditionNums = 1:4;
    %
    for c = conditionNums

        if (c==1)
            observation = 's';
            eventType = 'stim';
            congruency = 'i';
            accuracy = 1;
            responded = 1;
            validRt = 1;
        elseif (c==2)
            observation = 's';
            eventType = 'stim';
            congruency = 'c';
            accuracy = 1;
            responded = 1;
            validRt = 1;
        elseif (c==3)
            observation = 'ns';
            eventType = 'stim';
            congruency = 'i';
            accuracy = 1;
            responded = 1;
            validRt = 1;
        elseif (c==4)
            observation = 'ns';
            eventType = 'stim';
            congruency = 'c';
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
save('thrive_Stim_erps_min_8t_60acc.mat','erpDat_data', 'erpDat_subIds')

%% Plot ERPs!!

%load the mat file that has the erps and subject list
load('thrive_Resp_erps_min_8t_60acc.mat')
%load('thrive_Stim_erps_min_8t_60acc.mat')

%make a copy/rename the erp matrix
allData = erpDat_data;

%load in one of the participants EEGLAB-formatted data; this is to load
%parameters needed for plotting (sampling rate, chanlocs, etc).
EEG = pop_loadset( 'filename', datafile_names{1}, 'filepath', datafile_paths{1});
EEG = eeg_checkset(EEG);
eeglab redraw

%round EEG.times to nearest whole ms to make easier to work with
EEG.times = round(EEG.times);

%setup for baseline correcting the ERP data (always done before plotting or extracting
%erps, not done to the data previously to allow use of different baselines
%as a function of review comments)
startTime = -200; %(in ms)
endTime = 0 ; %(in ms)

%find closest values in (rounded) EEG.times to the specified start/stop
[temp,startIdx] = min(abs(EEG.times-startTime));
[temp2,endIdx] = min(abs(EEG.times-endTime));

%baseline corrections
Range = startIdx:endIdx;
allBase = squeeze(mean(allData(:,:,:,Range),4));
allBase = mean(allData(:,:,:,Range),4);

for i=1:size(allData,4)
    newData(:,:,:,i) = allData(:,:,:,i) - allBase;
end

%select channel(s) to plot: frontocentral cluster
chan = (newData(:,:,[1, 2, 5, 37, 34],:));
chan = (newData(:,:,[53, 55],:));
chan = (newData(:,:,[17, 49, 50, 19, 18],:));
chan = (newData(:,:,[1, 33, 17],:));

chan = mean(chan,3);

%pull out four conditions of interest for all subs
s_resp_incon_error = chan(:,1,:,:);
s_resp_incon_corr = chan(:,2,:,:);
ns_resp_incon_error = chan(:,3,:,:);
ns_resp_incon_corr = chan(:,4,:,:);

%average across subs
s_resp_incon_error_Mean = squeeze(mean(s_resp_incon_error,1));
s_resp_incon_corr_Mean = squeeze(mean(s_resp_incon_corr,1));
ns_resp_incon_error_Mean = squeeze(mean(ns_resp_incon_error,1));
ns_resp_incon_corr_Mean = squeeze(mean(ns_resp_incon_corr,1));

%label for plot and define colors for plot
blue = [0  0 1];
red = [1 0 0];

%plot the two response-related erps
figure;
hold on
plot(EEG.times, s_resp_incon_error_Mean, 'color', red, 'LineWidth', 1.5);
plot(EEG.times, s_resp_incon_corr_Mean, 'color', blue, 'LineWidth', 1.5);
plot(EEG.times, ns_resp_incon_error_Mean, 'color', red, 'LineWidth', 1.5, 'LineStyle', ':');
plot(EEG.times, ns_resp_incon_corr_Mean, 'color', blue, 'LineWidth', 1.5, 'LineStyle', ':');

%title(sprintf('Fz'), 'FontSize', 30);
legendHandle = legend('Social-Error', 'Social-Correct', 'Alone-Error', 'Alone-Correct');
set(legendHandle, 'box', 'off', 'FontSize', 26);
hold off;

% set parameters
plotStartTime = -400; %(in ms)
plotEndTime = 1800 ; %(in ms)
set(gcf, 'Color', [1 1 1]);
set(gca, 'YLim', [-10 20]);
set(gca, 'XLim', [plotStartTime plotEndTime]);
set(gca, 'FontSize', 20);
set(get(gca, 'YLabel'), 'String', 'Amplitude in  \muV', 'FontSize', 26);
set(get(gca, 'XLabel'), 'String', 'Time Relative to Response (ms)', 'FontSize', 26);
set(gca, 'Box', 'off');
set(gcf, 'Position', [0 0 1440 900]);

%% Plot Individual subject ERPs!!

%load the mat file that has the erps and subject list
load('thrive_Resp_erps_min_8t_60acc.mat')
%load('thrive_Stim_erps_min_8t_60acc.mat')

%make a copy/rename the erp matrix
allData = erpDat_data;

%load in one of the participants EEGLAB-formatted data; this is to load
%parameters needed for plotting (sampling rate, chanlocs, etc).
EEG = pop_loadset( 'filename', datafile_names{1}, 'filepath', datafile_paths{1});
EEG = eeg_checkset(EEG);
eeglab redraw

%round EEG.times to nearest whole ms to make easier to work with
EEG.times = round(EEG.times);

%setup for baseline correcting the ERP data (always done before plotting or extracting
%erps, not done to the data previously to allow use of different baselines
%as a function of review comments)
startTime = -200; %(in ms)
endTime = 0 ; %(in ms)

%find closest values in (rounded) EEG.times to the specified start/stop
[temp,startIdx] = min(abs(EEG.times-startTime));
[temp2,endIdx] = min(abs(EEG.times-endTime));

%baseline corrections
Range = startIdx:endIdx;
allBase = squeeze(mean(allData(:,:,:,Range),4));
allBase = mean(allData(:,:,:,Range),4);

for i=1:size(allData,4)
    newData(:,:,:,i) = allData(:,:,:,i) - allBase;
end

%select channel(s) to plot: frontocentral cluster
chan = (newData(:,:,[1, 2, 5, 37, 34],:));
%chan = (newData(:,:,[53, 55],:));
%chan = (newData(:,:,[17, 49, 50, 19, 18],:));
%chan = (newData(:,:,[1, 33, 17],:));

chan = mean(chan,3);

%pull out four conditions of interest for all subs
s_resp_incon_error = chan(:,1,:,:);
s_resp_incon_corr = chan(:,2,:,:);
ns_resp_incon_error = chan(:,3,:,:);
ns_resp_incon_corr = chan(:,4,:,:);

all_s_resp_incon_error = chan(:,1,:,:);
all_s_resp_incon_corr = chan(:,2,:,:);
all_ns_resp_incon_error = chan(:,3,:,:);
all_ns_resp_incon_corr = chan(:,4,:,:);

for s = 1:2

    s_resp_incon_error = all_s_resp_incon_error(s,:,:,:);
    s_resp_incon_corr = all_s_resp_incon_corr(s,:,:,:);
    ns_resp_incon_error = all_ns_resp_incon_error(s,:,:,:);
    ns_resp_incon_corr = all_ns_resp_incon_corr(s,:,:,:);

    %average across subs
    s_resp_incon_error_Mean = squeeze(mean(s_resp_incon_error,1));
    s_resp_incon_corr_Mean = squeeze(mean(s_resp_incon_corr,1));
    ns_resp_incon_error_Mean = squeeze(mean(ns_resp_incon_error,1));
    ns_resp_incon_corr_Mean = squeeze(mean(ns_resp_incon_corr,1));

    %label for plot and define colors for plot
    blue = [0  0 1];
    red = [1 0 0];

    %plot the two response-related erps
    figure;
    hold on
    plot(EEG.times, s_resp_incon_error_Mean, 'color', red, 'LineWidth', 1.5);
    plot(EEG.times, s_resp_incon_corr_Mean, 'color', blue, 'LineWidth', 1.5);
    plot(EEG.times, ns_resp_incon_error_Mean, 'color', red, 'LineWidth', 1.5, 'LineStyle', ':');
    plot(EEG.times, ns_resp_incon_corr_Mean, 'color', blue, 'LineWidth', 1.5, 'LineStyle', ':');

    %title(sprintf('%d', s), 'FontSize', 30);
    title(erpDat_subIds{s}, 'FontSize', 30);

    legendHandle = legend('Social-Error', 'Social-Correct', 'Alone-Error', 'Alone-Correct');
    set(legendHandle, 'box', 'off', 'FontSize', 26);
    hold off;

    % set parameters
    plotStartTime = -400; %(in ms)
    plotEndTime = 800 ; %(in ms)
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'YLim', [-15 15]);
    set(gca, 'XLim', [plotStartTime plotEndTime]);
    set(gca, 'FontSize', 20);
    set(get(gca, 'YLabel'), 'String', 'Amplitude in  \muV', 'FontSize', 26);
    set(get(gca, 'XLabel'), 'String', 'Time Relative to Response (ms)', 'FontSize', 26);
    set(gca, 'Box', 'off');
    set(gcf, 'Position', [0 0 1440 900]);

end

%% Plot topos!!

%load the mat file that has the erps and subject list
load('thrive_Resp_erps_min_8t_60acc.mat')

%make a copy/rename the erp matrix
allData = erpDat_data;

%load in one of the participants EEGLAB-formatted data; this is to load
%parameters needed for plotting (sampling rate, chanlocs, etc).
EEG = pop_loadset( 'filename', datafile_names{1}, 'filepath', datafile_paths{1});
EEG = eeg_checkset(EEG);
eeglab redraw

%round EEG.times to nearest whole ms to make easier to work with
EEG.times = round(EEG.times);

%setup for baseline correcting the ERP data (always done before plotting or extracting
%erps, not done to the data previously to allow use of different baselines
%as a function of review comments)
startTime = -400; %(in ms)
endTime = -200 ; %(in ms)

%find closest values in (rounded) EEG.times to the specified start/stop
[temp,startIdx] = min(abs(EEG.times-startTime));
[temp2,endIdx] = min(abs(EEG.times-endTime));

%baseline corrections
Range = startIdx:endIdx;
allBase = squeeze(mean(allData(:,:,:,Range),4));
allBase = mean(allData(:,:,:,Range),4);

for i=1:size(allData,4)
    newData(:,:,:,i) = allData(:,:,:,i) - allBase;
end

%start and end time range for component of interest
compStartTime = 0; %(in ms)
compEndTime = 100 ; %(in ms)

%find closest values in (rounded) EEG.times to the specified start/stop
[temp,compStartIdx] = min(abs(EEG.times-compStartTime));
[temp2,compEndIdx] = min(abs(EEG.times-compEndTime));

%idxs of time range to plot topo for
compRange = compStartIdx:compEndIdx;

%pull out conditions of interest for all subs, and average over time
%range of interest
s_resp_incon_error = mean(newData(:,1,:,compRange),4);
s_resp_incon_corr = mean(newData(:,2,:,compRange),4);
ns_resp_incon_error = mean(newData(:,3,:,compRange),4);
ns_resp_incon_corr = mean(newData(:,4,:,compRange),4);

%average across subs
s_resp_incon_error_Mean = squeeze(mean(s_resp_incon_error,1));
s_resp_incon_corr_Mean = squeeze(mean(s_resp_incon_corr,1));
ns_resp_incon_error_Mean = squeeze(mean(ns_resp_incon_error,1));
ns_resp_incon_corr_Mean = squeeze(mean(ns_resp_incon_corr,1));

%compute difference topo
s_resp_incon_corr_Mean_diff = s_resp_incon_error_Mean - s_resp_incon_corr_Mean;
ns_resp_incon_corr_Mean_diff = ns_resp_incon_error_Mean - ns_resp_incon_corr_Mean;
sDiff_resp_incon_error_Mean_diff = s_resp_incon_error_Mean - ns_resp_incon_error_Mean;

%plot topos
figure
topoplot(ns_resp_incon_corr_Mean_diff, EEG.chanlocs, 'maplimits', [-4 4], 'electrodes', 'on', 'gridscale', 300, 'plotrad', .6)
set(get(gca, 'title'), 'String', 'Alone Error Minus Correct (0-100 ms)', 'FontSize', 20);

figure
topoplot(s_resp_incon_corr_Mean_diff, EEG.chanlocs, 'maplimits', [-4 4], 'electrodes', 'on', 'gridscale', 300, 'plotrad', .6)
set(get(gca, 'title'), 'String', 'Social Error Minus Correct (0-100 ms)', 'FontSize', 20);

figure
topoplot(sDiff_resp_incon_error_Mean_diff, EEG.chanlocs, 'maplimits', [-1 1], 'electrodes', 'on', 'gridscale', 300, 'plotrad', .6)
set(get(gca, 'title'), 'String', 'Social Minus Alone Error (0-100 ms)', 'FontSize', 20);

%% Plot individual subject topos!!

%load the mat file that has the erps and subject list
load('thrive_Resp_erps_min_8t_60acc.mat')

%make a copy/rename the erp matrix
allData = erpDat_data;

%load in one of the participants EEGLAB-formatted data; this is to load
%parameters needed for plotting (sampling rate, chanlocs, etc).
EEG = pop_loadset( 'filename', datafile_names{1}, 'filepath', datafile_paths{1});
EEG = eeg_checkset(EEG);
eeglab redraw

%round EEG.times to nearest whole ms to make easier to work with
EEG.times = round(EEG.times);

%setup for baseline correcting the ERP data (always done before plotting or extracting
%erps, not done to the data previously to allow use of different baselines
%as a function of review comments)
startTime = -400; %(in ms)
endTime = -200 ; %(in ms)

%find closest values in (rounded) EEG.times to the specified start/stop
[temp,startIdx] = min(abs(EEG.times-startTime));
[temp2,endIdx] = min(abs(EEG.times-endTime));

%baseline corrections
Range = startIdx:endIdx;
allBase = squeeze(mean(allData(:,:,:,Range),4));
allBase = mean(allData(:,:,:,Range),4);

for i=1:size(allData,4)
    newData(:,:,:,i) = allData(:,:,:,i) - allBase;
end

%start and end time range for component of interest
compStartTime = 0; %(in ms)
compEndTime = 100 ; %(in ms)

%find closest values in (rounded) EEG.times to the specified start/stop
[temp,compStartIdx] = min(abs(EEG.times-compStartTime));
[temp2,compEndIdx] = min(abs(EEG.times-compEndTime));

%idxs of time range to plot topo for
compRange = compStartIdx:compEndIdx;

%pull out conditions of interest for all subs, and average over time
%range of interest
all_s_resp_incon_error = mean(newData(:,1,:,compRange),4);
all_s_resp_incon_corr = mean(newData(:,2,:,compRange),4);
all_ns_resp_incon_error = mean(newData(:,3,:,compRange),4);
all_ns_resp_incon_corr = mean(newData(:,4,:,compRange),4);

for s = [1,2]

    %average across subs
    s_resp_incon_error_Mean = squeeze(all_s_resp_incon_error(s,:,:));
    s_resp_incon_corr_Mean = squeeze(all_s_resp_incon_corr(s,:,:));
    ns_resp_incon_error_Mean = squeeze(all_ns_resp_incon_error(s,:,:));
    ns_resp_incon_corr_Mean = squeeze(all_ns_resp_incon_corr(s,:,:));

    %compute difference topo
    s_resp_incon_corr_Mean_diff = s_resp_incon_error_Mean - s_resp_incon_corr_Mean;
    ns_resp_incon_corr_Mean_diff = ns_resp_incon_error_Mean - ns_resp_incon_corr_Mean;
    sDiff_resp_incon_error_Mean_diff = s_resp_incon_error_Mean - ns_resp_incon_error_Mean;

    figure
    % First subplot
    subplot(1, 2, 1);
    topoplot(ns_resp_incon_corr_Mean_diff, EEG.chanlocs, 'maplimits', [-4 4], 'electrodes', 'on', 'gridscale', 50, 'plotrad', .6)
    title([erpDat_subIds{s} ' - Alone Error Minus Correct (0-100 ms)'], 'FontSize', 14);
    % Second subplot
    subplot(1, 2, 2);
    topoplot(s_resp_incon_corr_Mean_diff, EEG.chanlocs, 'maplimits', [-4 4], 'electrodes', 'on', 'gridscale', 50, 'plotrad', .6)
    title([erpDat_subIds{s} ' - Social Error Minus Correct (0-100 ms)'], 'FontSize', 14);

end

%% Extract mean component amplitudes (for statistics)

%clear output
output = [];

%create variable names for count trials output and write to disk
outputHeader = {'id, ERN, CRN, ERN_min_CRN_diff, PE_error, PE_corr, PE_err_min_corr_diff '};
dlmwrite(strcat('thrive_compMeans_', date, '.csv'), outputHeader, 'delimiter', '', '-append');

%each row corresponds to a different component. within each row, electrodes
%to include in the cluster for a given component are defined
clustCell= { [1, 2, 5, 37, 34];
    [15 14]};
%clustCell= {[17 21 22]};

%each row corresponds to a different component. within each row, the start/end
%times for a given component are defined
timeCell = {[0 100]; % ERN cluster
    [300 500]}; % PE cluster
%timeCell = {[0 100]};

%load the mat file that has the erps and subject list
load('thrive_Resp_erps_min_8t_60acc.mat')

%make a copy/rename the erp matrix
allData = erpDat_data;

%load in one of the participants EEGLAB-formatted data; this is to load
%parameters needed for plotting (sampling rate, chanlocs, etc).
EEG = pop_loadset( 'filename', datafile_names{1}, 'filepath', datafile_paths{1});
EEG = eeg_checkset(EEG);
% eeglab redraw

%round EEG.times to nearest whole ms to make easier to work with
EEG.times = round(EEG.times);

%setup for baseline correcting the ERP data (always done before plotting or extracting
%erps, not done to the data previously to allow use of different baselines
%as a function of review comments)
startTime = -400; %(in ms)
endTime = -200 ; %(in ms)

%find closest values in (rounded) EEG.times to the specified start/stop
[temp,startIdx] = min(abs(EEG.times-startTime));
[temp2,endIdx] = min(abs(EEG.times-endTime));

%baseline corrections
Range = startIdx:endIdx;
allBase = squeeze(mean(allData(:,:,:,Range),4));
allBase = mean(allData(:,:,:,Range),4);

for i=1:size(allData,4)
    newData(:,:,:,i) = allData(:,:,:,i) - allBase;
end
subjects = [50, 53] % these are just ids to put to the resulting csv, they are not used for computations
%write sub numbers to ouput
output(:,1) = subjects;

%initialize index var at 2 because i=1 is the column for subject ids
i = 2;

for comp = 1:length(clustCell) %loop through component clusters

    cluster= clustCell{comp};
    times = timeCell{comp};

    %start and end time range for component of interest
    compStartTime = times(1); %(in ms)
    compEndTime = times(2) ; %(in ms)

    %find closest values in (rounded) EEG.times to the specified start/stop
    [temp,compStartIdx] = min(abs(EEG.times-compStartTime));
    [temp2,compEndIdx] = min(abs(EEG.times-compEndTime));

    %idxs of time range to plot topo for
    compRange = compStartIdx:compEndIdx;

    %pull out conditions of interest for all subs, and average over time
    %range of interest
    resp_incon_error_avgTime = mean(newData(:,1,:,compRange),4);
    resp_incon_corr_avgTime = mean(newData(:,2,:,compRange),4);

    %average cluster of interest
    resp_incon_error_avgTimeClust = mean(resp_incon_error_avgTime(:,:,cluster),3);
    resp_incon_corr_avgTimeClust = mean(resp_incon_corr_avgTime(:,:,cluster),3);

    %compute difference scores
    resp_incon_error_avgTimeClust_diff = resp_incon_error_avgTimeClust - resp_incon_corr_avgTimeClust;

    output(:,i) = resp_incon_error_avgTimeClust;
    output(:,i+1) = resp_incon_corr_avgTimeClust;
    output(:,i+2) = resp_incon_error_avgTimeClust_diff;
    i= i+3;

    %lend cluster loop
end

%write component means to disk
dlmwrite(strcat('thrive_compMeans_', date, '.csv'), output, 'delimiter', ',', '-append');