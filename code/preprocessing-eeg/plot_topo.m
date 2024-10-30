clear % clear matlab workspace
clc % clear matlab command window

erp=true;
csd=false;

%running in "EEG_training" folder on your computer
main_dir = '/Users/fzaki001/Downloads/EEG_training';

%Location of MADE and ADJUSTED-ADJUST scripts
addpath(genpath([main_dir filesep 'MADE-EEG-preprocessing-pipeline']));% enter the path of the EEGLAB folder in this line
% addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeg_preprocessing'));% enter the path of the folder in this line

%Location of "EEG
addpath(genpath([main_dir filesep 'eeglab13_4_4b']));% enter the path of the EEGLAB folder in this line
% addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b'));% enter the path of the EEGLAB folder in this line

%remove path to octave functions inside matlab to prevent errors when
rmpath([main_dir filesep 'eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])
% rmpath(['/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])

%location of analysis folder
analysis_dir = '/Users/fzaki001/thrive-theta-ddm';
% analysis_dir = '/home/data/NDClab/analyses/thrive-theta-ddm';

%location of dataset folder
dataset_dir = '/Users/fzaki001/thrive-theta-ddm';
% dataset_dir = '/home/data/NDClab/datasets/thrive-dataset/derivatives/preprocessed/';
summary_csv_path = '/Users/fzaki001/thrive-theta-ddm/derivatives/behavior/summary.csv';
% summary_csv_path = '/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/behavior/summary-eeg.csv';

% 1. Enter the path of the folder that has the data to be analyzed
data_location = [dataset_dir filesep 'derivatives' filesep 'preprocessed'];

% 2. Enter the path of the folder where you want to save the postprocessing outputs
output_location = [analysis_dir filesep 'derivatives' filesep 'preprocessed/erp_check/topo/'];

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

sub_ind_soc = [1, 2, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, ...
    22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, ...
    37, 38, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56, 57, 58, ...
    59, 60, 61, 62, 63, 65, 67, 68, 69, 70, 71, 72, 73, 74, ...
    78, 79, 80, 81, 82, 83, 84, 85, 88, 89, 90, 91, 92, 93, 94, 95, 96, 98, 99, ...
    100, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, ...
    117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, ...
    133, 134, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 147, 148, 149, 151];
sub_ind_nonsoc = [1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, ...
    23, 24, 25, 26, 27, 29, 30, 32, 33, 34, 35, 36, ...
    37, 38, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56, 57, 58, ...
    59, 60, 61, 62, 63, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, ...
    78, 79, 80, 81, 82, 83, 84, 85, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, ...
    100, 101, 102, 103, 104, 105, 106, 107, 109, 110, 111, 112, 113, 115, 116, ...
    117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, ...
    133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 147, 148, 149, 151];

% path_to_mat = "/Users/fzaki001/Downloads/thrive_Resp_erps_csd_min_6t_60acc.mat"
path_to_mat = "/Users/fzaki001/thrive-theta-ddm/derivatives/preprocessed/erp_check/thrive_Resp_erps_min_6t_60acc.mat"
%load the mat file that has the erps and subject list
mat_arr = load(path_to_mat)

%make a copy/rename the erp matrix
allData = mat_arr.erpDat_data;

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

%start and end time range for component of interest
compStartTime = 0; %(in ms)
compEndTime = 100 ; %(in ms)

%find closest values in (rounded) EEG.times to the specified start/stop
[temp,compStartIdx] = min(abs(EEG.times-compStartTime));
[temp2,compEndIdx] = min(abs(EEG.times-compEndTime));

%idxs of time range to plot topo for
compRange = compStartIdx:compEndIdx;

if erp
    vmin = -4;
    vmax = 4;
elseif csd
    vmin = -10;
    vmax = 10;
end

%% Plot topos!!

%pull out conditions of interest for all subs, and average over time
%range of interest
s_resp_incon_error = mean(newData(sub_ind_soc,1,:,compRange),4);
s_resp_incon_corr = mean(newData(sub_ind_soc,2,:,compRange),4);
ns_resp_incon_error = mean(newData(sub_ind_nonsoc,3,:,compRange),4);
ns_resp_incon_corr = mean(newData(sub_ind_nonsoc,4,:,compRange),4);

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
fontsize = 16
figure('Visible', 'off')
subplot(1, 3, 1);
topoplot(ns_resp_incon_corr_Mean_diff, EEG.chanlocs, 'maplimits', [vmin vmax], 'electrodes', 'on', 'gridscale', 300, 'plotrad', .6)
set(get(gca, 'title'), 'String', 'Alone Error Minus Correct (0-100 ms)', 'FontSize', fontsize);

subplot(1, 3, 2);
topoplot(s_resp_incon_corr_Mean_diff, EEG.chanlocs, 'maplimits', [vmin vmax], 'electrodes', 'on', 'gridscale', 300, 'plotrad', .6)
set(get(gca, 'title'), 'String', 'Social Error Minus Correct (0-100 ms)', 'FontSize', fontsize);

subplot(1, 3, 3);
topoplot(sDiff_resp_incon_error_Mean_diff, EEG.chanlocs, 'maplimits', [vmin*0.25 vmax*0.25], 'electrodes', 'on', 'gridscale', 300, 'plotrad', .6)
set(get(gca, 'title'), 'String', 'Social Minus Alone Error (0-100 ms)', 'FontSize', fontsize);

set(gcf, 'Position', [0 0 1440 900]);
saveas(gcf, strcat(output_location, 'ern_avg_topo.png'));

%% Plot individual subject topos!!

%pull out conditions of interest for all subs, and average over time
%range of interest
all_s_resp_incon_error = mean(newData(:,1,:,compRange),4);
all_s_resp_incon_corr = mean(newData(:,2,:,compRange),4);
all_ns_resp_incon_error = mean(newData(:,3,:,compRange),4);
all_ns_resp_incon_corr = mean(newData(:,4,:,compRange),4);

for s = 1:length(mat_arr.sub_id)

    %average across subs
    s_resp_incon_error_Mean = squeeze(all_s_resp_incon_error(s,:,:));
    s_resp_incon_corr_Mean = squeeze(all_s_resp_incon_corr(s,:,:));
    ns_resp_incon_error_Mean = squeeze(all_ns_resp_incon_error(s,:,:));
    ns_resp_incon_corr_Mean = squeeze(all_ns_resp_incon_corr(s,:,:));

    %compute difference topo
    s_resp_incon_corr_Mean_diff = s_resp_incon_error_Mean - s_resp_incon_corr_Mean;
    ns_resp_incon_corr_Mean_diff = ns_resp_incon_error_Mean - ns_resp_incon_corr_Mean;
    sDiff_resp_incon_error_Mean_diff = s_resp_incon_error_Mean - ns_resp_incon_error_Mean;

    figure('Visible', 'off')
    subplot(1, 3, 1);
    topoplot(ns_resp_incon_corr_Mean_diff, EEG.chanlocs, 'maplimits', [vmin vmax], 'electrodes', 'on', 'gridscale', 50, 'plotrad', .6)
    title([mat_arr.sub_id{s} ' Alone Error Minus Correct (0-100 ms)'], 'FontSize', fontsize);
    
    subplot(1, 3, 2);
    topoplot(s_resp_incon_corr_Mean_diff, EEG.chanlocs, 'maplimits', [vmin vmax], 'electrodes', 'on', 'gridscale', 50, 'plotrad', .6)
    title([mat_arr.sub_id{s} ' Social Error Minus Correct (0-100 ms)'], 'FontSize', fontsize);
    
    subplot(1, 3, 3);
    topoplot(sDiff_resp_incon_error_Mean_diff, EEG.chanlocs, 'maplimits', [vmin*0.25 vmax*0.25], 'electrodes', 'on', 'gridscale', 50, 'plotrad', .6)
    title([mat_arr.sub_id{s} ' Social Minus Alone Error (0-100 ms)'], 'FontSize', fontsize);
    
    set(gcf, 'Position', [0 0 1440 900]);
    saveas(gcf, strcat(output_location, mat_arr.sub_id{s},'.png'));
end