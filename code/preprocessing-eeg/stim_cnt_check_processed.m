% STIM_CNT_CHECK_PROCESSED Count Trials in Processed Data
%
% Usage:
%   Run this script to scan processed .set files and count valid trials
%   per condition. Saves results to CSV.

%clear all;
%clc;

% to Run on FIU HPC
% create a local cluster object
cluster = parcluster('local');

% start matlabpool with max workers set in the slurm file
parpool(cluster, str2num(getenv('SLURM_CPUS_PER_TASK'))) % this should be same as --cpus-per-task

% temp test code; remove
pool = gcp('nocreate');  % Get the current parallel pool without creating a new one
if isempty(pool)
    disp('No parallel pool is currently running.');
else
    disp(['Parallel pool with ', num2str(pool.NumWorkers), ' workers is running.']);
end

addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeg_preprocessing'));% enter the path of the folder in this line
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b')); % enter the path of the EEGLAB folder in this line
rmpath(['/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal']);

% Define the dataset path and session
%dataset_path = '/home/data/NDClab/datasets/thrive-dataset/derivatives/preprocessed/'; % Modify if your EEG data is in another folder
dataset_path = '/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/preprocessed/csd_data/'; % Modify if your EEG data is in another folder
session = 's1_r1'; % Modify if using for another session
% Get the list of subject data paths
subject_data_paths = dir(fullfile(dataset_path, 'sub-*/', session));
to_remove = ismember({subject_data_paths.name}, {'.', '..','.DS_Store'});
subject_data_paths = subject_data_paths(~to_remove);
subject_data_paths = {subject_data_paths.folder};
subject_data_paths = sort(subject_data_paths);
subject_data_paths = unique(subject_data_paths);

% Define the pattern to extract the subject ID
pattern = 'sub-(\d{7})';

% Define the stimulus events to count
% --- Log File Setup ---
datetime_str = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
log_dir = 'logs';
% Ensure qa_logs directory exists
if ~exist(log_dir, 'dir')
   mkdir(log_dir)
end
%csv_file = fullfile(log_dir, sprintf('qa_log_eeg_deviation_stim_counts_%s_%s.csv', session, datetime_str));
csv_file = fullfile(log_dir, sprintf('qa_log_eeg_processed_trial_counts_%s.csv', session));
dfile = fullfile(log_dir, sprintf('qa_log_eeg_processed_trial_counts_%s_%s.txt', session, datetime_str));
diary(dfile);

% Pre-allocate a cell array to store results from each worker
% Each cell will hold a table for one subject
results_collector = cell(length(subject_data_paths), 1);

% Loop through each subject data path
parfor i = 1:length(subject_data_paths)
    sub_path = subject_data_paths{i};
    fprintf('Processing: %s\n', sub_path);
    
    no_data = 0;

    % Extract the subject ID using the pattern
    sub_match = regexp(sub_path, pattern, 'tokens');
    if isempty(sub_match)
        fprintf('Could not extract subject ID from path: %s. Skipping.\n', sub_path);
        continue;
    end
    sub = sub_match{1}{1};

    % Define the subject folder path
    subject_folder = fullfile(dataset_path, sprintf('sub-%s/%s/eeg/', sub, session));

    % Get all files in the subject folder
    sub_files = dir(fullfile(subject_folder, '*'));
    sub_files = {sub_files.name};
    sub_files = sub_files(~ismember(sub_files, {'.', '..','.DS_Store'}));

    % Check for no-data.txt
    if any(contains(sub_files, 'no-data.txt'))
        no_data = 1;
        fprintf('sub-%s has NO DATA! Skipping.\n', sub);
        continue; % Skip this subject
    end

    % --- Subject has deviation, proceed to count stim markers ---
    % fprintf('sub-%s HAS deviation.txt. Processing EEG files...\n', sub);
    
    % Find all .vhdr files
    eeg_files = dir(fullfile(subject_folder, '*all_eeg_processed*.set'));

    if isempty(eeg_files)
        fprintf('sub-%s has NO .set files. Skipping.\n', sub);
        continue;
    end
    
    % Create a temporary table for this subject's results
    num_files = length(eeg_files);
    sub_table = table('Size', [num_files, 10], ...
                      'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
                      'VariableNames', {'Subject', 'Filename', 'StimCountNonSoc', 'StimCountSoc', 'resp_s_i_0', 'resp_s_i_1', 'resp_s_c_1', 'resp_ns_i_0', 'resp_ns_i_1', 'resp_ns_c_1'});

    % Loop through each .vhdr file
    for f = 1:num_files
        fname = eeg_files(f).name;
        fpath = eeg_files(f).folder;
        stim_count_nonsoc = 0; % Default
        stim_count_soc = 0; % Default
        
        fprintf('sub-%s: Loading file %s...\n', sub, fname);
        
        try
            % Load EEG data
            EEG = pop_loadset(fname, fpath);
            EEG = eeg_checkset(EEG);
            EEG = pop_selectevent(EEG, 'latency', '-.1 <= .1', 'deleteevents', 'on');

            stim_count_nonsoc = sum(ismember(string({EEG.event.observation}), 'ns'));
            stim_count_soc = sum(ismember(string({EEG.event.observation}), 's'));
            fprintf('sub-%s: File %s has %d NONSOC and %d SOC stimulus events.\n', sub, fname, stim_count_nonsoc, stim_count_soc);

            mask = ...
                string({EEG.event.eventType}) == 'resp' & ...
                string({EEG.event.observation}) == 's' & ...
                string({EEG.event.congruency}) == 'i' & ...
                [EEG.event.accuracy] == 0 & ...
                [EEG.event.responded] == 1  & ...
                [EEG.event.validRt] == 1 & ...
                [EEG.event.extraResponse] == 0;
            resp_s_i_0 = sum(mask);
            
            mask = ...
                string({EEG.event.eventType}) == 'resp' & ...
                string({EEG.event.observation}) == 's' & ...
                string({EEG.event.congruency}) == 'i' & ...
                [EEG.event.accuracy] == 1 & ...
                [EEG.event.responded] == 1  & ...
                [EEG.event.validRt] == 1 & ...
                [EEG.event.extraResponse] == 0;
            resp_s_i_1 = sum(mask);
            
            mask = ...
                string({EEG.event.eventType}) == 'resp' & ...
                string({EEG.event.observation}) == 's' & ...
                string({EEG.event.congruency}) == 'c' & ...
                [EEG.event.accuracy] == 1 & ...
                [EEG.event.responded] == 1  & ...
                [EEG.event.validRt] == 1 & ...
                [EEG.event.extraResponse] == 0;
            resp_s_c_1 = sum(mask);
            
            mask = ...
                string({EEG.event.eventType}) == 'resp' & ...
                string({EEG.event.observation}) == 'ns' & ...
                string({EEG.event.congruency}) == 'i' & ...
                [EEG.event.accuracy] == 0 & ...
                [EEG.event.responded] == 1  & ...
                [EEG.event.validRt] == 1 & ...
                [EEG.event.extraResponse] == 0;
            resp_ns_i_0 = sum(mask);
            
            mask = ...
                string({EEG.event.eventType}) == 'resp' & ...
                string({EEG.event.observation}) == 'ns' & ...
                string({EEG.event.congruency}) == 'i' & ...
                [EEG.event.accuracy] == 1 & ...
                [EEG.event.responded] == 1  & ...
                [EEG.event.validRt] == 1 & ...
                [EEG.event.extraResponse] == 0;
            resp_ns_i_1 = sum(mask);
            
            mask = ...
                string({EEG.event.eventType}) == 'resp' & ...
                string({EEG.event.observation}) == 'ns' & ...
                string({EEG.event.congruency}) == 'c' & ...
                [EEG.event.accuracy] == 1 & ...
                [EEG.event.responded] == 1  & ...
                [EEG.event.validRt] == 1 & ...
                [EEG.event.extraResponse] == 0;
            resp_ns_c_1 = sum(mask);

        catch ME
            fprintf('sub-%s: FAILED to load or process %s. Error: %s\n', sub, fname, ME.message);
            stim_count_nonsoc = NaN; % Use NaN to indicate error
            stim_count_soc = NaN; % Use NaN to indicate error
            resp_s_i_0 = NaN; % Use NaN to indicate error
            resp_s_i_1 = NaN; % Use NaN to indicate error
            resp_s_c_1 = NaN; % Use NaN to indicate error
            resp_ns_i_0 = NaN; % Use NaN to indicate error
            resp_ns_i_1 = NaN; % Use NaN to indicate error
            resp_ns_c_1 = NaN; % Use NaN to indicate error
        end
        
        % Add data to the subject's temporary table
        sub_table(f, :) = {string(sub), string(fname), stim_count_nonsoc, stim_count_soc, resp_s_i_0, resp_s_i_1, resp_s_c_1, resp_ns_i_0, resp_ns_i_1, resp_ns_c_1};
    end
    
    % Store this subject's table in the main collector cell
    results_collector{i} = sub_table;

end  

% --- Post-processing: Combine results and write to CSV ---

% Combine all non-empty tables from the collector
final_table = vertcat(results_collector{:});

% Write the final table to CSV
if ~isempty(final_table)
    writetable(final_table, csv_file);
    fprintf('\nSuccessfully wrote results to %s\n', csv_file);
else
    fprintf('\nNo files were found or processed. CSV not written.\n');
end

diary off
