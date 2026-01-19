% STIM_CNT_CHECK_CSV Count Stimulus Markers from Raw Data
%
% Usage:
%   Run this script to scan raw VHDR files (specifically those with deviations)
%   and count stimulus markers to verify data integrity. Saves results to CSV.

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
dataset_path = '/home/data/NDClab/datasets/thrive-dataset/sourcedata/checked/'; % Modify if your behavioral data is in another folder
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
stim_events = {'S 41', 'S 42', 'S 43', 'S 44', 'S 51', 'S 52', 'S 53', 'S 54'};
stim_events_nonsoc = {'S 41', 'S 42', 'S 43', 'S 44'};
stim_events_soc = {'S 51', 'S 52', 'S 53', 'S 54'}; 
% --- Log File Setup ---
datetime_str = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
log_dir = 'logs';
% Ensure qa_logs directory exists
if ~exist(log_dir, 'dir')
   mkdir(log_dir)
end
%csv_file = fullfile(log_dir, sprintf('qa_log_eeg_deviation_stim_counts_%s_%s.csv', session, datetime_str));
csv_file = fullfile(log_dir, sprintf('qa_log_eeg_deviation_stim_counts_%s.csv', session));
dfile = fullfile(log_dir, sprintf('qa_log_eeg_deviation_stim_counts_%s_%s.txt', session, datetime_str));
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

    % Check for *deviation*.txt file
    has_deviation = false;
    for k = 1:length(sub_files)
        if contains(sub_files{k}, 'deviation.txt')
            has_deviation = true;
            break;
        end
    end

    % --- CORE LOGIC: Only process subjects with a deviation file ---
    if ~has_deviation
        fprintf('sub-%s has no deviation.txt. Skipping.\n', sub);
        continue; % Skip subject as requested
    end

    % --- Subject has deviation, proceed to count stim markers ---
    fprintf('sub-%s HAS deviation.txt. Processing EEG files...\n', sub);
    
    % Find all .vhdr files
    eeg_files = dir(fullfile(subject_folder, '*.vhdr'));

    if isempty(eeg_files)
        fprintf('sub-%s has deviation.txt but NO .vhdr files. Skipping.\n', sub);
        continue;
    end
    
    % Create a temporary table for this subject's results
    num_files = length(eeg_files);
    sub_table = table('Size', [num_files, 4], ...
                      'VariableTypes', {'string', 'string', 'double', 'double'}, ...
                      'VariableNames', {'Subject', 'Filename', 'StimCountNonSoc', 'StimCountSoc'});

    % Loop through each .vhdr file
    for f = 1:num_files
        fname = eeg_files(f).name;
        fpath = eeg_files(f).folder;
        stim_count_nonsoc = 0; % Default
        stim_count_soc = 0; % Default
        
        fprintf('sub-%s: Loading file %s...\n', sub, fname);
        
        try
            % Load EEG data
            EEG = pop_loadbv(fpath, fname);
            EEG = eeg_checkset(EEG);
            
            % Select stimulus events
            % [~, selected_events_nonsoc] = pop_selectevent(EEG, 'type', stim_events_nonsoc);
            % [~, selected_events_soc] = pop_selectevent(EEG, 'type', stim_events_soc);
            
            % Get the count
            % stim_count_nonsoc = length(selected_events_nonsoc);
            stim_count_nonsoc = sum(ismember(string({EEG.event.type}), stim_events_nonsoc));
            % stim_count_soc = length(selected_events_soc);
            stim_count_soc = sum(ismember(string({EEG.event.type}), stim_events_soc));
            fprintf('sub-%s: File %s has %d NONSOC and %d SOC stimulus events.\n', sub, fname, stim_count_nonsoc, stim_count_soc);

        catch ME
            fprintf('sub-%s: FAILED to load or process %s. Error: %s\n', sub, fname, ME.message);
            stim_count_nonsoc = NaN; % Use NaN to indicate error
            stim_count_soc = NaN; % Use NaN to indicate error
        end
        
        % Add data to the subject's temporary table
        sub_table(f, :) = {string(sub), string(fname), stim_count_nonsoc, stim_count_soc};
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
    fprintf('\nNo subjects with deviation files were found or processed. CSV not written.\n');
end

diary off
