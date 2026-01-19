
% STIM_CNT_CHECK_PROCESSED_CSD Count Trials in Processed CSD Data
%
% Usage:
%   Run this script to scan processed CSD .set files and count valid trials.
%   Saves results to CSV.

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

addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeg_preprocessing'));
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b')); 
rmpath(['/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal']);

% --- Configuration ---
% Define the dataset path and session
dataset_path = '/home/data/NDClab/analyses/thrive-theta-ddm/derivatives/preprocessed/csd_data/'; 
session = 's1_r1'; 

% Define where the files actually are
data_dir = fullfile(dataset_path, session);

% --- Subject Discovery ---
% Find all .set files in the data directory
fprintf('Scanning for files in: %s\n', data_dir);
all_files = dir(fullfile(data_dir, '*sub-*.set'));

if isempty(all_files)
    error('No .set files found in %s', data_dir);
end

% Extract unique Subject IDs from the filenames
pattern = 'sub-(\d{7})';
subject_ids = {};
for k = 1:length(all_files)
    tok = regexp(all_files(k).name, pattern, 'tokens');
    if ~isempty(tok)
        subject_ids{end+1} = tok{1}{1}; 
    end
end

% Create a sorted, unique list of Subject IDs
subject_list = unique(subject_ids);
subject_list = sort(subject_list);
fprintf('Found %d unique subjects.\n', length(subject_list));

% --- Log File Setup ---
datetime_str = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
log_dir = 'logs';
if ~exist(log_dir, 'dir')
   mkdir(log_dir)
end
csv_file = fullfile(log_dir, sprintf('qa_log_eeg_processed_trial_counts_csv_%s.csv', session));
dfile = fullfile(log_dir, sprintf('qa_log_eeg_processed_trial_counts_csd_%s_%s.txt', session, datetime_str));
diary(dfile);

% Pre-allocate results
results_collector = cell(length(subject_list), 1);

% --- Processing Loop ---
parfor i = 1:length(subject_list)
    sub = subject_list{i}; 
    fprintf('Processing Subject: %s\n', sub);
    
    % Find files for this specific subject in the data_dir
    % Pattern: *sub-ID*.set
    file_pattern = sprintf('*sub-%s*.set', sub);
    eeg_files = dir(fullfile(data_dir, file_pattern));
    
    if isempty(eeg_files)
        fprintf('sub-%s: No .set files found (unexpected).\n', sub);
        continue;
    end
    
    num_files = length(eeg_files);
    
    % Prepare table for this subject
    sub_table = table('Size', [num_files, 10], ...
                      'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
                      'VariableNames', {'Subject', 'Filename', 'StimCountNonSoc', 'StimCountSoc', 'resp_s_i_0', 'resp_s_i_1', 'resp_s_c_1', 'resp_ns_i_0', 'resp_ns_i_1', 'resp_ns_c_1'});

    for f = 1:num_files
        fname = eeg_files(f).name;
        fpath = eeg_files(f).folder;
        
        % Initialize counts as NaN (in case of error) or 0
        stim_count_nonsoc = 0; 
        stim_count_soc = 0; 
        resp_s_i_0 = 0; resp_s_i_1 = 0; resp_s_c_1 = 0; 
        resp_ns_i_0 = 0; resp_ns_i_1 = 0; resp_ns_c_1 = 0;

        fprintf('sub-%s: Loading file %s...\n', sub, fname);
        
        try
            % Load EEG data
            EEG = pop_loadset(fname, fpath);
            EEG = eeg_checkset(EEG);
            EEG = pop_selectevent(EEG, 'latency', '-.1 <= .1', 'deleteevents', 'on');

            % Count Stim Events
            stim_count_nonsoc = sum(ismember(string({EEG.event.observation}), 'ns'));
            stim_count_soc = sum(ismember(string({EEG.event.observation}), 's'));
            
            % Count Response Events (using logical masks)
            % s_i_0
            mask = string({EEG.event.eventType}) == 'resp' & string({EEG.event.observation}) == 's' & string({EEG.event.congruency}) == 'i' & [EEG.event.accuracy] == 0 & [EEG.event.responded] == 1  & [EEG.event.validRt] == 1 & [EEG.event.extraResponse] == 0;
            resp_s_i_0 = sum(mask);
            
            % s_i_1
            mask = string({EEG.event.eventType}) == 'resp' & string({EEG.event.observation}) == 's' & string({EEG.event.congruency}) == 'i' & [EEG.event.accuracy] == 1 & [EEG.event.responded] == 1  & [EEG.event.validRt] == 1 & [EEG.event.extraResponse] == 0;
            resp_s_i_1 = sum(mask);
            
            % s_c_1
            mask = string({EEG.event.eventType}) == 'resp' & string({EEG.event.observation}) == 's' & string({EEG.event.congruency}) == 'c' & [EEG.event.accuracy] == 1 & [EEG.event.responded] == 1  & [EEG.event.validRt] == 1 & [EEG.event.extraResponse] == 0;
            resp_s_c_1 = sum(mask);
            
            % ns_i_0
            mask = string({EEG.event.eventType}) == 'resp' & string({EEG.event.observation}) == 'ns' & string({EEG.event.congruency}) == 'i' & [EEG.event.accuracy] == 0 & [EEG.event.responded] == 1  & [EEG.event.validRt] == 1 & [EEG.event.extraResponse] == 0;
            resp_ns_i_0 = sum(mask);
            
            % ns_i_1
            mask = string({EEG.event.eventType}) == 'resp' & string({EEG.event.observation}) == 'ns' & string({EEG.event.congruency}) == 'i' & [EEG.event.accuracy] == 1 & [EEG.event.responded] == 1  & [EEG.event.validRt] == 1 & [EEG.event.extraResponse] == 0;
            resp_ns_i_1 = sum(mask);
            
            % ns_c_1
            mask = string({EEG.event.eventType}) == 'resp' & string({EEG.event.observation}) == 'ns' & string({EEG.event.congruency}) == 'c' & [EEG.event.accuracy] == 1 & [EEG.event.responded] == 1  & [EEG.event.validRt] == 1 & [EEG.event.extraResponse] == 0;
            resp_ns_c_1 = sum(mask);

            fprintf('sub-%s: Processed %s successfully.\n', sub, fname);

        catch ME
            fprintf('sub-%s: FAILED to process %s. Error: %s\n', sub, fname, ME.message);
            % Set to NaN on failure
            stim_count_nonsoc = NaN; stim_count_soc = NaN; 
            resp_s_i_0 = NaN; resp_s_i_1 = NaN; resp_s_c_1 = NaN; 
            resp_ns_i_0 = NaN; resp_ns_i_1 = NaN; resp_ns_c_1 = NaN;
        end
        
        sub_table(f, :) = {string(sub), string(fname), stim_count_nonsoc, stim_count_soc, resp_s_i_0, resp_s_i_1, resp_s_c_1, resp_ns_i_0, resp_ns_i_1, resp_ns_c_1};
    end
    
    results_collector{i} = sub_table;
end  

% --- Post-processing ---
final_table = vertcat(results_collector{:});

if ~isempty(final_table)
    writetable(final_table, csv_file);
    fprintf('\nSuccessfully wrote results to %s\n', csv_file);
else
    fprintf('\nNo files were processed. CSV not written.\n');
end

diary off
