function [EEG] =  load_and_merge(path, force_cut_list)
% LOAD_AND_MERGE Load and merge raw EEG files.
%
% Usage:
%   EEG = load_and_merge(path, force_cut_list)
%
% Inputs:
%   path           - Directory containing .vhdr/.vmrk files.
%   force_cut_list - Cell array of conditions to forcibly remove ('social', 'nonsocial').
%
% Outputs:
%   EEG            - Merged EEGLAB data structure.

    if nargin < 2
        force_cut_list = {};
    end
    stim_events_nonsoc = {'S 41', 'S 42', 'S 43', 'S 44'};
    stim_events_soc = {'S 51', 'S 52', 'S 53', 'S 54'};
    stim_count_thresh = 360;
    max_stim_count = 400;
    % Step 1: Get all .vmrk files in the specified path
    vmrk_files = dir(fullfile(path, '*.vmrk'));

    % Step 2: Extract dates from each .vmrk file
    dates = cell(1, length(vmrk_files));
    for f = 1:length(vmrk_files)
	% Read the .vmrk file
	fid = fopen(fullfile(path, vmrk_files(f).name), 'r');
	vmrk_text = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
	fclose(fid);

	% Find the line starting with "Mk1=New Segment"
	date_line = '';
	for i = 1:length(vmrk_text{1})
	    line = vmrk_text{1}{i};
	    if startsWith(line, 'Mk1=New Segment')
		date_line = line;
		break;
	    end
	end

	% Extract the date string (last comma-separated value)
	date_str = '';
	if ~isempty(date_line)
	    parts = strsplit(date_line, ',');
	    date_str = parts{end};
	end

	% Convert the date string to a datetime object
	date = datetime(date_str, 'InputFormat', 'yyyyMMddHHmmssSSS');
	dates{f} = date;
    end

    % Step 3: Sort the dates and reorder the files accordingly
    dates = cellfun(@datenum, dates);
    [sorted_dates, sort_idx] = sort(dates);
    sorted_vmrk_files = vmrk_files(sort_idx);

    % Step 4: Load and merge the corresponding .vhdr files in chronological order
    EEG = [];
    for f = 1:length(sorted_vmrk_files)
	% Get the corresponding .vhdr file name
	[~, vmrk_name] = fileparts(sorted_vmrk_files(f).name);
	vhdr_file = fullfile(path, [vmrk_name '.vhdr']);

	% Load the EEG data
	EEG_temp = pop_loadbv(path, [vmrk_name '.vhdr']);

	% Merge with the existing EEG structure
	if isempty(EEG)
	    EEG = EEG_temp;
	else
	    EEG = pop_mergeset(EEG, EEG_temp);
	end
    end
    EEG = eeg_checkset(EEG);
    
    % Step 5: Decide whether to remove blocks / conditions based on stimulus marker counts
    % New helper variables for the cutting logic
    buffer_s = 5;                % Keep 5s of data before/after the cut to save EEG
    block_break_s = 6;           % A gap > 5s defines a "Block Break"
    trials_per_block = 40;

    % Calculate marker counts
    stim_count_soc = sum(ismember(string({EEG.event.type}), stim_events_soc));
    stim_count_nonsoc = sum(ismember(string({EEG.event.type}), stim_events_nonsoc));

    % Get indices in the event structure to get latencies
    soc_mask = ismember(string({EEG.event.type}), stim_events_soc);
    nonsoc_mask = ismember(string({EEG.event.type}), stim_events_nonsoc);
    soc_indices = find(soc_mask);
    nonsoc_indices = find(nonsoc_mask);

    % Initialize the rejection matrix: [StartSample, EndSample]
    regions_to_cut = [];

    % --- STEP 2: Evaluate SOCIAL Condition ---
    if stim_count_soc < stim_count_thresh || ismember('social', force_cut_list)
	% CASE: Too few trials -> Remove ENTIRE condition
        if stim_count_soc > 0
	    fprintf('Social count (%d) < threshold. Removing entire condition.\n', stim_count_soc);
	
	% Define time window: (First Event - Buffer) to (Last Event + Buffer)
	    t_start = EEG.event(soc_indices(1)).latency - (buffer_s * EEG.srate);
	    t_end = EEG.event(soc_indices(end)).latency + (buffer_s * EEG.srate);
	    regions_to_cut = [regions_to_cut; t_start, t_end];
        else
            fprintf('Social count is 0. Nothing to cut.\n');
        end

    elseif stim_count_soc >= stim_count_thresh && stim_count_soc < max_stim_count
	fprintf('Social count (%d) in warning zone. Searching for incomplete block...\n', stim_count_soc);
	
	lats = [EEG.event(soc_indices).latency];
	itis = diff(lats);
	
	% Find all indices where a block break occurs
	break_indices = find(itis > (block_break_s * EEG.srate));
	
	% Define block boundaries: [Start_Index, End_Index] (indices into 'lats')
	block_bounds = [];
	current_start = 1;
	
	for k = 1:length(break_indices)
	    current_end = break_indices(k);
	    block_bounds = [block_bounds; current_start, current_end];
	    current_start = current_end + 1;
	end
	block_bounds = [block_bounds; current_start, length(lats)]; % Add the final block
	
	% Check each block size
	found_bad_block = false;
	for b = 1:size(block_bounds, 1)
	    n_trials_in_block = block_bounds(b,2) - block_bounds(b,1) + 1;
	    
	    % If this block is significantly smaller than expected full block
	    if n_trials_in_block < trials_per_block
		fprintf(' -> Found incomplete Social block (Block %d: %d trials). Marking for removal.\n', b, n_trials_in_block);
		
		t_start = lats(block_bounds(b,1)) - (buffer_s * EEG.srate);
		t_end   = lats(block_bounds(b,2)) + (buffer_s * EEG.srate);
		regions_to_cut = [regions_to_cut; t_start, t_end];
		found_bad_block = true;
	    end
	end

        if ~found_bad_block
            fprintf('Social condition is incomplete (%d trials), but all detected blocks are full. Keeping condition as is.\n', stim_count_soc);
        end	
    end

    % --- STEP 3: Evaluate NONSOCIAL Condition ---
    if stim_count_nonsoc < stim_count_thresh || ismember('nonsocial', force_cut_list)
	% CASE: Too few trials -> Remove ENTIRE condition
        if stim_count_nonsoc > 0
	    fprintf('Nonsocial count (%d) < threshold. Removing entire condition.\n', stim_count_nonsoc);
	
	    t_start = EEG.event(nonsoc_indices(1)).latency - (buffer_s * EEG.srate);
	    t_end   = EEG.event(nonsoc_indices(end)).latency + (buffer_s * EEG.srate);
	    regions_to_cut = [regions_to_cut; t_start, t_end];
        else
            fprintf('Nonsocial count is 0. Nothing to cut.\n');
        end

    elseif stim_count_nonsoc >= stim_count_thresh && stim_count_nonsoc < max_stim_count
        fprintf('Nonsocial count (%d) in warning zone. Searching for incomplete block...\n', stim_count_nonsoc);
    
	lats = [EEG.event(nonsoc_indices).latency];
	itis = diff(lats);
	break_indices = find(itis > (block_break_s * EEG.srate));
	
	block_bounds = [];
	current_start = 1;
	for k = 1:length(break_indices)
	    current_end = break_indices(k);
	    block_bounds = [block_bounds; current_start, current_end];
	    current_start = current_end + 1;
	end
	block_bounds = [block_bounds; current_start, length(lats)]; 
	
	found_bad_block = false;
	for b = 1:size(block_bounds, 1)
	    n_trials_in_block = block_bounds(b,2) - block_bounds(b,1) + 1;
	    
	    if n_trials_in_block < trials_per_block
		fprintf(' -> Found incomplete Nonsocial block (Block %d: %d trials). Marking for removal.\n', b, n_trials_in_block);
		t_start = lats(block_bounds(b,1)) - (buffer_s * EEG.srate);
		t_end   = lats(block_bounds(b,2)) + (buffer_s * EEG.srate);
		regions_to_cut = [regions_to_cut; t_start, t_end];
		found_bad_block = true;
	    end
	end

        if ~found_bad_block
	    fprintf('Nonsocial condition is incomplete (%d trials), but all detected blocks are full. Keeping condition as is.\n', stim_count_nonsoc);
	end
    end

    % --- STEP 4: Execute Cuts ---
    if ~isempty(regions_to_cut)
	regions_to_cut(:, 1) = max(1, regions_to_cut(:, 1));
	regions_to_cut(:, 2) = min(EEG.pnts, regions_to_cut(:, 2));
	
	fprintf('Cutting %d regions from data...\n', size(regions_to_cut, 1));
	EEG = eeg_eegrej(EEG, regions_to_cut);
	EEG = eeg_checkset(EEG);
    else
	disp('Counts look good. No cuts needed.');
    end
end
