function [stim_count_nonsoc, stim_count_soc] = stim_cnt_check_thrive(path)

    stim_events_nonsoc = {'S 41', 'S 42', 'S 43', 'S 44'};
    stim_events_soc = {'S 51', 'S 52', 'S 53', 'S 54'};
    stim_count_nonsoc = 0;
    stim_count_soc = 0;

    % fprintf('sub-%s: Loading files for deviated EEG...\n', sub);
    vhdr_files = dir(fullfile(path, '*.vhdr'));
    for f = 1:length(vhdr_files)
	fname = vhdr_files(f).name;
	fpath = vhdr_files(f).folder;
	try
	    % Load EEG data
	    EEG = pop_loadbv(fpath, fname);
	    EEG = eeg_checkset(EEG);
	    % Select stimulus events
	    stim_count_nonsoc = stim_count_nonsoc + sum(ismember(string({EEG.event.type}), stim_events_nonsoc));
	    stim_count_soc = stim_count_soc + sum(ismember(string({EEG.event.type}), stim_events_soc));	
	catch ME
	    fprintf('FAILED to load or process files for stim count check!. Error: %s\n', ME.message);
	    stim_count_nonsoc = NaN; % Use NaN to indicate error
	    stim_count_soc = NaN; % Use NaN to indicate error
	end
    end
    fprintf('Total %d NONSOC and %d SOC stimulus events found.\n', stim_count_nonsoc, stim_count_soc);
end
