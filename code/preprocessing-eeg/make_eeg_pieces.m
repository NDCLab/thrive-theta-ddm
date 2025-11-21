function [eeg_pieces] = make_eeg_pieces(EEG)

    % Step 1: Extract boundary event latencies
    [~, boundary_indices] = pop_selectevent(EEG, 'type', 'boundary');
    boundary_events = EEG.event(boundary_indices);
    latencies = [boundary_events.latency] / 1000;
    disp(latencies);

    % Step 2: Identify merging points
    % Initialize a vector to store unique merging points
    unique_latencies = [];

    % Iterate through the sorted latencies
    for i = 1:length(latencies)
	if i == 1
	    % Always keep the first point
	    unique_latencies = [unique_latencies, latencies(i)];
	else
	    % Check if the current point is too close to the previous point
	    if abs(latencies(i) - latencies(i-1)) >= 1*1.0e-3
		unique_latencies = [unique_latencies, latencies(i)];
	    end
	end
    end

    merging_points = unique_latencies(unique_latencies > 1*1.0e-3);
    disp(merging_points);

    % Step 3: Split the EEG data into pieces
    eeg_pieces = {};

    % The first piece is from the start to the first merging point
    start_latency = 0;
    end_latency = merging_points(1);
    eeg_pieces{1} = pop_select(EEG, 'time', [start_latency, end_latency]);
%   eeg_pieces{1} = eeg_checkset(eeg_pieces{1});
    % Loop through the remaining merging points
    for i = 2:length(merging_points)
	start_latency = merging_points(i-1);
	end_latency = merging_points(i);
	eeg_pieces{i} = pop_select(EEG, 'time', [start_latency, end_latency]);
 %       eeg_pieces{i} = eeg_checkset(eeg_pieces{i});
    end

    % The last piece is from the last merging point to the end of the EEG
    last_start_latency = merging_points(end);
    eeg_pieces{end+1} = pop_select(EEG, 'time', [last_start_latency, EEG.times(end) / 1000]);
 %   eeg_pieces{end} = eeg_checkset(eeg_pieces{end});
    % Display the number of EEG pieces
    disp(['Number of EEG pieces: ', num2str(length(eeg_pieces))]);

    cum_length = 0;
    for i=1:length(eeg_pieces)
	cum_length = cum_length + length(eeg_pieces{i}.times) / 1000;
    end
    assert((length(EEG.times) / 1000 - cum_length) <= 1, ...
	   'WARNING: The difference between the EEG duration and cumulative piece length exceeds 1 second!');

end
