% this function is run inside MADE, so all EEGLab dependencies should be already imported and initialized

function [bad_channels] = preprocess_eeg_piece(eeg_piece, channel_locations, stimulus_timeoffset)
    down_sample = 1;
    sampling_rate = 1000;
    adjust_time_offset = 1;
    delete_outerlayer = 0;
    outerlayer_channel = {'16','15','12','13','8','31','26','25','30','32','60','64','61','62','56','57','63','41','46','45','48'}; % list of channels
    highpass = .1; % High-pass frequency
    lowpass  = 49; % Low-pass frequency. We recommend low-pass filter at/below line noise frequency (see manuscript for detail)
    stimulus_markers = {'S  1', 'S  2', 'S  3', 'S  4', 'S 41', 'S 42', 'S 43', ...
    'S 44', 'S 51', 'S 52', 'S 53', 'S 54'}; % enter the stimulus markers that need to be adjusted for time offset % fine only if we dont adjust for onset, not even used further in the code
    response_markers = {}; % enter the response makers that need to be adjusted for time offset % same as line above !!!

    %% Initialize EEG structure, output variables, and report table
    EEG=[]; %initialize eeg structure
    report_table = []; %report table that will be created and written to disk (appended) after processing completes for this participant
    reference_used_for_faster=[]; % reference channel used for running faster to identify bad channel/s
    faster_bad_channels=[]; % number of bad channel/s identified by faster
    ica_preparation_bad_channels=[]; % number of bad channel/s due to channel/s exceeding xx% of artifacted epochs
    length_ica_data=[]; % length of data (in second) fed into ICA decomposition
    total_ICs=[]; % total independent components (ICs)
    ICs_removed=[]; % number of artifacted ICs
    total_epochs_before_artifact_rejection=[];
    total_epochs_after_artifact_rejection=[];
    total_channels_interpolated=[]; % total_channels_interpolated=faster_bad_channels+ica_preparation_bad_channels

    %fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names{subject});

    %% STEP 1: Import EEG data file and relevant information

    %load in raw data
    %EEG = pop_loadbv(rawdata_location, eeg_piece);
    EEG = eeg_checkset(eeg_piece);

    %% STEP 4: Change sampling rate
    if down_sample==1
	if floor(sampling_rate) > EEG.srate
	    error ('Sampling rate cannot be higher than recorded sampling rate');
	elseif floor(sampling_rate) ~= EEG.srate
	    EEG = pop_resample( EEG, sampling_rate);
	    EEG = eeg_checkset( EEG );
	end
    end

    %make a copy of GSR and sync channels, then delete from eeg structure

    gsrChan = EEG.data(64, :);
    EEG = pop_select( EEG,'nochannel', 64);
    EEG = eeg_checkset( EEG );

    syncChan = EEG.data(64, :);
    EEG = pop_select( EEG,'nochannel', 64);
    EEG = eeg_checkset( EEG );

    %add in ref channels
    origData = EEG.data;
    [origData_NumRows, origData_NumCols] = size(origData);
    EEG.data = NaN(origData_NumRows+1, origData_NumCols);
    EEG.data(1,:) = 0; %add ref as zeros
    disp('processs eeg piece DEBUG 1');
    EEG.data(2:end,:) = origData; %copy over orig EEG data
    %delete ground from newChanLocs
    modNewChanlocs = channel_locations(2:end);

    %replace chanlocs with
    EEG.chanlocs = modNewChanlocs;
    EEG.nbchan = EEG.nbchan+1;
    %%%%
    EEG = eeg_checkset( EEG );

    EEG = eeg_checkset( EEG );

    % Check whether the channel locations were properly imported. The EEG signals and channel numbers should be same.
    if size(EEG.data, 1) ~= length(EEG.chanlocs)
	error('The size of the data does not match with channel numbers.');
    end

    %% STEP 1b: convert all type field markers to string (if not already)

    %loop through all the type markes, if numeric, convert to string
    % (Given that this script assumes that "type" field markers are strings, we need to
    % convert all type field markers to string, in case they are not
    % already)
    for atm=1:length({EEG.event.type})
	if isnumeric(EEG.event(atm).type)
	    EEG.event(atm).type = num2str(EEG.event(atm).type);
	end
    end

    %% STEP 3: Adjust anti-aliasing and task related time offset
    if adjust_time_offset==1
    %    %%%%%% adjust anti-aliasing filter time offset
    %    if filter_timeoffset~=0
    %        for aafto=1:length(EEG.event)
    %            EEG.event(aafto).latency=EEG.event(aafto).latency+(filter_timeoffset/1000)*EEG.srate;
    %        end
    %    end
	% adjust stimulus time offset
	if stimulus_timeoffset~=0
	    for sto=1:length(EEG.event)
		for sm=1:length(stimulus_markers)
		    if strcmp(EEG.event(sto).type, stimulus_markers{sm})
			EEG.event(sto).latency=EEG.event(sto).latency+(stimulus_timeoffset/1000)*EEG.srate;
		    end
		end
	    end
	end
    %    % adjust response time offset
    %    if response_timeoffset~=0
    %        for rto=1:length(EEG.event)
    %            for rm=1:length(response_markers)
    %                if strcmp(EEG.event(rto).type, response_markers{rm})
    %                    EEG.event(rto).latency=EEG.event(rto).latency-(response_timeoffset/1000)*EEG.srate;
    %                end
    %            end
    %        end
    %    end
    end

    %% STEP 5: Delete outer layer of channels
    chans_labels=cell(1,EEG.nbchan);
    for i=1:EEG.nbchan
	chans_labels{i}= EEG.chanlocs(i).labels;
    end

    if delete_outerlayer==1
	[chans,chansidx] = ismember(outerlayer_channel, chans_labels);
	outerlayer_channel_idx = chansidx(chansidx ~= 0);
	if isempty(outerlayer_channel_idx)==1
	    error(['None of the outer layer channels present in channel locations of data.'...
		' Make sure outer layer channels are present in channel labels of data (EEG.chanlocs.labels).']);
	else
	    EEG = pop_select( EEG,'nochannel', outerlayer_channel_idx);
	    EEG = eeg_checkset( EEG );
	end
    end

    %% STEP 6: Filter data
    % Calculate filter order using the formula: m = dF / (df / fs), where m = filter order,
    % df = transition band width, dF = normalized transition width, fs = sampling rate
    % dF is specific for the window type. Hamming window dF = 3.3

    high_transband = highpass; % high pass transition band
    low_transband = 10; % low pass transition band

    hp_fl_order = 3.3 / (high_transband / EEG.srate);
    lp_fl_order = 3.3 / (low_transband / EEG.srate);

    % Round filter order to next higher even integer. Filter order is always even integer.
    if mod(floor(hp_fl_order),2) == 0
	hp_fl_order=floor(hp_fl_order);
    elseif mod(floor(hp_fl_order),2) == 1
	hp_fl_order=floor(hp_fl_order)+1;
    end

    if mod(floor(lp_fl_order),2) == 0
	lp_fl_order=floor(lp_fl_order)+2;
    elseif mod(floor(lp_fl_order),2) == 1
	lp_fl_order=floor(lp_fl_order)+1;
    end

    % Calculate cutoff frequency
    high_cutoff = highpass/2;
    low_cutoff = lowpass + (low_transband/2);

    % Performing high pass filtering
    EEG = eeg_checkset( EEG );
    EEG = pop_firws(EEG, 'fcutoff', high_cutoff, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', hp_fl_order, 'minphase', 0);
    EEG = eeg_checkset( EEG );

    % Performing low pass filtering
    EEG = eeg_checkset( EEG );
    EEG = pop_firws(EEG, 'fcutoff', low_cutoff, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', lp_fl_order, 'minphase', 0);
    EEG = eeg_checkset( EEG );

    %% STEP 7: Run faster to find bad channels
    % First check whether reference channel (i.e. zeroed channels) is present in data
    % reference channel is needed to run faster
    ref_chan=[]; FASTbadChans=[]; all_chan_bad_FAST=0;
    ref_chan=find(any(EEG.data, 2)==0);
    if numel(ref_chan)>1
	error(['There are more than 1 zeroed channel (i.e. zero value throughout recording) in data.'...
	    ' Only reference channel should be zeroed channel. Delete the zeroed channel/s which is not reference channel.']);
    elseif numel(ref_chan)==1
	list_properties = channel_properties(EEG, 1:EEG.nbchan, ref_chan); % run faster
	FASTbadIdx=min_z(list_properties);
	FASTbadChans=find(FASTbadIdx==1);
	FASTbadChans=FASTbadChans(FASTbadChans~=ref_chan);
	reference_used_for_faster={EEG.chanlocs(ref_chan).labels};
	% EEG = pop_select( EEG,'nochannel', ref_chan); % a bug [kia
	% removed it as George said]
	EEG = eeg_checkset(EEG);
	channels_analysed=EEG.chanlocs; % keep full channel locations to use later for interpolation of bad channels
    elseif numel(ref_chan)==0
	warning('Reference channel is not present in data. channel 1 will be used as reference channel.');
	ref_chan=find(strcmp({EEG.chanlocs.labels}, '1')); % find Cz channel index
	EEG_copy=[];
	EEG_copy=EEG; % make a copy of the dataset
	EEG_copy = pop_reref( EEG_copy, ref_chan,'keepref','on'); % rerefer to Cz in copied dataset
	EEG_copy = eeg_checkset(EEG_copy);
	list_properties = channel_properties(EEG_copy, 1:EEG_copy.nbchan, ref_chan); % run faster on copied dataset
	FASTbadIdx=min_z(list_properties);
	FASTbadChans=find(FASTbadIdx==1);
	channels_analysed=EEG.chanlocs;
	reference_used_for_faster={EEG.chanlocs(ref_chan).labels};
    end
    bad_channels = FASTbadChans
end
