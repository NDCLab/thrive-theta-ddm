%%to Run on FIU HPC%
% create a local cluster object
%cluster = parcluster('local');

% start matlabpool with max workers set in the slurm file
%parpool(cluster, str2num(getenv('SLURM_CPUS_ON_NODE')))
%Compute TF, ITPS, ICPS, and wPLI measures for EEG data
%Maureen Bowers 6/29/2021 based on scripts by Ranjan Debnath

%clear all;
%clc;

%%
%%%%% Setting paths %%%%%

main_dir = '/home/data/NDClab/analyses/thrive-theta-ddm'; %directory on the HPC
%session = 's2_r1';
tf_output_dir = 'TF';
itps_output_dir = 'ITPS';
icps_output_dir = 'ICPS';
wpli_output_dir = 'wPLI';
%1. Data Location
data_location = [main_dir filesep 'derivatives' filesep 'preprocessed' filesep 'csd_data' filesep session filesep];

%2. Save Data Location
save_location = [main_dir filesep 'derivatives' filesep 'preprocessed' filesep 'TF_outputs' filesep session filesep 'resp' filesep 'seed_1' filesep];
%disp(save_location)
% Create output folders to save data
if exist(save_location, 'dir') == 0
    mkdir(save_location);
end

% Create output folders to save data
if exist([save_location filesep tf_output_dir], 'dir') == 0
    mkdir([save_location filesep tf_output_dir]);
end

% Create output folders to save data
if exist([save_location filesep itps_output_dir], 'dir') == 0
    mkdir([save_location filesep itps_output_dir]);
end

% Create output folders to save data
if exist([save_location filesep icps_output_dir], 'dir') == 0
    mkdir([save_location filesep icps_output_dir]);
end

% Create output folders to save data
if exist([save_location filesep wpli_output_dir], 'dir') == 0
    mkdir([save_location filesep wpli_output_dir]);
end

%3. Scripts Location
scripts_location = [main_dir filesep 'code' filesep 'preprocessing-eeg' filesep 'tf'];

% addpath(genpath([main_dir filesep 'MADE-EEG-preprocessing-pipeline']));% enter the path of the EEGLAB folder in this line
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeg_preprocessing'));% enter the path of the folder in this line

%Location of "EEG
% addpath(genpath([main_dir filesep 'eeglab13_4_4b']));% enter the path of the EEGLAB folder in this line
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b'));% enter the path of the EEGLAB folder in this line
addpath(genpath('/home/data/NDClab/analyses/thrive-theta-ddm/code/matlab'));
%remove path to octave functions inside matlab to prevent errors when
% rmpath([main_dir filesep 'eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])
rmpath(['/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])

%%
%%%%% Information about dataset and analysis procedures %%%%%%
%5. Number of channels
nbchan = 64;

%6. What kind of data? Resting State or Event-Related Data
RestorEvent = 0; %1 = rest, 0 = event

%7. What are your conditions of interest if using Event-Related Data? This
%naming convention should come from the Edit_events.m script provided.
%Note: not needed for resting state data.
%Super_Conds = {
%{'resp_s_i_0'},...
%{'resp_s_i_1'},...
%{'resp_s_c_1'},...
%{'resp_ns_i_0'},...
%{'resp_ns_i_1'},...
%{'resp_ns_c_1'},...
%}
% RESPONSE CONDITIONS
Conds = {
%'resp_s_i_0',...
%'resp_s_i_1',...
%'resp_s_c_1',...
%'resp_ns_i_0',...
%'resp_ns_i_1',...
%'resp_ns_c_1',...
condition_to_compute
}; % resp_i_0 = incong error response; resp_i_1 = incong correct response

% STIMULUS CONDITIONS
%Conds = {
%'stim_s_i_0',...
%'stim_s_i_1',...
%'stim_s_c_1',...
%'stim_ns_i_0',...
%'stim_ns_i_1',...
%'stim_ns_c_1',...
%}; % resp_i_0 = incong error stim; resp_i_1 = incong correct response

%8. Minimum number of trials to analyze
mintrialnum = 6; %If the participant does not have enough trials in a condition based on this cutoff, a "notenoughdata.mat" file will be saved into save_location.

%9. Would you like to baseline correct your data? NOTE: ICPS calculated over time will not be baseline corrected.
BaselineCorrect = 1; %1=Yes, 0=No

%10. What time period would you like to use to baseline correct
BaselineTime = [-400 -200];
%Put in time in ms for event-related data. For example, if you want to baseline correct from -100 to 0ms before the 
% event for evernt-related paradigms, put [-100 0].

%11. Would you like to downsample the output to 125Hz? This is done after the time-frequency computations and will have minor impacts on the resolution. 
%We recommend downsampling to reduce file size for ease of storage. All data will be initially downsampled to 250Hz. 
Downsample = 1; %0=No, 1=Yes  

%12. Dataset Name - will be appended to the saved files
DatasetName = '';

%%
%%%%% Settings for Time-Frequency Measures %%%%%

%13. Minimum and Maximum Frequency, Number of Frequency Bins, and range cycles to calculate complex Morlet wavelet decomposition
min_freq = 1;
max_freq = 30; % check
num_frex = 59; % number of frequency bins between minimum and maximum frequency
range_cycles = [3 10]; % wavelet cycles: min 3 max 10

%%
%%%%% Questions about phase-based measures %%%%%

%14. Would you like to calculate inter-trial phase synchrony (ITPS) in addition to TF?
ITPS_calc = 1; %(1=yes,0=no)

%15. Would you like to subsample trials? This is recommended for event-related paradigms, 
% especially when there are uneven numbers of trials in conditions. 
Subsample = 1; %1=Yes, 0=No

%How many trials to pull for each subsample? 
NumTrialsPulled = 6;

%How many times to do subsampling? We recommend doing at least 10 subsamples and to have the possiblility of using all your data
% (e.g., if you have 150 trials, do 15 subsamples of 10 trials).
NumSubsamples = 100; %

%16. Would you like to calculate inter-channel phase synchrony (ICPS or wPLI) in addition to TF?
ICPS_calc = 1; %(1=yes,0=no) 

%17. Select which connectivity measures to run (set both to 1 to run both)
Compute_ICPS = 1; % Coherence
Compute_wPLI = 1; % Weighted Phase Lag Index

%18. Inter-channel phase synchrony over trials or connectivity over time?
% Kia: over trials results in a frequency by time matrix and over time
% results in a frequency by trials, which removes the time dimension. 
TimeOrTrials = 1; %0 = over time, 1 = over trials
% NOTE: over time calculations will not be able to be subsampled or
% baseline corrected

%19. Type of Connectivity to compute 
ConnectType = 1; %(0= all-to-all connectivity; 1=seed-based connectivity)
%If seed based, choose seed and which electrodes to compute connectivity:
Seed = '1';
Elecs4Connect = { '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '34' '35' '36' '37' '38' '39' '40' '41' '42' '43' '44' '45' '46' '47' '48' '49' '50' '51' '52' '53' '54' '55' '56' '57' '58' '59' '60' '61' '62' '63' '64'};

% 20. Create List of subjects to loop through
subnum = dir([data_location filesep '*.set']); % Use regex to find your files 
subject = {subnum.name};
for ii=1:length(subject)
    subject_list{ii}=subject{ii};   
end
% already_processed = dir([save_location filesep 'TF' filesep '*' Conds{1} '*.mat']);
already_processed = dir([save_location filesep wpli_output_dir filesep '*' Conds{1} '*.mat']);
already_processed = {already_processed.name};

% Initialize an empty cell array to store the extracted parts
processed_ids = {};

% Loop through each filename and extract the sub-xxxxxxx part
for i = 1:length(already_processed)
    filename = already_processed{i};
    % Use regular expression to find the sub-xxxxxxx part
    match = regexp(filename, '^sub-\d+', 'match');
    if ~isempty(match)
        processed_ids{end+1} = match{1}; % Append the match to the cell array
    end
end

%TrialNums = struct('subject', {}, 'condition', {}, 'TrialNum', {});
%TrialNums(34*length(Conds)).subject = '';  % Pre-allocate for max possible size
%TrialNums(34*length(Conds)).condition = '';
%TrialNums(34*length(Conds)).TrialNum = 0;
% We recommend checking that your subject list looks like you would expect:
% subject_list

%eeglab % Loading EEGLAB
TrialNums = zeros(length(subject_list), length(Conds));
%rng(2, 'twister'); % to fix seed

%%%%%%%%%%%%%%%%%%%% COMPUTATIONS BEGIN BELOW HERE %%%%%%%%%%%%%%%%
%% loop through all subject
parfor sub=1:length(subject_list)
%parfor sub=1:2
   % rng(2, 'twister'); % if you want to fix the seed put it inside the loop as well    
    local_TrialNums = zeros(1, length(Conds));
       
    % Initialize objects for this participant:
    timefreqs_data = [];
    phase_data=[];
    ITPS_all=[];
    ICPS_all=[];
    wPLI_all=[];
    EEG=[];
    
    subject = subject_list{sub};
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', sub, subject);

    % Check if subject was already processed
    match = regexp(subject, '^sub-\d+', 'match');
    
    % Ensure match is a string
    if ~isempty(match)
        match = match{1}; % Extract the string from the cell array
    else
        match = ''; % Handle cases where the pattern is not found
    end
    
    if ismember(match, processed_ids)
        fprintf('%s: is already processed. Skipping...\n', match);
        continue;
    end

    % Load data
    EEG=pop_loadset('filename', [subject], 'filepath', data_location);
    EEG = pop_selectevent(EEG, 'latency', '-.1 <= .1', 'deleteevents', 'on');
    
    EEG = pop_editeventfield(EEG, 'indices',  strcat('1:', int2str(length(EEG.event))), 'Condition', 'NaN');
    EEG = eeg_checkset(EEG);
    
    for t=1:length(EEG.event)

    % Condition that we want to concatenate into one variable - THIS WILL
    % NEED TO BE CHANGED FOR EACH STUDY !!!!
        eventType = {EEG.event.eventType};
        observation = {EEG.event.observation};
        congruency = {EEG.event.congruency};
        accuracy = {EEG.event.accuracy};
    
        % Selecting specific variables for this trial
        vars_to_join = { % "resp_s_i_1" or "stim_ns_c_0"
            eventType{t}, ...
            observation{t}, ...
            congruency{t}, ...
            num2str(accuracy{t}) ...
            };
        EEG.event(t).Condition = strjoin(vars_to_join, '_');
        % resulting pattern will look like stim_s_i_0 or resp_ns_c_1
    end
    
    %Downsample to 250Hz
    EEG = pop_resample(EEG, 250);
    EEG = eeg_checkset(EEG);
    
    % Keep only markers of interest
    try
        EEG = pop_selectevent(EEG, 'Condition', Conds, 'deleteevents', 'on', 'deleteepochs', 'on', 'invertepochs', 'off');
        EEG = eeg_checkset(EEG);

    catch ME
        warning('Error occurred when selecting only markers of interest for subject %d: %s', subject, ME.message);
        continue; % This will skip to the next iteration of the loop (next subject)
    end
    
    % Trials that are labeled as responded and validRt (trials that have larger than 150 ms RT) 
    % will only be included.
    try 
        EEG = pop_selectevent(EEG, 'validRt', 1, 'extraResponse', 0, 'deleteevents', 'on', 'deleteepochs', 'on', 'invertepochs', 'off');
        EEG = eeg_checkset(EEG);

    catch ME
        warning('Error occurred when selecting only valid trials for subject %d: %s', subject, ME.message);
        continue; % This will skip to the next iteration of the loop (next subject)
    end

    try 
        EEG = pop_selectevent(EEG, 'responded', 1, 'deleteevents', 'on', 'deleteepochs', 'on', 'invertepochs', 'off');
        EEG = eeg_checkset(EEG);

    catch ME
        warning('Error occurred when selecting only valid trials (responded) for subject %d: %s', subject, ME.message);
        continue; % This will skip to the next iteration of the loop (next subject)
    end
    % If no error occurred during selecting only valid trials, the code
    % will continue from here.
    
    %save times
    time = EEG.times;
    
    %save channel locations
    channel_location = EEG.chanlocs;
    
    %% Compute complex morlet wavelet time frequency decomposition
    cd(scripts_location)
    % timefreq_script;
    % INPUT
    data = EEG.data; % data = data matrix; channel x time/data points x trials
    srate = EEG.srate; % srate = sampling rate;
    % times = time vector / data time window;
    % baseline = [xx xx]; time range/window to be used for baseline normalization
    % min_freq = lowest frequency to analysis;
    % max_freq = highest frequency to analysis;
    
    % OUTPUT
    % timefreq_data = frequency x time x channel matrix;
    % times = time vector;
    % frex = frequency vector;
    
    % frequencies vector
    frex = linspace(min_freq, max_freq, num_frex);
    frequency = frex;
    
    % cycles vector
    cylvec = linspace(range_cycles(1),range_cycles(end), num_frex)./ (2*pi*frex);
    
    %% wavelet parameters
    wavtime = -1:1/srate:1; % length of wavelet
    half_wave = (length(wavtime)-1)/2;
    
    %% FFT parameters
    nWave = length(wavtime);
    nData = size(data, 2);
    nConv = nWave + nData - 1;
    %% Prepare data for wavelet
    %Resting Data - only one trial type
    %All to all connectivity
    if ConnectType == 0
        %Loop through epochs and pull only epochs with appropriate number of
        %channels
        epochs2remove = [];
        for e=1:size(data, 3) %begin loop through epochs
            %Only keep epoch if none of channels are NaNs
            for elec = 1:nbchan %begin loop through electrodes
                %check if the electrode is NaN within that epoch
                if sum(isnan(EEG.data(elec, :, e)))
                    %add that epoch to vector
                    epochs2remove = [epochs2remove e];
                else
                    %do nothing
                end %end if loop for if channel within epoch is NaN
            end %end loop through electrodes
        end %end loop through epochs
        %remove all epochs with NaNs for channels of interest
        data(:,:,unique(epochs2remove))=[];
        %Seed-based connectivity
    elseif ConnectType == 1
        epochs2remove=[];
        for e=1:size(data,3) %begin loop through epochs
            %Channel indices
            %concatenate Elecs4Connect and Seed to have vector of all
            %channels of interest
            Elecs4SeedBased = [Seed Elecs4Connect];
            Chan_idx = zeros(1, length(Elecs4SeedBased));
            %loop through and find indices of these electrodes
            for i=1:length(Elecs4SeedBased)
                %fprintf('FINE HERE');
                Chan_idx(i)= find(strcmp({EEG.chanlocs.labels}, Elecs4SeedBased{i}));
                %fprintf('FINE HERE STILL!');
            end
            %fprintf('FINE HERE STILL!');
            %do any of these channels have NaNs in them?
            if min(sum(isnan(EEG.data(Chan_idx,:,e)))) > 0
                %if so, add that epoch to vector of epochs to be removed
                epochs2remove = [epochs2remove e];
            else
                %do nothing
            end %end if loop for if channel within epoch is NaN
            %fprintf('FINE HERE STILL!');
        end %end loop through epochs
        %remove all epochs with NaNs for channels of interest
        data(:,:,unique(epochs2remove))=[];
        %keep only channels of interest
        data = data(Chan_idx,:,:);
        %fprintf('FINE HERE STILL!');
    end % end if all-to-all connectivity for rest
    
    fprintf('\n\n\n*** Calculating TF for subject %d (%s) ***\n\n\n', sub, subject);
    
    %% Run wavelet convolution
    if RestorEvent ==1
        
        %TrialNums.subject = subject;
        %TrialNums.TrialNum = size(data,3);
        
        if size(data,3)>= mintrialnum
            timefreqs = zeros(length(frex), nData, size(data, 3));
            timefreq_strial = zeros(length(frex), nData, size(data, 3), size(data,1));
            for ch=1:size(data, 1) % Loop through all channels
                for fi=1:length(frex) % loop through all frequencies
                    trial_conv = zeros(nData, size(data, 3));
                    
                    %% Create wavelet
                    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*cylvec(fi)^2));
                    waveletX = fft(wavelet, nConv); % fft of wavelet
                    waveletX = waveletX ./ max(waveletX); % normalize fft of wavelet
                    
                    %% Loop through all trials
                    for trl=1:size(data, 3)
                        temp_data = fft(squeeze(data(ch,:,trl)), nConv);
                        
                        %% run convolution
                        temp_conv = ifft(waveletX .* temp_data);
                        temp_conv = temp_conv(half_wave+1:end-half_wave);
                        
                        %% data for phase analysis
                        trial_conv (:,trl) = temp_conv;
                        
                        %% compute power
                        trial_temppow (:,trl) = abs(temp_conv(half_wave+1:end-half_wave)).^2;
                    end
                    %Build dataset for phase data
                    phase_data(fi,:,:) = trial_conv;
                            
                    %Build dataset for TF
                    timefreqs(fi,:,:) = trial_temppow;
    
                    
                end
                %Build dataset for phase data with channels
                phase_data_strial_ch(:,:,:,ch) = phase_data;
                
                %Build dataset for TF with channels
                timefreqs_strial_ch(:,:,:,ch) = timefreqs;
            end %end loop through channels
            
        %Trial average timefrequency data
        timefreqs_data = squeeze(mean(timefreqs_strial_ch,3));
        
        else
            data2save='';
            save_data =[save_location, subject(1:end-4), '_notenoughdata.mat'];
            %save([save_location subject(1:end-4) '_notenoughdata.mat'],data2save)
            parsave(save_data, 'data2save', data2save);

        end %end if there's enough trials
        
        %Baseline Correction
        if BaselineCorrect == 1
            
            %% baseline time indices
            basetimeidx   = dsearchn(EEG.times', BaselineTime');
            
            Baseline = squeeze(mean(timefreqs_data(:,basetimeidx(1):basetimeidx(end),:),2));
            %loop through samples
            for t=1:size(timefreqs_data,2)
                timefreqs_baselinecorr(:,t,:) = 10*log10( squeeze(timefreqs_data(:,t,:)) ./ Baseline);
            end %end loop through frequencies
            
            if Downsample ==0
                %save out trial averaged and baseline corrected TF data for this subject for this condition
                save_data =[save_location, tf_output_dir, filesep, subject(1:end-4),DatasetName,'_TF_baselinecorrected'];
                %save(save_data, 'timefreqs_baselinecorr', 'frequency', 'time','channel_location', '-v7.3');
                 parsave(save_data, ...
                 'timefreqs_baselinecorr', timefreqs_baselinecorr, ...
                 'frequency', frequency, ...
                 'time', time, ...
                 'channel_location', channel_location);
            elseif Downsample ==1
                
                %Downsample to 125Hz
                timefreqs_baselinecorr = timefreqs_baselinecorr(:,1:2:size(timefreqs_baselinecorr,2),:);
                
                %Downsample time variable
                ds_time = downsample(time,2);
                
                %save out trial averaged, downsampled, baseline corrected TF data for this subject for this condition
                save_data =[save_location, tf_output_dir, filesep, subject(1:end-4),DatasetName,'_TF_baselinecorrected'];
                %save(save_data, 'timefreqs_baselinecorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                 parsave(save_data, ...
                 'timefreqs_baselinecorr', timefreqs_baselinecorr, ...
                 'frequency', frequency, ...
                 'ds_time', ds_time, ...
                 'channel_location', channel_location);
            end
            
        elseif BaselineCorrect == 0
            
            if Downsample ==0
                %save out trial averaged and non-baseline corrected TF data
                save_data =[save_location, tf_output_dir, filesep, subject(1:end-4),DatasetName,'_TF_nobaselinecorrection'];
                %save(save_data, 'timefreqs_data', 'frequency', 'time','channel_location', '-v7.3');
                 parsave(save_data, ...
                 'timefreqs_data', timefreqs_data, ...
                 'frequency', frequency, ...
                 'time', time, ...
                 'channel_location', channel_location);
            elseif Downsample ==1
                %Downsample to 125Hz
                timefreqs_data = timefreqs_data(:,1:2:size(timefreqs_data,2),:);
                
                %Downsample time variable
                ds_time = downsample(time,2);
                
                %save out trial averaged and non-baseline corrected TF data
                save_data =[save_location, tf_output_dir, filesep, subject(1:end-4),DatasetName,'_TF_nobaselinecorrection'];
                %save(save_data, 'timefreqs_data', 'frequency', 'ds_time','channel_location', '-v7.3');
                 parsave(save_data, ...
                 'timefreqs_data', timefreqs_data, ...
                 'frequency', frequency, ...
                 'ds_time', ds_time, ...
                 'channel_location', channel_location);
            end
        end %end if you want baseline correction
           
    elseif RestorEvent==0
        for cond=1:length(Conds)
            
            %pull out condition data
	try
            cond_data = []; cond_EEG=[];
            cond_EEG = pop_selectevent(EEG, 'Condition',(Conds{cond}), 'deleteevents','on','deleteepochs','on','invertepochs','off');
            cond_data = eeg_checkset(cond_EEG.data);%add pulling out condition specific data
	    if isempty(cond_EEG.data) || size(cond_EEG.data, 3) == 0
                fprintf('\nWARNING: No epochs found for subject %s, condition %s. Skipping...\n', subject, Conds{cond});
                continue; % Skip to next condition
            end
	catch ME
            fprintf('\nERROR processing subject %s, condition %s: %s\n', subject, Conds{cond}, ME.message);
	end
            fprintf('\n* Condition %s data shape: %d\n *\n', Conds{cond}, size(cond_data, 3));
            %TrialNums(sub, cond) = size(cond_data, 3);
            %Create TrialNums structure to be saved out
           % idx = sub+(cond-1)+((sub-1)*(length(Conds)-1));
            %TrialNums(idx).subject = subject;
            %TrialNums(idx).condition = Conds{cond};
            %TrialNums(idx).TrialNum = size(cond_data,3);
   %         local_TrialNums(cond).subject = subject;
    %        local_TrialNums(cond).condition = Conds{cond};
     %       local_TrialNums(cond).TrialNum = size(cond_data,3);
  %          TrialNums_cell{sub} = local_TrialNums;             
	    local_TrialNums(cond) = size(cond_data, 3);
            %initialize matrices
            timefreqs = zeros(length(frex), nData, size(cond_data, 3));
            phase_data = zeros(length(frex), nData, size(cond_data, 3));
            phase_data_strial_ch = zeros(length(frex),  nData, size(cond_data, 3), size(cond_data,1));
            timefreqs_strial_ch = zeros(length(frex),  nData, size(cond_data, 3), size(cond_data,1));
            timefreqs_data = zeros(length(frex),  nData,  size(cond_data,1));
%            fprintf('FINE HERE STILL'); 
            if size(cond_data,3)< mintrialnum
                data2save='';
		save_data = [save_location subject(1:end-4) '_notenoughdata.mat']
                %save([save_location subject(1:end-4) '_notenoughdata.mat'],data2save)
		parsave(save_data, 'data2save', data2save);
                break %go onto next subject if one condition doesn't have enough trials
            else
                    for ch=1:size(cond_data, 1) % Loop through all channels
                        for fi=1:length(frex) % loop through all frequencies
                            trial_conv = zeros(nData, size(cond_data, 3));
                            trial_temppow = zeros(nData, size(cond_data, 3));
                            
                            %% Create wavelet
                            wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*cylvec(fi)^2));
                            waveletX = fft(wavelet, nConv); % fft of wavelet
                            waveletX = waveletX ./ max(waveletX); % normalize fft of wavelet
                            
                            %% Loop through all trials
                            for trl=1:size(cond_data, 3)
                                temp_data = fft(squeeze(cond_data(ch,:,trl)), nConv);
                                
                                %% run convolution
                                temp_conv = ifft(waveletX .* temp_data);
                                temp_conv2 = temp_conv(half_wave+1:end-half_wave);
                                
                                %% data for phase analysis
                                trial_conv (:,trl) = temp_conv2;
                                
                                %% compute power
                                trial_temppow (:,trl) = abs(temp_conv(half_wave+1:end-half_wave)).^2;
                            end
                            
                            %Build dataset for phase data
                            phase_data(fi,:,:) = trial_conv;
                            
                            %Build dataset for TF
                            timefreqs(fi,:,:) = trial_temppow;
                            
                        end
                        
                        %Build dataset for phase data with channels
                        phase_data_strial_ch(:,:,:,ch) = phase_data;
                        
                        %Build dataset for TF with channels
                        timefreqs_strial_ch(:,:,:,ch) = timefreqs;
                        
                        
                    end %end loop through channels
                    
                    %Trial average timefrequency data
                    timefreqs_data = squeeze(mean(timefreqs_strial_ch,3));
                        
                        %Baseline Correction
                        if BaselineCorrect == 1
                            
                            %% baseline time indices
                            basetimeidx   = dsearchn(EEG.times', BaselineTime');
                                
                            Baseline = squeeze(mean(timefreqs_data(:,basetimeidx(1):basetimeidx(end),:),2));
                            %loop through samples    
                            for t=1:size(timefreqs_data,2)
                                timefreqs_baselinecorr(:,t,:) = 10*log10( squeeze(timefreqs_data(:,t,:)) ./ Baseline);
                            end %end loop through frequencies
                            if Downsample ==0
                                %save out trial averaged and baseline corrected TF data for this subject for this condition    
                                save_data =[save_location, tf_output_dir, filesep, subject(1:end-4),DatasetName,'_TF_baselinecorrected_', 'condition_',Conds{cond}];
                                %save(save_data, 'timefreqs_baselinecorr', 'frequency', 'time','channel_location', '-v7.3');
                                parsave(save_data, ...
                                'timefreqs_baselinecorr', timefreqs_baselinecorr, ...
                                'frequency', frequency, ...
                                'time', time, ...
                                'channel_location', channel_location);
                            elseif Downsample ==1
                                
                                %Downsample to 125Hz
                                timefreqs_baselinecorr = timefreqs_baselinecorr(:,1:2:size(timefreqs_baselinecorr,2),:);
                                
                                %Downsample time variable
                                ds_time = downsample(time,2);
                                
                                %save out trial averaged, downsampled, baseline corrected TF data for this subject for this condition
                                save_data =[save_location, tf_output_dir, filesep, subject(1:end-4),DatasetName,'_TF_baselinecorrected_', 'condition_',Conds{cond}];
                                %save(save_data, 'timefreqs_baselinecorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                parsave(save_data, ...
                                'timefreqs_baselinecorr', timefreqs_baselinecorr, ...
                                'frequency', frequency, ...
                                'ds_time', ds_time, ...
                                'channel_location', channel_location);
                                %fprintf('FINE HERE STILL');
                            end
                                
                        elseif BaselineCorrect == 0
                            
                            if Downsample ==0
                                %save out trial averaged and non-baseline corrected TF data
                                save_data =[save_location, tf_output_dir, filesep, subject(1:end-4),DatasetName,'_TF_nobaselinecorrection_', 'condition_',Conds{cond}];
                                %save(save_data, 'timefreqs_data', 'frequency', 'time','channel_location', '-v7.3');
                                parsave(save_data, ...
                                'timefreqs_data', timefreqs_data, ...
                                'frequency', frequency, ...
                                'time', time, ...
                                'channel_location', channel_location);
                            elseif Downsample ==1
                                %Downsample to 125Hz
                                timefreqs_data = timefreqs_data(:,1:2:size(timefreqs_data,2),:);
                                
                                %Downsample time variable
                                ds_time = downsample(time,2);
                                
                                %save out trial averaged and non-baseline corrected TF data
                                save_data =[save_location, tf_output_dir, filesep, subject(1:end-4),DatasetName,'_TF_nobaselinecorrection_', 'condition_',Conds{cond}];
                                %save(save_data, 'timefreqs_data', 'frequency', 'ds_time','channel_location', '-v7.3');
                                parsave(save_data, ...
                                'timefreqs_data', timefreqs_data, ...
                                'frequency', frequency, ...
                                'ds_time', ds_time, ...
                                'channel_location', channel_location);
                            end
                        end %end if you want baseline correction
                        
                        %Blank out power data over trials after saved to save memory space
                         timefreqs_data=[]; timefreqs_baselinecorr = [];
                        
                        %Inter-trial phase synchrony
                        if ITPS_calc == 1
                                cd(scripts_location)
                                % ITPS_script;
                                % INPUT
                                % data = single trial wavelet decomposed signal; frequency x time x trial x channel
                                data=[];
                                data = phase_data_strial_ch;
                                                                               
                                % OUTPUT: ITPS_all
                                %frequency x time x channels 
                                
                                %% Initialize output matrices depending on type of connectivity
                                ITPS_all       = zeros(size(data, 1), size(data, 2), size(data, 4));
                                         
                                %% ITPS Computations
                                %with subsampling
                                if Subsample == 1
                                    %initialize matrix with subsamples
                                    ITPS_ch = zeros(size(data,1), size(data,2), size(data,4));
                                    ITPS_subsamples = zeros(NumSubsamples, size(data,1), size(data,2), size(data,4));
                                    %Begin subsampling
                                    for samp=1:NumSubsamples
                                        fprintf('\n\n\n*** Calculating ITPS for subject %d (%s) subsample %d ***\n\n\n', sub, subject, samp);    
                                        for chan=1:size(data, 4)
                                            ITPS=[];
                                            for freq=1:size(data, 1)
                                                subtrials=[]; data_temp=[];
                                                %Get indices of trials for this subsample
                                                subtrials = randsample(1:size(data,3),NumTrialsPulled,false);
                                                %Index into those trials and pull them out for wPLI analyses
                                                data_temp = squeeze(data(freq,:,subtrials,chan));
                                                ITPS(freq,:) = abs(mean( exp(1i*angle(data_temp)),2));
                                            end %end loop through frequencies
                                            ITPS_ch(:,:,chan) = ITPS;
                                        end %end loop through channels
                                        %create matrix of subsamples for ITPS - samp x freq x time
                                        ITPS_subsamples(samp,:,:,:) = ITPS_ch;
                                    end %end loop through subsamples
                                    %average over subsamples - this is final ITPS
                                    ITPS_all = squeeze(mean(ITPS_subsamples,1));
                                     
                                    
                                % WITHOUT subsampling
                                elseif Subsample == 0
                                fprintf('\n\n\n*** Calculating ITPS for subject %d (%s) ***\n\n\n', sub, subject);    
                                    for chan=1:size(data, 4)
                                        ITPS=[];
                                        %loop through frequencies
                                        for freq=1:size(data, 1)
                                            % pull out data needed
                                            data_temp=[];
                                            data_temp = squeeze(data(freq,:,:,chan));
                                            ITPS(freq,:) = abs(mean( exp(1i*angle(data_temp)),2));
                                        end %end loop through frequencies
                                        ITPS_all(:,:,chan) = ITPS;
                                     end %end loop through channels
                                end %end loop through subsampling
                                
                                     %Baseline Correction
                                     if BaselineCorrect == 1
                                         
                                         %% baseline time indices
                                         basetimeidx   = dsearchn(EEG.times', BaselineTime');
                                         
                                         %Initialize ITPS_blncorr
                                         ITPS_blncorr = zeros(size(ITPS_all));
                                            
                                         for chanj=1:size(ITPS_all,3)
                                                for fi = 1:size(ITPS_all,1)
                                                    ITPS_blncorr(fi,:,chanj) = ( ITPS_all(fi,:,chanj) - mean(ITPS_all(fi,basetimeidx(1):basetimeidx(end),chanj)));
                                                end
                                         end
                                            
                                         if RestorEvent==1
                                             
                                             if Downsample==0
                                                 %save out trial averaged and baseline corrected ITPS data for this subject
                                                 save_data =[save_location, itps_output_dir, filesep, subject(1:end-4),DatasetName,'_ITPS_baselinecorrected'];
                                                 %save(save_data, 'ITPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ITPS_blncorr', ITPS_blncorr, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                             elseif Downsample==1
                                                 %Downsample
                                                 ITPS_blncorr = ITPS_blncorr(:,1:2:size(ITPS_blncorr,2),:);
                                                 
                                                 %Downsample time variable
                                                 ds_time = downsample(time,2);
                                                 
                                                 %save out trial averaged, downsampled, and baseline corrected ITPS data for this subject
                                                 save_data =[save_location, itps_output_dir, filesep, subject(1:end-4),DatasetName,'_ITPS_baselinecorrected'];
                                                 %save(save_data, 'ITPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ITPS_blncorr', ITPS_blncorr, ...
                                                 'frequency', frequency, ...
                                                 'ds_time', ds_time, ...
                                                 'channel_location', channel_location);
                                             end
                                         elseif RestorEvent==0
                                             if Downsample ==0
                                                 %save out trial averaged and baseline corrected ITPS data for this subject
                                                 save_data =[save_location, itps_output_dir, filesep, subject(1:end-4),DatasetName,'_ITPS_baselinecorrected_','condition_',Conds{cond}];
                                                 %save(save_data, 'ITPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ITPS_blncorr', ITPS_blncorr, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                             elseif Downsample ==1
                                                 %Downsample
                                                 ITPS_blncorr = ITPS_blncorr(:,1:2:size(ITPS_blncorr,2),:);
                                                 
                                                 %Downsample time variable
                                                 ds_time = downsample(time,2);
                                                 
                                                 %save out trial averaged and baseline corrected ITPS data for this subject
                                                 save_data =[save_location, itps_output_dir, filesep, subject(1:end-4),DatasetName,'_ITPS_baselinecorrected_','condition_',Conds{cond}];
                                                 %save(save_data, 'ITPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ITPS_blncorr', ITPS_blncorr, ...
                                                 'frequency', frequency, ...
                                                 'ds_time', ds_time, ...
                                                 'channel_location', channel_location);
                                             end
                                         end
                                         
                                     elseif BaselineCorrect == 0
                                         
                                         if RestorEvent==1
                                             if Downsample ==0
                                                 %save out trial averaged and baseline corrected ITPS data for this subject
                                                 save_data =[save_location, itps_output_dir, filesep, subject(1:end-4),DatasetName,'_ITPS_nobaselinecorrection'];
                                                 %save(save_data, 'ITPS_all', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ITPS_all', ITPS_all, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                             elseif Downsample == 1
                                                 %Downsample
                                                 ITPS_all = ITPS_all(:,1:2:size(ITPS_all,2),:);
                                                 
                                                 %Downsample time variable
                                                 ds_time = downsample(time,2);
                                                 
                                                 %save out trial averaged and baseline corrected ITPS data for this subject
                                                 save_data =[save_location, itps_output_dir, filesep, subject(1:end-4),DatasetName,'_ITPS_nobaselinecorrection'];
                                                 %save(save_data, 'ITPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ITPS_all', ITPS_all, ...
                                                 'frequency', frequency, ...
                                                 'ds_time', ds_time, ...
                                                 'channel_location', channel_location);
                                             end
                                         elseif RestorEvent==0 %if event-related, then save condition name in save file
                                             if Downsample ==0
                                                 %save out trial averaged and baseline corrected ITPS data for this subject
                                                 save_data =[save_location, itps_output_dir, filesep, subject(1:end-4),DatasetName,'_ITPS_nobaselinecorrection_','condition_',Conds{cond}];
                                                 %save(save_data, 'ITPS_all', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ITPS_all', ITPS_all, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                             elseif Downsample ==1
                                                 %Downsample
                                                 ITPS_all = ITPS_all(:,1:2:size(ITPS_all,2),:);
                                                 
                                                 %Downsample time variable
                                                 ds_time = downsample(time,2);
                                                 
                                                 %save out trial averaged, downsampled, and baseline corrected ITPS data for this subject
                                                 save_data =[save_location, itps_output_dir, filesep, subject(1:end-4),DatasetName,'_ITPS_nobaselinecorrection_','condition_',Conds{cond}];
                                                 %save(save_data, 'ITPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ITPS_all', ITPS_all, ...
                                                 'frequency', frequency, ...
                                                 'ds_time', ds_time, ...
                                                 'channel_location', channel_location);
                                             end %end if Downsampling should be done
                                         end %end if resting or event
                                     end %end if you want baseline correction
                                
                                
                                subtrials=[];
                                data_temp=[];
                                ITPS=[];
                                ITPS_ch=[];
                                ITPS_subsamples=[];
                                ITPS_all=[];
                                ITPS_blncorr=[];
                                ITPS_all=[]; %blank out matrix when done
                        else
                            %do nothing
                        end
                        
                        %Inter-channel calculations
                        if ICPS_calc == 1
                                if Compute_ICPS == 1 %calculate coherence
                                    % INPUT
                                    % data = single trial wavelet decomposed signal; for resting, frequency x time x trial x channel
                                                                                   % for event-related, condition x frequency x time x trial x channel
                                    data=[];
                                    data = phase_data_strial_ch;
                                                                                   
                                    % OUTPUT: ICPS_all
                                    %Resting, over time, all-to-all connectivity; frequency x trials x channel x channel
                                    %Resting, over trials, all-to-all connectivity; frequency x time x channel x channel
                                    %Resting, over time, seed-based; frequency x trials x non-seed channel 
                                    %Resting, over trials, seed-based; frequency x time x non-seed channel 
                                    
                                    %Event-Related, over time, all-to-all connectivity; condition x frequency x trials x channel x channel
                                    %Event-Related, over trials, all-to-all connectivity; condition x frequency x time x channel x channel
                                    %Event-Related, over time, seed-based; condition x frequency x trials x non-seed channel
                                    %Event-Related, over trials, seed-based; condition x frequency x time x non-seed channel

                                    % Convert all channel labels in EEG.chanlocs to uppercase.
                                    % EEG.chanlocs is a struct array where each entry corresponds to one
                                    % electrode. The .labels field contains the string name of that electrode.
                                    % This loop overwrites each label in place with its uppercase equivalent.
                                    % this is not needed for THRIVE (because numbered electrode labels)
                                    % and breaks the parfor loop, therefore commenting out
                                    % for i = 1:length(EEG.chanlocs)
                                    %     EEG.chanlocs(i).labels = upper(EEG.chanlocs(i).labels);
                                    % end
                                    
                                    % Convert the user-defined seed electrode name to uppercase.
                                    % Seed is a single string (e.g., 'Fz') defined earlier in the parameter
                                    % section. This ensures it will match the now-uppercased EEG.chanlocs.labels
                                    % when find(strcmp(...)) is called later to locate the seed channel index.
                                    % this is not needed for THRIVE (because numbered electrode labels)
                                    % and breaks the parfor loop, therefore commenting out
                                    % Seed = upper(Seed);
                                    
                                    % Convert all non-seed channel names in Elecs4Connect to uppercase.
                                    % Elecs4Connect is a cell array of strings defining which channels to
                                    % compute connectivity for (e.g., {'Fz','Cz','Pz'}).
                                    % cellfun applies upper() to each cell element individually.
                                    % 'UniformOutput', false is required because upper() returns a string
                                    % (non-scalar), so the output cannot be collapsed into a standard array.
                                    % this is not needed for THRIVE (because numbered electrode labels)
                                    % and breaks the parfor loop, therefore commenting out
                                    % Elecs4Connect = cellfun(@upper, Elecs4Connect, 'UniformOutput', false);
                                    
                                    %% Initialize output matrices depending on type of connectivity
                                    if TimeOrTrials == 0 %over time 
                                        if ConnectType == 0
                                            %Resting, over time, all-to-all
                                            ICPS_all       = zeros(size(data, 1), size(data, 3), size(data, 4), size(data, 4));
                                        elseif ConnectType == 1
                                            %Resting, over time, seed-based
                                            ICPS_all       = zeros(size(data, 1), size(data, 3), length(Elecs4Connect));
                                        end
                                    elseif TimeOrTrials ==1 %over trials
                                        if ConnectType == 0
                                            %Resting, over trials, all-to-all
                                            ICPS_all       = zeros(size(data, 1), size(data, 2), size(data, 4), size(data, 4));
                                        elseif ConnectType == 1
                                            %Resting, over trials, seed-based
                                            ICPS_all       = zeros(size(data, 1), size(data, 2), length(Elecs4Connect));
                                        end
                                    end
                                    
                                    %% Connectivity Computations
                                    %%Compute all-to-all connectivity
                                    if ConnectType == 0
                                        if TimeOrTrials==0
                                            fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) ***\n\n\n', sub, subject);
                                            for chani=1:size(data, 4)
                                                for chanj=chani:size(data, 4)
                                                    %% take cross-spectral density
                                                    crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                                                    %% ICPS
                                                    for freq=1:size(data, 1)
                                                        ICPS_all(freq,:,chani,chanj) = abs( mean( abs(crossspecden(freq,:,:)).*sign(crossspecden(freq,:,:)),2))./mean(abs(crossspecden(freq,:,:)),2);
                                                        ICPS_all(freq,:,chanj,chani) = abs( mean( abs(crossspecden(freq,:,:)).*sign(crossspecden(freq,:,:)),2))./mean(abs(crossspecden(freq,:,:)),2);
                                                    end %end loop through frequencies
                                                end %end second loop through channels
                                            end %end first loop through channels
                                            
                                        elseif TimeOrTrials==1
                                            %with subsampling
                                            if Subsample == 1
                                                ICPS_subsamples = zeros(NumSubsamples, size(data,1), size(data,2), size(data,4), size(data,4));
                                                %Begin subsampling
                                                for samp=1:NumSubsamples
                                                    fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) subsample %d ***\n\n\n', sub, subject, samp);
                                                    crossspecden=[]; subtrials=[]; crossspecden_temp=[]; ICPS_temp=[];
                                                    for chani=1:size(data, 4)
                                                        for chanj=chani:size(data, 4)
                                                            %% take cross-spectral density
                                                            crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                                                            %Get indices of trials for this subsample
                                                            subtrials = randsample(1:size(data,3),NumTrialsPulled,false);
                                                            %Index into those trials and pull them out for wPLI analyses
                                                            crossspecden_temp = crossspecden(:,:,subtrials);
                                                            for freq=1:size(data, 1)
                                                                ICPS_temp(freq,:,chani,chanj) = abs( mean( abs(crossspecden_temp(freq,:,:)).*sign(crossspecden_temp(freq,:,:)),3))./mean(abs(crossspecden_temp(freq,:,:)),3);
                                                                ICPS_temp(freq,:,chanj,chani) = abs( mean( abs(crossspecden_temp(freq,:,:)).*sign(crossspecden_temp(freq,:,:)),3))./mean(abs(crossspecden_temp(freq,:,:)),3);
                                                            end
                                                            
                                                        end %end second loop through channels
                                                    end %end first loop through channels
                                                    %create matrix of subsamples for wPLI - samp x freq x time x channel x channel
                                                    ICPS_subsamples(samp,:,:,:,:) = ICPS_temp;
                                                end %end loop through subsamples
                                                %average over subsamples - this is final ICPS
                                                ICPS_all = squeeze(mean(ICPS_subsamples,1));
                                                
                                                % WITHOUT subsampling
                                            elseif Subsample == 0
                                                fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) ***\n\n\n', sub, subject);
                                                for chani=1:size(data, 4)
                                                    for chanj=chani:size(data, 4)
                                                        %% take cross-spectral density
                                                        crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                                                        %% ICPS
                                                        for freq=1:size(data, 1)
                                                            ICPS_all(freq,:,chani,chanj) = abs( mean( abs(crossspecden(freq,:,:)).*sign(crossspecden(freq,:,:)),3))./mean(abs(crossspecden(freq,:,:)),3);
                                                            ICPS_all(freq,:,chanj,chani) = abs( mean( abs(crossspecden(freq,:,:)).*sign(crossspecden(freq,:,:)),3))./mean(abs(crossspecden(freq,:,:)),3);
                                                        end %end loop through frequencies
                                                    end %end second loop through channels
                                                end %end first loop through channels
                                            end %end if statement for subsample
                                        end %end over time or over trials
                                        
                                        %Baseline Correct, Downsample, and Save Data
                                        if TimeOrTrials == 0 %over time does not baseline correct or downsample
                                            if RestorEvent==1
                                                %save out trial averaged and baseline corrected ICPS data for this subject
                                                save_data =[save_location, icps_output_dir, filesep, subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection'];
                                                %save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ICPS_all', ICPS_all, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                            elseif RestorEvent==0
                                                %save out trial averaged and baseline corrected ICPS data for this subject
                                                save_data =[save_location, icps_output_dir, filesep,  subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection_','condition_',Conds{cond}];
                                                %save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ICPS_all', ICPS_all, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                            end
                                        elseif TimeOrTrials ==1
                                            %Baseline Correction
                                            if BaselineCorrect == 1
                                                
                                                %% baseline time indices
                                                basetimeidx   = dsearchn(EEG.times', BaselineTime');
                                                
                                                Baseline = squeeze(mean(ICPS_all(:,basetimeidx(1):basetimeidx(end),:,:),2));
                                                %loop through samples
                                                for t=1:size(ICPS_all,2)
                                                    ICPS_blncorr(:,t,:,:) = squeeze(ICPS_all(:,t,:,:)) - Baseline;
                                                end %end loop through frequencies
                                                
                                                if RestorEvent==1
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location, icps_output_dir, filesep,  subject(1:end-4),DatasetName,'_ICPS_baselinecorrected'];
                                                        %save(save_data, 'ICPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_blncorr', ICPS_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                                        ICPS_blncorr=ICPS_blncorr(:,1:2:size(ICPS_blncorr,2),:,:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location, icps_output_dir, filesep,  subject(1:end-4),DatasetName,'_ICPS_baselinecorrected'];
                                                        %save(save_data, 'ICPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_blncorr', ICPS_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                elseif RestorEvent==0 %if event-related, add condition name to saved file
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location, icps_output_dir, filesep,  subject(1:end-4),DatasetName,'_ICPS_baselinecorrected_','condition_',Conds{cond}];
                                                        %save(save_data, 'ICPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_blncorr', ICPS_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample==1
                                                        %Downsample
                                                        ICPS_blncorr=ICPS_blncorr(:,1:2:size(ICPS_blncorr,2),:,:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged, downsampled, and baseline corrected ICPS data for this subject
                                                        save_data =[save_location,  icps_output_dir, filesep, subject(1:end-4),DatasetName,'_ICPS_baselinecorrected_','condition_',Conds{cond}];
                                                        %save(save_data, 'ICPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_blncorr', ICPS_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                end
                                            elseif BaselineCorrect == 0
                                                if RestorEvent==1 %Resting
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location, icps_output_dir, filesep,  subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection'];
                                                        %save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_all', ICPS_all, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                                        ICPS_all = ICPS_all(:,1:2:size(ICPS_all,2),:,:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged, downsampled, baseline corrected ICPS data for this subject
                                                        save_data =[save_location,  icps_output_dir, filesep, subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection'];
                                                        %save(save_data, 'ICPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_all', ICPS_all, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                elseif RestorEvent==0 %Event-Related
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location,  icps_output_dir, filesep, subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection_','condition_',Conds{cond}];
                                                        %save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_all', ICPS_all, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                                        ICPS_all = ICPS_all(:,1:2:size(ICPS_all,2),:,:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location, icps_output_dir, filesep,  subject(1:end-4),DatasetName,'_ICPS_nobaselinecorrection_','condition_',Conds{cond}];
                                                        %save(save_data, 'ICPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_all', ICPS_all, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end %end if downsample
                                                end %end if resting or event
                                            end %end if you want baseline correction
                                        end %end if time or trials
                                        
                                        
                                        
                                        %SEED-BASED
                                    elseif ConnectType == 1
                                        if TimeOrTrials == 0
                                            fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) ***\n\n\n', sub, subject);
                                                
                                                % Kia's CHANGE: Added isempty() checks after each find() call for Seed_idx and
                                                % Elecs_idx. In the original, if a channel label was not found, find()
                                                % returned [] and the subsequent array assignment would either crash with
                                                % an unhelpful index error or silently assign nothing. The improved version
                                                % throws an explicit, informative error identifying which channel was missing.
                                                % Also replaced strcmpi() with strcmp() for Seed_idx lookup, since all
                                                % labels are now pre-normalized to uppercase at the top of the script.
                                        
                                                %Find index of seed electrode
                                                Seed_idx = find(strcmp({EEG.chanlocs.labels}, Seed));
                                                if isempty(Seed_idx)
                                                    error('Seed channel not found in dataset.');
                                                end
                                        
                                                Elecs_idx = zeros(1, length(Elecs4Connect));
                                                % Find indices of the non-seed channels
                                                for i = 1:length(Elecs4Connect)
                                                    idx = find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
                                                    if isempty(idx)
                                                        error('Channel "%s" from Elecs4Connect was not found in the dataset.', Elecs4Connect{i});
                                                    end
                                                    Elecs_idx(i) = idx;
                                                end
                                             %loop through channels
                                             for chanj=1:length(Elecs4Connect)
                                                % take cross-spectral density between two channels - one being the seed for ICPS
                                                crossspecden = squeeze(data(:,:,:,Seed_idx) .* conj(data(:,:,:,Elecs_idx(chanj))));
                                                %Initialize matrix
                                                ICPS_seed = zeros(size(data, 1), size(data, 3));
                                                %Loop through frequencies
                                                for freq=1:size(data, 1)
                                                    %wPLI for nonsubsampled seed-based over trials
                                                    ICPS_seed(freq,:) = abs(mean(exp(1i*angle(crossspecden(freq,:,:))),2));
                                                end %end loop through frequencies
                                                ICPS_all(:,:,chanj) = ICPS_seed;
                                            end % end loop through channels
                                        elseif TimeOrTrials==1
                                            %WITH subsampling
                                            if Subsample == 1
                                                %initialize matrix with subsamples
                                                ICPS_ch = zeros(size(data,1), size(data,2),length(Elecs4Connect));
                                                ICPS_subsamples = zeros(NumSubsamples, size(data,1), size(data,2),length(Elecs4Connect));
                                                %Start subsampling
                                                for samp=1:NumSubsamples
                                                    fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) subsample %d ***\n\n\n', sub, subject, samp);
                                                    crossspecden=[]; crossspecden_imag=[]; subtrials=[]; crossspecden_imag_temp=[]; 
                                                    
                                                        % Kia's CHANGE: Added isempty() checks after each find() call for Seed_idx and
                                                        % Elecs_idx. In the original, if a channel label was not found, find()
                                                        % returned [] and the subsequent array assignment would either crash with
                                                        % an unhelpful index error or silently assign nothing. The improved version
                                                        % throws an explicit, informative error identifying which channel was missing.
                                                        % Also replaced strcmpi() with strcmp() for Seed_idx lookup, since all
                                                        % labels are now pre-normalized to uppercase at the top of the script.
                                            
                                                        %Find index of seed electrode
                                                        Seed_idx = find(strcmp({EEG.chanlocs.labels}, Seed));
                                                        if isempty(Seed_idx)
                                                            error('Seed channel not found in dataset.');
                                                        end
                                            
                                                        Elecs_idx = zeros(1, length(Elecs4Connect));
                                                        % Find indices of the non-seed channels
                                                        for i = 1:length(Elecs4Connect)
                                                            idx = find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
                                                            if isempty(idx)
                                                                error('Channel "%s" from Elecs4Connect was not found in the dataset.', Elecs4Connect{i});
                                                            end
                                                            Elecs_idx(i) = idx;
                                                        end
                                                    
                                                    %loop through channels
                                                    for chanj=1:length(Elecs4Connect)
                                                        % take cross-spectral density between two channels - one being the seed for ICPS
                                                        crossspecden = squeeze(data(:,:,:,Seed_idx) .* conj(data(:,:,:,Elecs_idx(chanj)))); %
                                                        %subsampling the crossspecden - with replacement across subsample, but without replacement within subsamples
                                                        %Get indices of trials for this subsample
                                                        subtrials = randsample(1:size(data,3),NumTrialsPulled,false);
                                                        %Index into those trials and pull them out for wPLI analyses
                                                        crossspecden_temp = crossspecden(:,:,subtrials);
                                                        %Initialize matrix
                                                        ICPS_seed = zeros(size(data, 1), size(data, 2));
                                                        %Loop through frequencies
                                                        for freq=1:size(data, 1)
                                                            ICPS_seed(freq,:) = abs(mean(exp(1i*angle(crossspecden_temp(freq,:,:))),3));
                                                        end %end loop through frequencies
                                                        %create matrix with channelssamp x freq x time
                                                        ICPS_ch(:,:,chanj) = ICPS_seed;
                                                    end %end loop through channels
                                                    ICPS_subsamples(samp,:,:,:) = ICPS_ch;
                                                end%end loop through subsamples
                                                %average over subsamples - wPLI for subsampled seed-based connectivity
                                                ICPS_all = squeeze(mean(ICPS_subsamples,1));
                                                %SEED-BASED WITHOUT subsampling
                                            elseif Subsample == 0
                                                fprintf('\n\n\n*** Calculating ICPS for subject %d (%s) ***\n\n\n', sub, subject);
                                                % take cross-spectral density
                                                for chanj=1:length(Elecs4Connect)

                                                    % Kia's CHANGE: Added isempty() checks after each find() call for Seed_idx and
                                                    % Elecs_idx. In the original, if a channel label was not found, find()
                                                    % returned [] and the subsequent array assignment would either crash with
                                                    % an unhelpful index error or silently assign nothing. The improved version
                                                    % throws an explicit, informative error identifying which channel was missing.
                                                    % Also replaced strcmpi() with strcmp() for Seed_idx lookup, since all
                                                    % labels are now pre-normalized to uppercase at the top of the script.
                                        
                                                    %Find index of seed electrode
                                                    Seed_idx = find(strcmp({EEG.chanlocs.labels}, Seed));
                                                    if isempty(Seed_idx)
                                                        error('Seed channel not found in dataset.');
                                                    end
                                        
                                                    Elecs_idx = zeros(1, length(Elecs4Connect));
                                                    % Find indices of the non-seed channels
                                                    for i = 1:length(Elecs4Connect)
                                                        idx = find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
                                                        if isempty(idx)
                                                            error('Channel "%s" from Elecs4Connect was not found in the dataset.', Elecs4Connect{i});
                                                        end
                                                        Elecs_idx(i) = idx;
                                                    end
                                                    
                                                    % take cross-spectral density between two channels - one being the seed for ICPS
                                                    crossspecden = squeeze(data(:,:,:,Seed_idx) .* conj(data(:,:,:,Elecs_idx(chanj))));
                                                    %Initialize matrix
                                                    ICPS_seed = zeros(size(data,1), size(data,2));
                                                    for freq=1:size(data, 1)
                                                        %wPLI for nonsubsampled seed-based over trials
                                                        ICPS_seed(freq,:) = abs(mean(exp(1i*angle(crossspecden(freq,:,:))),3));
                                                    end %end loop through frequencies
                                                    ICPS_all(:,:,chanj) = ICPS_seed;
                                                end % end loop through channels
                                            end % end if subsampling statement
                                        end %end if TimeOrTrials
                                        
                                        %Baseline Correct, Downsample, and Save Data
                                        if TimeOrTrials == 0 %over time does not baseline correct or downsample
                                            if RestorEvent==1
                                                %save out trial averaged and baseline corrected ICPS data for this subject
                                                save_data =[save_location, icps_output_dir, filesep,  subject(1:end-4),DatasetName,'_ICPS_overtime_nobaselinecorrection'];
                                                %save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ICPS_all', ICPS_all, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                            elseif RestorEvent==0
                                                %save out trial averaged and baseline corrected ICPS data for this subject
                                                save_data =[save_location, icps_output_dir, filesep,  subject(1:end-4),DatasetName,'_ICPS_overtime_nobaselinecorrection_','condition_',Conds{cond}];
                                                %save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'ICPS_all', ICPS_all, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                            end
                                        elseif TimeOrTrials ==1
                                            %Baseline Correction
                                            if BaselineCorrect == 1
                                                
                                                %% baseline time indices
                                                basetimeidx   = dsearchn(EEG.times', BaselineTime');
                                                
                                                %Initialize ICPS_blncorr
                                                ICPS_blncorr = zeros(size(ICPS_all));
                                                
                                                %Baseline correct
                                                for chanj=1:size(ICPS_all,3)
                                                    for fi = 1:size(ICPS_all,1)
                                                        ICPS_blncorr(fi,:,chanj) = ( ICPS_all(fi,:,chanj) - mean(ICPS_all(fi,basetimeidx(1):basetimeidx(end),chanj)));
                                                    end
                                                end
                                    
                                                if RestorEvent==1
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location, icps_output_dir, filesep,  subject(1:end-4),DatasetName,'_ICPS_overtrials_baselinecorrected'];
                                                        %save(save_data, 'ICPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_blncorr', ICPS_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                                        ICPS_blncorr = ICPS_blncorr(:,1:2:size(ICPS_blncorr,2),:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location,  icps_output_dir, filesep, subject(1:end-4),DatasetName,'_ICPS_overtrials_baselinecorrected'];
                                                        %save(save_data, 'ICPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_blncorr', ICPS_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                elseif RestorEvent==0
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location,  icps_output_dir, filesep, subject(1:end-4),DatasetName,'_ICPS_overtrials_baselinecorrected_','condition_',Conds{cond}];
                                                        %save(save_data, 'ICPS_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_blncorr', ICPS_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample

                                                        % Kia's CHANGE (BUG FIX): the trailing indexing was corrected from ',:,:' to
                                                        % just ',:' because seed-based wPLI_blncorr is 3D (freq x time x channel),
                                                        % not 4D like the all-to-all matrix. Using ',:,:' on a 3D array would
                                                        % cause a MATLAB dimension mismatch error.
                                                        
                                                        ICPS_blncorr = ICPS_blncorr(:,1:2:size(ICPS_blncorr,2),:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location,  icps_output_dir, filesep, subject(1:end-4),DatasetName,'_ICPS_overtrials_baselinecorrected_','condition_',Conds{cond}];
                                                        %save(save_data, 'ICPS_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_blncorr', ICPS_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                end
                                            elseif BaselineCorrect == 0
                                                if RestorEvent==1 %resting
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location,  icps_output_dir, filesep, subject(1:end-4),DatasetName,'_ICPS_overtrials_nobaselinecorrection'];
                                                        %save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_all', ICPS_all, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                                        ICPS_all = ICPS_all(:,1:2:size(ICPS_all,2),:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged, downsampled, and baseline corrected ICPS data for this subject
                                                        save_data =[save_location,  icps_output_dir, filesep, subject(1:end-4),DatasetName,'_ICPS_overtrials_nobaselinecorrection'];
                                                        %save(save_data, 'ICPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_all', ICPS_all, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                elseif RestorEvent==0 %Event-Related
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location,  icps_output_dir, filesep, subject(1:end-4),DatasetName,'_ICPS_overtrials_nobaselinecorrection_','condition_',Conds{cond}];
                                                        parsave(save_data, ...
                                                        'ICPS_all', ICPS_all, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                        %save(save_data, 'ICPS_all', 'frequency', 'time','channel_location', '-v7.3');
                                                    elseif Downsample ==1
                                                        %Downsample
                                                        ICPS_all = ICPS_all(:,1:2:size(ICPS_all,2),:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged and baseline corrected ICPS data for this subject
                                                        save_data =[save_location, icps_output_dir, filesep,  subject(1:end-4),DatasetName,'_ICPS_overtrials_nobaselinecorrection_','condition_',Conds{cond}];
                                                        %save(save_data, 'ICPS_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'ICPS_all', ICPS_all, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end %end if downsampling
                                                end %end if resting or event
                                            end %end if you want baseline correction
                                        end %end if time or trials
                                        
                                    end %end if connectivity type loops
                                    
                                    ICPS_all=[];
                                    ICPS_blncorr=[];
                                    Baseline=[];
                                    crossspecden=[]; 
                                    subtrials=[];
                                    crossspecden_temp=[];
                                    ICPS_temp=[];
                                    ICPS_seed=[];
                                    ICPS_ch=[];
                                    ICPS_subsamples=[];                                    
                                    ICPS_all=[]; %blank out matrix when done
                                end

                                if Compute_wPLI == 1 % calculate wPLI
                                    % wPLI_script;
                                    % INPUT
                                    % data = single trial wavelet decomposed signal; frequency x time x trial x channel
                                    %The input will come from the output of the strial_timefreq_decomposition_script.m
                                    data=[];
                                    data = phase_data_strial_ch;
                                                                                   
                                    % OUTPUT: wPLI_all
                                    %over time, all-to-all connectivity; frequency x trials x channel x channel
                                    %over trials, all-to-all connectivity; frequency x time x channel x channel
                                    %over time, seed-based; frequency x trials x non-seed channel 
                                    %over trials, seed-based; frequency x time x non-seed channel 

                                    % Convert all channel labels in EEG.chanlocs to uppercase.
                                    % EEG.chanlocs is a struct array where each entry corresponds to one
                                    % electrode. The .labels field contains the string name of that electrode.
                                    % This loop overwrites each label in place with its uppercase equivalent.
                                    % this is not needed for THRIVE (because numbered electrode labels)
                                    % and breaks the parfor loop, therefore commenting out
                                    % for i = 1:length(EEG.chanlocs)
                                    %    EEG.chanlocs(i).labels = upper(EEG.chanlocs(i).labels);
                                    % end
                                    
                                    % Convert the user-defined seed electrode name to uppercase.
                                    % Seed is a single string (e.g., 'Fz') defined earlier in the parameter
                                    % section. This ensures it will match the now-uppercased EEG.chanlocs.labels
                                    % when find(strcmp(...)) is called later to locate the seed channel index.
                                    % this is not needed for THRIVE (because numbered electrode labels)
                                    % and breaks the parfor loop, therefore commenting out
                                    % Seed = upper(Seed);
                                    
                                    % Convert all non-seed channel names in Elecs4Connect to uppercase.
                                    % Elecs4Connect is a cell array of strings defining which channels to
                                    % compute connectivity for (e.g., {'Fz','Cz','Pz'}).
                                    % cellfun applies upper() to each cell element individually.
                                    % 'UniformOutput', false is required because upper() returns a string
                                    % (non-scalar), so the output cannot be collapsed into a standard array.
                                    % this is not needed for THRIVE (because numbered electrode labels)
                                    % and breaks the parfor loop, therefore commenting out
                                    % Elecs4Connect = cellfun(@upper, Elecs4Connect, 'UniformOutput', false);
                                    
                                    %% Initialize output matrices depending on type of connectivity
                                    if TimeOrTrials == 0 %over time
                                        if ConnectType == 0 % all-to-all
                                            %over time, all-to-all
                                            wPLI_all       = zeros(size(data, 1), size(data, 3), size(data, 4), size(data, 4));
                                        elseif ConnectType == 1 % seed-based
                                            %over time, seed-based
                                            wPLI_all       = zeros(size(data, 1), size(data, 3), length(Elecs4Connect));
                                        end
                                    elseif TimeOrTrials ==1 %over trials
                                        if ConnectType == 0 % all-to-all
                                            %over trials, all-to-all
                                            wPLI_all       = zeros(size(data, 1), size(data, 2), size(data, 4), size(data, 4));
                                        elseif ConnectType == 1 % seed-based
                                            %over trials, seed-based
                                            wPLI_all       = zeros(size(data, 1), size(data, 2), length(Elecs4Connect));
                                        end
                                    end
                                    
                                    %% Connectivity Computations
                                    %%Compute all-to-all connectivity
                                    if ConnectType == 0 % all-to-all
                                        if TimeOrTrials ==0 % over time
                                            fprintf('\n\n\n*** Calculating wPLI for subject %d (%s) ***\n\n\n', sub, subject);
                                            for chani=1:size(data, 4)
                                                for chanj=chani:size(data, 4)
                                                    %% take cross-spectral density
                                                    crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                                                    % take imaginary part of signal only
                                                    crossspecden_imag = imag(crossspecden);
                                                    %% weighted phase-lag index (shortcut, as implemented in Cohen's book)
                                                    for freq=1:size(data, 1)
                                                        wPLI_all(freq,:,chani,chanj) = abs( mean( abs(crossspecden_imag(freq,:,:)).*sign(crossspecden_imag(freq,:,:)),2))./mean(abs(crossspecden_imag(freq,:,:)),2);
                                                        wPLI_all(freq,:,chanj,chani) = abs( mean( abs(crossspecden_imag(freq,:,:)).*sign(crossspecden_imag(freq,:,:)),2))./mean(abs(crossspecden_imag(freq,:,:)),2);
                                                    end % end loop through frequencies 
                                                end % end second loop through channels 
                                            end % end first loop through channels
                                            
                                        elseif TimeOrTrials==1 % over trials
                                            %with subsampling
                                            if Subsample == 1
                                                %initialize matrix with subsamples
                                                wPLI_subsamples = zeros(NumSubsamples, size(data,1), size(data,2), size(data,4), size(data,4));
                                                %Begin subsampling
                                                for samp=1:NumSubsamples
                                                    fprintf('\n\n\n*** Calculating wPLI for subject %d (%s) subsample %d ***\n\n\n', sub, subject, samp);
                                                    crossspecden=[]; crossspecden_imag=[]; subtrials=[]; crossspecden_imag_temp=[]; weighted_phaselagidx_temp=[];
                                                    for chani=1:size(data, 4)
                                                        for chanj=chani:size(data, 4)
                                                            %% take cross-spectral density
                                                            crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                                                            % take imaginary part of signal only
                                                            crossspecden_imag = imag(crossspecden);
                                                            % Get indices of trials for this subsample
                                                            subtrials = randsample(1:size(data,3),NumTrialsPulled,false);
                                                            % Index into those trials and pull them out for wPLI analyses
                                                            crossspecden_imag_temp = crossspecden_imag(:,:,subtrials);
                                                            for freq=1:size(data, 1)
                                                                weighted_phaselagidx_temp(freq,:,chani,chanj) = abs( mean( abs(crossspecden_imag_temp(freq,:,:)).*sign(crossspecden_imag_temp(freq,:,:)),3))./mean(abs(crossspecden_imag_temp(freq,:,:)),3); 
                                                                weighted_phaselagidx_temp(freq,:,chanj,chani) = abs( mean( abs(crossspecden_imag_temp(freq,:,:)).*sign(crossspecden_imag_temp(freq,:,:)),3))./mean(abs(crossspecden_imag_temp(freq,:,:)),3); 
                                                            end % end loop through frequencies
                                                        end % end second loop through channels
                                                    end % end first loop through channels
                                                % create matrix of subsamples for wPLI - samp x freq x time
                                                wPLI_subsamples(samp,:,:,:,:) = weighted_phaselagidx_temp;
                                                end %end loop through subsamples
                                                %average over subsamples - this is final wPLI
                                                wPLI_all = squeeze(mean(wPLI_subsamples,1));
                                                
                                                % all-to-all connectivity WITHOUT subsampling
                                            elseif Subsample == 0
                                                fprintf('\n\n\n*** Calculating wPLI for subject %d (%s)***\n\n\n', sub, subject);
                                                for chani=1:size(data, 4)
                                                    for chanj=chani:size(data, 4)
                                                        %% take cross-spectral density
                                                        crossspecden = squeeze(data(:,:,:,chani) .* conj(data(:,:,:,chanj)));
                                                        % take imaginary part of signal only
                                                        crossspecden_imag = imag(crossspecden);
                                                        %% weighted phase-lag index (shortcut, as implemented in Cohen's book)
                                                        for freq=1:size(data, 1)
                                                            wPLI_all(freq,:,chani,chanj) = abs( mean( abs(crossspecden_imag(freq,:,:)).*sign(crossspecden_imag(freq,:,:)),3))./mean(abs(crossspecden_imag(freq,:,:)),3);
                                                            wPLI_all(freq,:,chanj,chani) = abs( mean( abs(crossspecden_imag(freq,:,:)).*sign(crossspecden_imag(freq,:,:)),3))./mean(abs(crossspecden_imag(freq,:,:)),3);
                                                        end %end loop through frequencies
                                                    end %end second loop through channels
                                                end %end first loop through channels
                                            end %end if subsampling
                                        end %end if time or trials
                                        
                                        %Baseline Correct, Downsample, and Save Data
                                        if TimeOrTrials == 0 %over time does not baseline correct or downsample
                                            if RestorEvent==1 % rest
                                                %save out trial averaged and no baseline corrected wPLI data for this subject
                                                save_data =[save_location,  wpli_output_dir, filesep, subject(1:end-4),DatasetName,'_wPLI_overtime_nobaselinecorrection'];
                                                %save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'wPLI_all', wPLI_all, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                            elseif RestorEvent==0 % event
                                                %save out trial averaged and no baseline corrected wPLI data for this subject
                                                save_data =[save_location,wpli_output_dir, filesep, subject(1:end-4),DatasetName,'_wPLI_overtime_nobaselinecorrection_','condition_',Conds{cond}];
                                                %save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'wPLI_all', wPLI_all, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                            end
                                        elseif TimeOrTrials ==1 % over trials
                                            %Baseline Correction
                                            if BaselineCorrect == 1
                                                
                                                %% baseline time indices
                                                basetimeidx   = dsearchn(EEG.times', BaselineTime');
                                    
                                                %Initialize wPLI_blncorr
                                                wPLI_blncorr = zeros(size(wPLI_all)); % Kia's CHANGE: pre-initialize before loop
                                                
                                                Baseline = squeeze(mean(wPLI_all(:,basetimeidx(1):basetimeidx(end),:,:),2));
                                                %loop through samples 
                                                for t=1:size(wPLI_all,2)
                                                    wPLI_blncorr(:,t,:,:) = squeeze(wPLI_all(:,t,:,:)) - Baseline;
                                                end %end loop through frequencies
                                                
                                                if RestorEvent==1 % rest
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected'];
                                                        %save(save_data, 'wPLI_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_blncorr', wPLI_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                    
                                                        % Kia's CHANGE (BUG FIX): The original referenced 'wPLI_baselinecorr' which is
                                                        % an undefined variable (the correct name used throughout the script is
                                                        % 'wPLI_blncorr'). This would have caused a MATLAB runtime error whenever
                                                        % baseline correction + downsampling was selected for all-to-all over-trials.
                                    
                                                        wPLI_blncorr=wPLI_blncorr(:,1:2:size(wPLI_blncorr,2),:,:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected'];
                                                        %save(save_data, 'wPLI_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_blncorr', wPLI_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                elseif RestorEvent==0 %if event-related, add condition name to saved file
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected_','condition_',Conds{cond}];
                                                        %save(save_data, 'wPLI_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_blncorr', wPLI_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample==1
                                                        %Downsample
                                                        wPLI_blncorr=wPLI_blncorr(:,1:2:size(wPLI_blncorr,2),:,:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged, downsampled, and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected_','condition_',Conds{cond}];
                                                        %save(save_data, 'wPLI_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_blncorr', wPLI_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                end
                                            elseif BaselineCorrect == 0
                                                if RestorEvent==1 %Resting
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection'];
                                                        %save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_all', wPLI_all, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                                        wPLI_all = wPLI_all(:,1:2:size(wPLI_all,2),:,:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged, downsampled, baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection'];
                                                        %save(save_data, 'wPLI_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_all', wPLI_all, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                elseif RestorEvent==0 %Event-Related
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection_','condition_',Conds{cond}];
                                                        %save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_all', wPLI_all, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                                        wPLI_all = wPLI_all(:,1:2:size(wPLI_all,2),:,:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection_','condition_',Conds{cond}];
                                                        %save(save_data, 'wPLI_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_all', wPLI_all, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end %end if downsampling
                                                end %end if resting or event
                                            end %end if you want baseline correction
                                        end %end over time or trials                                        
%SEED-BASED
                                    elseif ConnectType == 1 
                                        if TimeOrTrials== 0 % over time
                                            fprintf('\n\n\n*** Calculating wPLI for subject %d (%s) ***\n\n\n', sub, subject);
                                            
                                            % Kia's CHANGE: Added isempty() checks after each find() call for Seed_idx and
                                            % Elecs_idx. In the original, if a channel label was not found, find()
                                            % returned [] and the subsequent array assignment would either crash with
                                            % an unhelpful index error or silently assign nothing. The improved version
                                            % throws an explicit, informative error identifying which channel was missing.
                                            % Also replaced strcmpi() with strcmp() for Seed_idx lookup, since all
                                            % labels are now pre-normalized to uppercase at the top of the script.
                                    
                                            %Find index of seed electrode
                                            Seed_idx = find(strcmp({EEG.chanlocs.labels}, Seed));
                                            if isempty(Seed_idx)
                                                error('Seed channel not found in dataset.');
                                            end
                                    
                                            Elecs_idx = zeros(1, length(Elecs4Connect));
                                            % Find indices of the non-seed channels
                                            for i = 1:length(Elecs4Connect)
                                                idx = find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
                                                if isempty(idx)
                                                    error('Channel "%s" from Elecs4Connect was not found in the dataset.', Elecs4Connect{i});
                                                end
                                                Elecs_idx(i) = idx;
                                            end
                                    
                                            % take cross-spectral density
                                            for chanj=1:length(Elecs4Connect)
                                                % take cross-spectral density between two channels - one being the seed for wPLI
                                                crossspecden = squeeze(data(:,:,:,Seed_idx) .* conj(data(:,:,:,Elecs_idx(chanj))));
                                                % take imaginary part of signal only
                                                crossspecden_imag = imag(crossspecden);
                                                %initialize matrix
                                                weighted_phaselagidx_seed = zeros(size(data, 1), size(data, 3));
                                                for freq=1:size(data, 1)
                                                    %wPLI for nonsubsampled seed-based over time
                                                    
                                                    % Kia's CHANGE: The original seed-based formula used abs(mean(exp(1i*angle(...))))
                                                    % which computes PLI (phase-lag index), not wPLI (weighted PLI).
                                                    % The correct wPLI formula weights each observation by the magnitude of
                                                    % the imaginary cross-spectral density: it is the mean of |imag| * sign(imag)
                                                    % divided by the mean of |imag|. This matches Cohen (2014) and is now
                                                    % consistent with the all-to-all computation branches above.
                                    
                                                    weighted_phaselagidx_seed(freq,:) = abs( mean( abs(crossspecden_imag(freq,:,:)).*sign(crossspecden_imag(freq,:,:)),2))./mean(abs(crossspecden_imag(freq,:,:)),2);
                                               
                                                end %end loop through frequencies
                                                wPLI_all(:,:,chanj) = weighted_phaselagidx_seed;
                                            end % end loop through channels
                                        elseif TimeOrTrials==1 % over trials
                                            %subsampling seed-based connectivity
                                            if Subsample == 1
                                                %initialize subsamples matrix
                                                wPLI_subsamples = zeros(NumSubsamples, size(data,1), size(data,2),length(Elecs4Connect));
                                                
                                                % Kia's CHANGE: Added isempty() checks after each find() call for Seed_idx and
                                                % Elecs_idx. In the original, if a channel label was not found, find()
                                                % returned [] and the subsequent array assignment would either crash with
                                                % an unhelpful index error or silently assign nothing. The improved version
                                                % throws an explicit, informative error identifying which channel was missing.
                                                % Also replaced strcmpi() with strcmp() for Seed_idx lookup, since all
                                                % labels are now pre-normalized to uppercase at the top of the script.
                                    
                                                %Find index of seed electrode
                                                Seed_idx = find(strcmp({EEG.chanlocs.labels}, Seed));
                                                if isempty(Seed_idx)
                                                    error('Seed channel not found in dataset.');
                                                end
                                    
                                                Elecs_idx = zeros(1, length(Elecs4Connect));
                                                % Find indices of the non-seed channels
                                                for i = 1:length(Elecs4Connect)
                                                    idx = find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
                                                    if isempty(idx)
                                                        error('Channel "%s" from Elecs4Connect was not found in the dataset.', Elecs4Connect{i});
                                                    end
                                                    Elecs_idx(i) = idx;
                                                end
                                    
                                                        
                                                %Start subsampling
                                                for samp=1:NumSubsamples
                                                    fprintf('\n\n\n*** Calculating wPLI for subject %d (%s) subsample %d ***\n\n\n', sub, subject, samp);
                                                    crossspecden=[]; crossspecden_imag=[]; subtrials=[]; crossspecden_imag_temp=[];
                                                    %Initialize matrix
                                                    wPLI_ch = zeros(size(data,1), size(data,2),length(Elecs4Connect));
                                                    %loop through channels
                                                    for chanj=1:length(Elecs4Connect)
                                                        %Take cross spectral density
                                                        crossspecden = squeeze(data(:,:,:,Seed_idx) .* conj(data(:,:,:,Elecs_idx(chanj))));
                                                        % take imaginary part of signal only
                                                        crossspecden_imag = imag(crossspecden);
                                                        %subsampling the crossspecden - with replacement across subsample, but without replacement within subsamples
                                                        %Get indices of trials for this subsample
                                                        subtrials = randsample(1:size(data,3),NumTrialsPulled,false);
                                                        %Index into those trials and pull them out for wPLI analyses
                                                        crossspecden_imag_temp = crossspecden_imag(:,:,subtrials);
                                                        %Initialize matrix
                                                        weighted_phaselagidx_seed = zeros(size(data, 1), size(data, 2));
                                                        %Loop through frequencies
                                                        for freq=1:size(data, 1)
                                                            % wPLI
                                                            
                                                            % Kia's CHANGE: The original seed-based formula used abs(mean(exp(1i*angle(...))))
                                                            % which computes PLI (phase-lag index), not wPLI (weighted PLI).
                                                            % The correct wPLI formula weights each observation by the magnitude of
                                                            % the imaginary cross-spectral density: it is the mean of |imag| * sign(imag)
                                                            % divided by the mean of |imag|. This matches Cohen (2014) and is now
                                                            % consistent with the all-to-all computation branches above.
                                    
                                                            weighted_phaselagidx_seed(freq,:) = abs( mean( abs(crossspecden_imag_temp(freq,:,:)).*sign(crossspecden_imag_temp(freq,:,:)),3))./mean(abs(crossspecden_imag_temp(freq,:,:)),3);
                                                        
                                                        end % end loop through frequencies 
                                                        wPLI_ch(:,:,chanj) = weighted_phaselagidx_seed;
                                                    end % end loop through channels
                                                    wPLI_subsamples(samp,:,:,:) = wPLI_ch;
                                                end % end loop through subsamples
                                                %average over subsamples - wPLI for subsampled seed-based connectivity
                                                wPLI_all = squeeze(mean(wPLI_subsamples,1));
                                                
                                                %SEED-BASED WITHOUT subsampling
                                            elseif Subsample == 0
                                                fprintf('\n\n\n*** Calculating wPLI for subject %d (%s) ***\n\n\n', sub, subject);
                                                
                                                
                                                % Kia's CHANGE: Added isempty() checks after each find() call for Seed_idx and
                                                % Elecs_idx. In the original, if a channel label was not found, find()
                                                % returned [] and the subsequent array assignment would either crash with
                                                % an unhelpful index error or silently assign nothing. The improved version
                                                % throws an explicit, informative error identifying which channel was missing.
                                                % Also replaced strcmpi() with strcmp() for Seed_idx lookup, since all
                                                % labels are now pre-normalized to uppercase at the top of the script.
                                    
                                    
                                                %Find index of seed electrode
                                                Seed_idx = find(strcmp({EEG.chanlocs.labels}, Seed));
                                                if isempty(Seed_idx)
                                                    error('Seed channel not found in dataset.');
                                                end
                                    
                                                Elecs_idx = zeros(1, length(Elecs4Connect));
                                                % Find indices of the non-seed channels
                                                for i = 1:length(Elecs4Connect)
                                                    idx = find(strcmp({EEG.chanlocs.labels}, Elecs4Connect{i}));
                                                    if isempty(idx)
                                                        error('Channel "%s" from Elecs4Connect was not found in the dataset.', Elecs4Connect{i});
                                                    end
                                                    Elecs_idx(i) = idx;
                                                end
                                                
                                                %loop through channels
                                                for chanj=1:length(Elecs4Connect)
                                                    % take cross-spectral density between two channels - one being the seed for wPLI
                                                    crossspecden = squeeze(data(:,:,:,Seed_idx) .* conj(data(:,:,:,Elecs_idx(chanj))));
                                                    % take imaginary part of signal only
                                                    crossspecden_imag = imag(crossspecden);
                                                    %Initialize matrix
                                                    weighted_phaselagidx_seed = zeros(size(data, 1), size(data, 2));
                                                    %Loop through frequencies
                                                    for freq=1:size(data, 1)
                                                        %wPLI for nonsubsampled seed-based over trials
                                                        
                                                        % Kia's CHANGE: The original seed-based formula used abs(mean(exp(1i*angle(...))))
                                                        % which computes PLI (phase-lag index), not wPLI (weighted PLI).
                                                        % The correct wPLI formula weights each observation by the magnitude of
                                                        % the imaginary cross-spectral density: it is the mean of |imag| * sign(imag)
                                                        % divided by the mean of |imag|. This matches Cohen (2014) and is now
                                                        % consistent with the all-to-all computation branches above.
                                    
                                                        weighted_phaselagidx_seed(freq,:) = abs( mean( abs(crossspecden_imag(freq,:,:)).*sign(crossspecden_imag(freq,:,:)),3))./mean(abs(crossspecden_imag(freq,:,:)),3);
                                                    
                                                    end %end loop through frequencies
                                                    wPLI_all(:,:,chanj) = weighted_phaselagidx_seed;
                                                end % end loop through channels
                                            end % end if subsampling statement
                                        end % end if time or trials
                                        
                                        % Baseline Correct, Downsample, and Save Data
                                        if TimeOrTrials == 0 %over time does not baseline correct or downsample
                                            if RestorEvent==1
                                                %save out trial averaged and baseline corrected wPLI data for this subject
                                                save_data =[save_location,wpli_output_dir, filesep, subject(1:end-4),DatasetName,'_wPLI_overtime_nobaselinecorrection'];
                                                %save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'wPLI_all', wPLI_all, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                            elseif RestorEvent==0
                                                %save out trial averaged and baseline corrected wPLI data for this subject
                                                save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtime_nobaselinecorrection_','condition_',Conds{cond}];
                                                %save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                                                 parsave(save_data, ...
                                                 'wPLI_all', wPLI_all, ...
                                                 'frequency', frequency, ...
                                                 'time', time, ...
                                                 'channel_location', channel_location);
                                            end
                                        elseif TimeOrTrials ==1
                                            % Baseline Correction
                                            if BaselineCorrect == 1
                                                
                                                %% baseline time indices
                                                basetimeidx   = dsearchn(EEG.times', BaselineTime');
                                                
                                                %Initialize wPLI_blncorr
                                                wPLI_blncorr = zeros(size(wPLI_all));
                                                
                                                %Baseline Correct
                                               for chanj=1:size(wPLI_all,3)
                                                    for fi = 1:size(wPLI_all,1)
                                                        wPLI_blncorr(fi,:,chanj) = ( wPLI_all(fi,:,chanj) - mean(wPLI_all(fi,basetimeidx(1):basetimeidx(end),chanj)));
                                                    end
                                                end
                                                
                                                if RestorEvent==1
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location,wpli_output_dir, filesep, subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected'];
                                                        %save(save_data, 'wPLI_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_blncorr', wPLI_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                    
                                                        % Kia's CHANGE (BUG FIX): Same wPLI_baselinecorr → wPLI_blncorr typo fix as
                                                        % above. Additionally, the trailing indexing was corrected from ',:,:' to
                                                        % just ',:' because seed-based wPLI_blncorr is 3D (freq x time x channel),
                                                        % not 4D like the all-to-all matrix. Using ',:,:' on a 3D array would
                                                        % cause a MATLAB dimension mismatch error.
                                    
                                                        wPLI_blncorr=wPLI_blncorr(:,1:2:size(wPLI_blncorr,2),:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected'];
                                                        %save(save_data, 'wPLI_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_blncorr', wPLI_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                elseif RestorEvent==0 %if event-related, add condition name to saved file
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected_','condition_',Conds{cond}];
                                                        %save(save_data, 'wPLI_blncorr', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_blncorr', wPLI_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample==1
                                                        %Downsample
                                                        wPLI_blncorr=wPLI_blncorr(:,1:2:size(wPLI_blncorr,2),:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged, downsampled, and baseline corrected wPLI data for this subject
                                                        save_data =[save_location,wpli_output_dir, filesep, subject(1:end-4),DatasetName,'_wPLI_overtrials_baselinecorrected_','condition_',Conds{cond}];
                                                        %save(save_data, 'wPLI_blncorr', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_blncorr', wPLI_blncorr, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                end
                                            elseif BaselineCorrect == 0
                                                if RestorEvent==1 %Resting
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection'];
                                                        %save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_all', wPLI_all, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                                        wPLI_all = wPLI_all(:,1:2:size(wPLI_all,2),:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged, downsampled, baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection'];
                                                        %save(save_data, 'wPLI_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_all', wPLI_all, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end
                                                elseif RestorEvent==0 %Event-Related
                                                    if Downsample ==0
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection_','condition_',Conds{cond}];
                                                        %save(save_data, 'wPLI_all', 'frequency', 'time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_all', wPLI_all, ...
                                                        'frequency', frequency, ...
                                                        'time', time, ...
                                                        'channel_location', channel_location);
                                                    elseif Downsample ==1
                                                        %Downsample
                                                        wPLI_all = wPLI_all(:,1:2:size(wPLI_all,2),:);
                                                        
                                                        %Downsample time variable
                                                        ds_time = downsample(time,2);
                                                        
                                                        %save out trial averaged and baseline corrected wPLI data for this subject
                                                        save_data =[save_location, wpli_output_dir, filesep,subject(1:end-4),DatasetName,'_wPLI_overtrials_nobaselinecorrection_','condition_',Conds{cond}];
                                                        %save(save_data, 'wPLI_all', 'frequency', 'ds_time','channel_location', '-v7.3');
                                                        parsave(save_data, ...
                                                        'wPLI_all', wPLI_all, ...
                                                        'frequency', frequency, ...
                                                        'ds_time', ds_time, ...
                                                        'channel_location', channel_location);
                                                    end %end if downsampling
                                                end %end if resting or event
                                            end %end if you want baseline correction
                                        end %end if time or trials
                                    end %end if connectivity type loops
                                    wPLI_all=[]; %blank out matrix when done
                                end
                        else
                            %do nothing
                        end
                        
            end %end if there is enough data for this condition
        end %end loop through conditions
    end %end if rest or event
%TrialNums_cell{sub} = local_TrialNums;
%all_TrialNums = vertcat(TrialNums_cell{:});
TrialNums(sub,:) = local_TrialNums;
%send(D, subject);
end
%delete(waitbar);
save_data = [save_location 'TrialNums_' datestr(now, 'mm_dd_yyyy_HH_MM_SS') '.mat'] 
parsave(save_data, 'TrialNums', TrialNums, 'subject_list', subject_list, 'Conds', Conds)    
%Save out table with trial numbers
%TrialNums_table = struct2table(TrialNums);
%TrialNums = [TrialNums_cell{:}];
%writetable(TrialNums_table,[save_location 'TrialNums.csv']);
%writetable(struct2table(all_TrialNums), [save_location 'TrialNums.csv']);
