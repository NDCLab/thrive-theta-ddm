
%clear % clear matlab workspace
cluster = parcluster('local');
%% Setting up other things
parpool(cluster, str2num(getenv('SLURM_CPUS_PER_TASK'))) % this should be same as --cpus-per-task
pool = gcp('nocreate');  % Get the current parallel pool without creating a new one
if isempty(pool)
    disp('No parallel pool is currently running.');
else
    disp(['Parallel pool with ', num2str(pool.NumWorkers), ' workers is running.']);
end
%Location of MADE and ADJUSTED-ADJUST scripts
% addpath(genpath([main_dir filesep 'MADE-EEG-preprocessing-pipeline']));% enter the path of the EEGLAB folder in this line
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeg_preprocessing'));% enter the path of the folder in this line

addpath(genpath('/home/data/NDClab/analyses/thrive-theta-ddm/code/matlab/'))

%Location of "EEG
% addpath(genpath([main_dir filesep 'eeglab13_4_4b']));% enter the path of the EEGLAB folder in this line
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b'));% enter the path of the EEGLAB folder in this line

%remove path to octave functions inside matlab to prevent errors when
% rmpath([main_dir filesep 'eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])
rmpath(['/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])

%% setup; run this section before any other section below
%modifying above, to account for files named differently
%specify parameters of data to process
task = 'all';
procStage = 'processed_data';
%session = 's1_r1'; 

% Derivations
visitDirName = session; %visit folder does not list "e1"
visitFileName = [session, '_e1']; %file names include "e1" designation

%location of analysis folder
analysis_dir = '/home/data/NDClab/analyses/thrive-theta-ddm';

%location of dataset folder
dataset_dir = '/home/data/NDClab/datasets/thrive-dataset';

% Setting up other things

% 1. Enter the path of the folder that has the data to be analyzed
data_location = [dataset_dir filesep 'derivatives' filesep 'preprocessed'];

% 2. Enter the path of the folder where you want to save the postprocessing outputs
save_location = [analysis_dir filesep 'derivatives' filesep 'preprocessed/csd_data' filesep visitDirName];
% 3. this is the correct channel location file BUT INCORRECT PATH!

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
if exist(save_location, 'dir') == 0
    mkdir(save_location);
end

for site = 1:64;
    trodes{site} = num2str(site);
end
Montage_64=ExtractMontage('/home/data/NDClab/analyses/thrive-theta-ddm/code/preprocessing-eeg/64ch_bv_montage_csd.csd', trodes');
[G, H] = GetGH(Montage_64);

diary(sprintf('csd_log_%s.log', datestr(now, 'mm-dd-yyyy_HH_MM_SS')))

% loop through each participant in the study
parfor subject = 1:length(datafile_names)
    try
	% extract participant number
        subNumText = datafile_names{subject}(5:11);
        output_file_path = fullfile(save_location, datafile_names{subject});
        
        % Check if the file already exists in the save location
        if exist(output_file_path, 'file')
            fprintf('Subject %s: File already processed. Skipping...\n', subNumText);
            continue;
        end
	%load the original data set
	EEG = pop_loadset('filename', datafile_names{subject}, 'filepath', datafile_paths{subject});
	EEG = eeg_checkset(EEG);
        
	EEG = pop_selectevent(EEG, 'latency','-.1 <= .1','deleteevents','on');
	EEG = eeg_checkset(EEG);
	fprintf('Subject %s: Processing %d events\n', subNumText, length(EEG.event));
        
        % Perform CSD transformation on each epoch
        data = zeros(size(EEG.data), 'single');
	for ne = 1:length(EEG.epoch)
            myEEG = single(EEG.data(:, :, ne));
	    MyResults = CSD(myEEG, G, H);            % compute CSD for <channels-by-samples> 2-D epoch
	    data(:, :, ne) = MyResults;
	end

	EEG.data = data;
	EEG = eeg_checkset(EEG);
	disp(size(EEG.event))
	    
	EEG = pop_editset(EEG, 'setname', datafile_names{subject});
	EEG = pop_saveset(EEG, 'filename', datafile_names{subject}, 'filepath', save_location);
        fprintf('Subject %s: Successfully processed and saved\n', subNumText);
%	clear data;
    catch ME
        fprintf('Error processing subject %s: %s\n', subNumText, ME.message);
        continue;
    end
end
