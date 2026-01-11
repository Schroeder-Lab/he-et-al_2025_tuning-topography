%% Folder definitions
% Define locations of data, plots, and code.
% folders.data:     pre-processed data; subfolders: [neurons/boutons]/[animal]/[date]
% folders.dataRawEphys: raw electrophysiology data (if available)
% folders.tools:    code repositories (npy-matlab, circstat-matlab)
% folders.repo:     current code repository
% folders.results:  to save intermediate results
% folders.plots:    to save plots

folders.data = 'C:\DataToPublish';
folders.dataRawEphys = 'Z:\RawData';
folders.tools = 'C:\dev\toolboxes';
folders.repo = 'C:\dev\he-et-al_2025_tuning-topography';
folders.results = 'C:\Results';
folders.plots = 'C:\Plots';

%% Global parameters
glob.figPositionDefault = [680 460 560 420];

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(genpath(fullfile(folders.tools, 'circstat-matlab')))
addpath(genpath(fullfile(folders.repo)))

%% 2P data

%% Fit response kernels (gratings + bars)
% Given the pre-processed calcium traces, fit response kernels, response
% amplitudes (per trial), and kernel shift (per stimulus, in case of bar 
% stimuli) to calcium responses to drifting gratings, static gratings, and
% moving bars.
main_fitStimResponseKernels(folders)

%% Determine direction and orientation selectivity
% Given the fitted response kernels (from main_fitStimResponseKernels.m),
% use vector averaging to determine preferred direction and orientation 
% (direction of final vector) as well as direction and orientation 
% selectivity (length of final vector). Use permutation test to determine
% significance of direction and orientation selectivity.
main_calculateDS_OS(folders)

%% Plot direction and orientation maps (at brain position)
% Plot the location of all units per recorded plane (using ROI masks) and
% color code preferred orientation or direction.
main_plotDirectionOrientationMapsPerPlane(folders)

%% Plot preference difference depending on brain distance
% Plot pairwise differences in preferred direction or orientation against
% the distance of ROIs in the brain. Consider only ROIs within the same
% plane or ROIs across all planes (ignoring distances in depth).
main_plotBrainDistanceVersusPreferenceDifference(folders)

%% Map receptive fields
% Given the calcium responses to the visual noise stimulus, fit spatial
% receptive fields (ON and OFF fields).
main_fitReceptiveFields(folders)

%% Fit retinotopy (relatioship: brain position - RF position)
% Find mapping between RF position and brain position (retinotopy). Then
% infer RF position of units where RF could not be mapped.
main_mapRFposToRetinotopy(folders)

%% Ephys data

%% Determine SC depth along probe
% Determine surface of SC based on visually evoked LFP, and SGS-SO border
% of SC from current-source-density
main_determineSCdepth(folders);

%% Determine depth of each unit relative to SC layers
% Given the site of the SC surface and its SGS-SO border, determine depth
% of each unit within SC (postive values -> microns below surface).
main_determineUnitDepths(folders);

%% Map receptive fields
% Given the spike responses to the visual noise or circle stimulus, fit spatial
% receptive fields (ON and OFF fields).
main_fitReceptiveFields_ephys(folders)

%% Determine direction and orientation selectivity
% Based on average firing rate during grating presentations,
% use vector averaging to determine preferred direction and orientation 
% (direction of final vector) as well as direction and orientation 
% selectivity (length of final vector). Use permutation test to determine
% significance of direction and orientation selectivity.
main_calculateDS_OS_ephys(folders)

%% Figures
% Distribution of direction & orientation preferences, and direction &
% orientation selectivity
Figure01(folders, glob)
Figure01S(folders, glob)

% Receptive fields 
Figure02(folders, glob)
Figure02S(folders)

% Match of preferences with topographic predictions
Figure03(folders, glob)
Figure03S(folders, glob)

% Pairwise difference in tuning preferences depending on horizontal distance
Figure04(folders, glob)
Figure04S(folders, glob)

% Global distribution of direction and orientation preferences across
% visual field
Figure05(folders, glob)
Figure05S(folders)

% Direction and orientation preferences across SC depth
Figure06(folders, glob)
Figure06S(folders, glob)

% Receptive fields and match of preferences with topographic predictions
Figure07(folders, glob)
