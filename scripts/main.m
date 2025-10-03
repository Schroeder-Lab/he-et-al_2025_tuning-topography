%% Folder definitions
% Define locations of data, plots, and code.
% folders.data:     raw data; subfolders: [neurons/boutons]/[animal]/[date]
% folders.tools:    code repositories
% folders.repo:     current repository
% folders.plots:    to save plots

% % Mac
% folders.data = '/Users/ss2302/Library/CloudStorage/OneDrive-UniversityofSussex/Projects/2023_OrientationColumns/DataToPublish';
% folders.results = '/Users/ss2302/Library/CloudStorage/OneDrive-UniversityofSussex/Projects/2023_OrientationColumns/Results_Sylvia';
% folders.tools = '/Users/ss2302/dev/toolboxes';
% folders.repo = '/Users/ss2302/dev/he_schroeder_columns';
% folders.plots = '/Users/ss2302/Library/CloudStorage/OneDrive-UniversityofSussex/Projects/2023_OrientationColumns/Results_Sylvia';

% % PC on campus (Sussex Desk)
% folders.data = 'C:\Users\Sylvia\OneDrive - University of Sussex\Projects\2023_OrientationColumns\DataToPublish';
% folders.dataRawEphys = 'Z:\RawData';
% % folders.results = 'D:\Results\OrientationColumns';
% folders.results = 'C:\Users\Sylvia\OneDrive - University of Sussex\Projects\2023_OrientationColumns\Results_Sylvia';
% folders.tools = 'C:\dev\toolboxes';
% folders.repo = 'C:\dev\workspaces\he_schroeder_columns';
% folders.plots = 'C:\Users\Sylvia\OneDrive - University of Sussex\Projects\2023_OrientationColumns\Results_Sylvia';

% PC on campus (The Machine)
folders.data = 'C:\Users\sylvi\OneDrive - University of Sussex\Projects\2023_OrientationColumns\DataToPublish';
folders.dataRawEphys = 'Z:\RawData';
% folders.results = 'D:\Results\OrientationColumns';
folders.results = 'C:\Users\sylvi\OneDrive - University of Sussex\Projects\2023_OrientationColumns\Results_Sylvia';
folders.tools = 'C:\dev\toolboxes';
folders.repo = 'C:\dev\workspaces\SchroederLab\he_schroeder_columns';
folders.plots = 'C:\Users\sylvi\OneDrive - University of Sussex\Projects\2023_OrientationColumns\Results_Sylvia';

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
% addpath(genpath(fullfile(folders.tools, 'CircStat2012a')))
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
main_determineSCdepth(folders);

%% Determine depth of each unit relative to SC layers
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
Figure01(folders)
Figure01S(folders)

% Pairwise difference in tuning preferences depending on horizontal distance
Figure02(folders)
Figure02S(folders)

% Receptive fields and pairwise tuning differences depending on RF distance
Figure03(folders)
Figure03S(folders)

% Global distribution of direction and orientation preferences across
% visual field
Figure04(folders)
Figure04S(folders)

% Match between direction/orientation preference and expected preference at
% RF location (as predicted by longitudinal and latitudinal geometry)
Figure05(folders)
Figure05S(folders)

% Direction and orientation preferences across SC depth
Figure06(folders)
Figure06S(folders)