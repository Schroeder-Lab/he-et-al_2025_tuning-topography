function data = getGratingInfo(folder, type)
%GETGRATINGINFO   Load stimulus data.

% INPUTS
% folder            path to data of recording session
% type              str, stimulus type:
%                   'grating'/'gratingsDrifting'/'gratingsStatic'/'bars'

% OUTPUTS
% data
%   .times          [trials x 2], time of stimulus onsets (1st column) and
%                   offsets (2nd column)
%   .ids            [trials x 1], ID of stimulus presented in each trial
%   .interval       [2 x 1], start and end time of stimulus paradigm
%   (.directions)   (only for gratings, gratingsDrifting, or bars) 
%                   [stim x 1], direction of movement for each stimulus
%   (.orientations) (only for gratingsStatic) 
%                   [stim x 1], orientation for each stimulus
%   (.phases)       (only for gratingsStatic) 
%                   [stim x 1], spatial phase for each stimulus

if nargin < 2
    type = 'grating';
end

% return empty array if stimulus data not available
if isfile(fullfile(folder, sprintf('_ss_%s.intervals.npy', type)))
    % stimulus onset and offset times, in sec
    data.times = readNPY(fullfile(folder, ...
        sprintf('_ss_%s.intervals.npy', type)));
else
    data = [];
    return
end
% stimulus IDs
data.ids = readNPY(fullfile(folder, sprintf('_ss_%s._ss_%sID.npy', type, type)));
% start and end of recording for specific stimulus paradigm
data.interval = readNPY(fullfile(folder, ...
    sprintf('_ss_recordings.%s_intervals.npy', type)));
% for each stimulus: direction or [orientation, spatial phase (for static gratings)]
switch type
    case {'gratingsDrifting', 'bars'}
        data.directions = double(readNPY(fullfile(folder, ...
            sprintf('_ss_%sID.directions.npy', type))));
    case 'gratingsStatic'
        data.orientations = readNPY(fullfile(folder, ...
            sprintf('_ss_%sID.orientations.npy', type)));
        data.phases = readNPY(fullfile(folder, ...
            sprintf('_ss_%sID.spatialPhases.npy', type)));
    otherwise
        data.directions = readNPY(fullfile(folder, ...
            sprintf('_ss_%sID.directions.npy', type)));
end