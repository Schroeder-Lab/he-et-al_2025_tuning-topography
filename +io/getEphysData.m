function data = getEphysData(folder)
%GETEPHYSDATA   Load ephys data. Only consider single- and multi-unit
%activity.

% INPUTS
% folder        path to data of recording session
% properties    .samplingRate (Hz)

% OUTPUTS
% data
%   .times
%   .amps
%   .clusters
%   .depths
%   .clusterIDs

data.times = readNPY(fullfile(folder, "spikes.times.npy"));
data.amps = readNPY(fullfile(folder, "spikes.amps.npy"));
data.clusters = readNPY(fullfile(folder, "spikes.clusters.npy"));
data.depths = readNPY(fullfile(folder, "spikes.depths.npy"));
data.clusterIDs = readNPY(fullfile(folder, "clusters.ids.npy"));
if isfile(fullfile(folder, "clusters.depths.npy"))
    data.clusterDepths = readNPY(fullfile(folder, "clusters.depths.npy"));
end