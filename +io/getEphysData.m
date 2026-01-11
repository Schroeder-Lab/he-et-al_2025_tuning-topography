function data = getEphysData(folder)
%GETEPHYSDATA   Load spiking data. Only consider single- and multi-unit
%activity.

% INPUTS
% folder        path to data of recording session

% OUTPUTS
% data
%   .times      [k], times of spikes of all considered units
%   .amps       [k], amplitudes of spikes of all considered units
%   .clusters   [k], cluster ID of spikes of all considered units
%   .depths     [k], depths of spikes of all considered units
%   .clusterIDs [units], IDs of considered units as returned after manual
%               spike sorting
%   .clusterDepths  [units], median depth of all spikes of each considered
%               unit

data.times = readNPY(fullfile(folder, "spikes.times.npy"));
data.amps = readNPY(fullfile(folder, "spikes.amps.npy"));
data.clusters = readNPY(fullfile(folder, "spikes.clusters.npy"));
data.depths = readNPY(fullfile(folder, "spikes.depths.npy"));
data.clusterIDs = readNPY(fullfile(folder, "clusters.ids.npy"));
if isfile(fullfile(folder, "clusters.depths.npy"))
    data.clusterDepths = readNPY(fullfile(folder, "clusters.depths.npy"));
end