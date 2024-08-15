function [tracesNew, timeNew] = alignPlaneTraces(traces, time, delays, planes)
%ALIGNPLANETRACES   Upsample and align traces of different imaging planes
%to a common time axis.

% INPUTS
% trace         [t x ROIs], traces of units
% time          [t], sample times of traces
% delays        [planes], time delay of each imaging plane relative to
%               first plane
% planes        [ROIs], planeID of each unit

% OUTPUTS
% tracesNew     [t_new x ROIs], aligned traces
% time_new      [t_new], new time axis

% upsample time trace to include delays of all planes
timeBin = median(diff(time));
delta_t = median(diff(delays));
upsample = round(timeBin / delta_t);
timeBin = timeBin / upsample;
timeNew = reshape((time + (0:upsample-1) * timeBin)', [], 1);
% interpolate traces to new time axis; take care of NaN values
tracesNew = NaN(length(timeNew), size(traces,2));
for d = 1:length(delays)
    indUnits = find(planes == d & ~all(isnan(traces),1)');
    for n = indUnits'
        nanInd1 = isnan(traces(:,n));
        tracesNew(:,n) = interp1(time(~nanInd1) + delays(d), ...
            traces(~nanInd1,n), timeNew, 'pchip');
        nanInd2 = reshape(repmat(nanInd1, 1, upsample)', [], 1);
        tracesNew(nanInd2,n) = NaN;
    end
end