function [newTraces, newTime] = alignSampling(traces, time, delayIDs, delays)
%ALIGNSAMPLING   Upsample and align traces of different imaging planes
%to a common time axis.

% INPUTS
% traces        [t x ROIs], traces of units
% time          [t], sample times of traces
% delayIDs      [ROIs], planeID of each unit
% delays        [planes], time delay of each imaging plane relative to
%               first plane

% OUTPUTS
% newTraces     [t_new x ROIs], aligned traces
% newTime       [t_new], new time axis

if length(delays) <= 1
    newTraces = traces;
    newTime = time;
    return
end

% upsample time trace to include delays of all planes
timeBin = median(diff(time));
delta_t = median(diff(delays));
upsample = round(timeBin / delta_t);
timeBin = timeBin / upsample;
newTime = reshape((time + (0:upsample-1) * timeBin)', [], 1);
% interpolate traces to new time axis; take care of NaN values
newTraces = NaN(length(newTime), size(traces,2));
for d = 1:length(delays)
    indUnits = find(delayIDs == d);
    for n = indUnits'
        if sum(isnan(traces(:,n)))/size(traces,1) > 0.8
            continue
        end
        nanInd1 = isnan(traces(:,n));
        newTraces(:,n) = interp1(time(~nanInd1) + delays(d), ...
            traces(~nanInd1,n), newTime, 'pchip');
        nanInd2 = reshape(repmat(nanInd1, 1, upsample)', [], 1);
        newTraces(nanInd2,n) = NaN;
    end
end