function [t, newTraces] = alignSampling(time, traces, delays, delayIDs)

dt = median(diff(time));
dd = median(diff(delays));
upsample = round(dt / dd);
newDt = dt / upsample;
t = reshape((time + (0:upsample-1) .* newDt)', [], 1);
newTraces = NaN(length(t), size(traces,2));
for d = 1:length(delays)
    ind = delayIDs == d;
    newTraces(:,ind) = interp1(time + delays(d), traces(:,ind), t, ...
        'linear', NaN);
end