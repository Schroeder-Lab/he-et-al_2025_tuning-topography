function [time, stim] = fillGapsInNoiseStim(t_stim, stim)

% find time gaps in stimulus presentation (usually when same visual noise
% stimulus was repeated several times)
stimBin = median(diff(t_stim));
indGap = find(diff(t_stim) > 2 * stimBin);
time = t_stim;
% fill gaps with zeros in stimulus matrix
if ~isempty(indGap)
    stim_new = [];
    time_new = [];
    k = 1;
    for g = 1:length(indGap)
        add = floor(diff(t_stim(indGap(g) + [0 1])) ./ stimBin);
        stim_new = [stim_new; stim(k:indGap(g),:); ...
            zeros(add, size(stim,2))];
        time_new = [time_new; time(k:indGap(g)); ...
            time(indGap(g)) + (1:add)' .* stimBin];
        k = indGap(g)+1;
    end
    stim_new = [stim_new; stim(k:end,:)];
    time_new = [time_new; time(k:end)];
    stim = stim_new;
    time = time_new;
end