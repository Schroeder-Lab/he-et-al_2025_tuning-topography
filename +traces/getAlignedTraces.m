function [aligned, t] = getAlignedTraces(trace, time, eventTimes, window)
%GETALIGNEDTRACES   Algin continuous traces to event times.

% INPUTS
% trace         [t x ROIs], traces of units
% time          [t], sample times of traces
% eventTimes    [events], event times
% window        [1 or events x 2], start and end of window relative to each
%               event time, specify one if window is the same for all
%               events, or speficy one window per event

% OUTPUTS
% aligned       [t_a x events x ROIs], aligned traces
% t             [t_a], time of aligned traces relative to event times

% initialize output
aligned = [];
t = [];

% get window samples (tInd) relative to event time
dt = median(diff(time));
w = round(window / dt);
mini = min(w(:,1));
maxi = max(w(:,2));
tmp = mini : maxi;
if size(window,1) == 1
    tInd = repmat(tmp, length(eventTimes), 1);
else
    if size(window,1) ~= length(eventTimes)
        disp('No. events and windows have to be the same!')
        return
    end
    tInd = repmat(tmp, size(window,1), 1);
    for j = 1:size(window,1)
        tInd(j, tmp<w(j,1) | tmp>w(j,2)) = NaN;
    end
end
% sample times relative to event times
t = tmp .* dt;

% collect aligned trace snippets
aligned = NaN(length(t), length(eventTimes), size(trace, 2));
for ev = 1:length(eventTimes)
    % find trace index of event start
    evInd = find(time > eventTimes(ev), 1);
    if isempty(evInd)
        continue
    end
    if evInd>1 && eventTimes(ev)-time(evInd-1) < time(evInd)-eventTimes(ev)
        evInd = evInd - 1;
    end
    % find window indices within trace
    inds = evInd + tInd(ev,:);
    % only consider samples within trace
    valid = inds > 0 & inds <= length(time);
    aligned(valid,ev,:) = trace(inds(valid),:);
end