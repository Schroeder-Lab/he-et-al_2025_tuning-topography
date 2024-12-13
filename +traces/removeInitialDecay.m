function newTraces = removeInitialDecay(traces, time, window_decay, ...
    threshold)
%REMOVEINITIALDECAY   Fit and remove exponential decay from trace.

% INPUTS
% traces        [t x ROIs], traces of units
% time          [t], sample times of traces
% decayWin      double, duration at start over which decay is considered
% decayThresh   double, threshold of minimum decay in terms of STDs of
%               trace

% OUTPUTS
% tracesNew     [t x ROIs], corrected traces

newTraces = traces;
timeBin = median(diff(time));
% identify traces with decay
indUnits = find(mean(traces(1:round(window_decay / timeBin),:),1,'omitnan') > ...
    mean(traces,1,'omitnan') + threshold .* std(traces,0,1,'omitnan'));
for iUnit = 1:length(indUnits)
    y = traces(:, indUnits(iUnit));
    y = fillmissing(y, 'linear');
    % fit double exponential to trace
    f = fit((1:length(y))', y, ...
        @(a,b,c,d,e,x) a + b .* exp(-x ./ c) + d .* exp(-x ./ e), ...
        'Lower', [0 0 0 0 0], ...
        'Upper', [max(y) max(y) length(y) max(y) length(y)], ...
        'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
    % remove fit
    newTraces(:, indUnits(iUnit)) = y - f(1:size(traces,1)) + f.a;
end