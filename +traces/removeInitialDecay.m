function tracesNew = removeInitialDecay(traces, time, decayWin, ...
    decayThresh)
%REMOVEINITIALDECAY   Fit and remove exponential decay from trace.

% INPUTS
% traces        [t x ROIs], traces of units
% time          [t], sample times of traces
% decayWin      double, duration at start over which decay is considered
% decayThresh   double, threshold of minimum decay in terms of STDs of
%               trace

% OUTPUTS
% tracesNew     [t x ROIs], corrected traces

tracesNew = traces;
timeBin = median(diff(time));
% identify traces with decay
indUnits = find(mean(traces(1:round(decayWin / timeBin),:),1,'omitnan') > ...
    mean(traces,1,'omitnan') + decayThresh .* std(traces,0,1,'omitnan'));
for iUnit = 1:length(indUnits)
    y = traces(:, indUnits(iUnit));
    y = fillmissing(y, 'linear');
    % fit double exponential to trace
    f = fit((1:length(y))', y, ...
        @(a,b,c,d,e,x) a + b .* exp(-x ./ c) + d .* exp(-x ./ e), ...
        'Lower', [0 0 0 0 0], ...
        'Upper', [max(y) max(y) 500 max(y) 500], ...
        'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
    % remove fit
    tracesNew(:, indUnits(iUnit)) = y - f(1:size(traces,1)) + f.a;
end