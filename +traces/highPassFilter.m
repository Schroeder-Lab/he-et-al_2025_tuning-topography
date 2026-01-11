function traces = highPassFilter(traces, time, smoothWin)
%HIGHPASSFILTER   Highpass traces by subtracting smoothed traces using 
%moving median.

% INPUTS
% traces    [t x ROIs], traces of units
% time      [t], sample times of traces
% smoothWin double, window size of moving median in seconds

% OUTPUTS
% traces    [t x ROIs], filtered traces

dt = median(diff(time));
winSamples = round(smoothWin / dt);
smoothed = smoothdata(traces, 1, "movmedian", winSamples, "omitnan");
traces = traces - smoothed;