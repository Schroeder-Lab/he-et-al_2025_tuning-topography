function data = getLFPMetaData(folder)
%GETLFPMETADATA   Load LFP meta-data.

% INPUTS
% folder            path to data of recording session

% OUTPUTS
% data
%   .samplingRate   double, sampling rate of LFP signal
%   .numChans       double, number of recorded channels
%   .Imax           double, dynamic range of signal
%   .Vmax           double, maximum voltage covered by dynamic range
%   .lfpGain        double, gain during recording

opts = detectImportOptions(folder, ...
    "FileType", "text", "ReadVariableNames", false, ...
    "NumHeaderLines", 0);
opts = setvartype(opts,'string');
meta = readtable(folder, opts);
data.samplingRate = str2double(meta.Var2( ...
    strcmp(meta.Var1, "imSampRate")));
data.numChans = str2double(meta.Var2( ...
    strcmp(meta.Var1, "nSavedChans")));
data.Imax = str2double(meta.Var2( ...
    strcmp(meta.Var1, "imMaxInt")));
data.Vmax = str2double(meta.Var2( ...
    strcmp(meta.Var1, "imAiRangeMax")));
lfpGain = split(meta.Var2(strcmp(meta.Var1, "~imroTbl")), ")(");
lfpGain = split(lfpGain(2));
data.lfpGain = str2double(lfpGain(5));