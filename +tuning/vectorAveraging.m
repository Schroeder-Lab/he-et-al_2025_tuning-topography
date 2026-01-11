function [direction, length, respSign] = vectorAveraging(amplitudes, stimDirections)
%VECTORAVERAGING   Determine angle and length of average vector.

% INPUT
% amplitudes        [repetitions x stim], response amplitudes to multiple
%                   repetitions of the same stimuli
% stimDirections    [stim x 1], direction of stimuli in degrees

% OUTPUTS
% direction         double, directin of average vector in degrees
% length            double, length of average vector (input vectors
%                   normalized by sum)
% respSign          1 or -1, sign of mean response amplitude (across all
%                   stimuli)

% determine mean response across repetitions
meanAmps = mean(amplitudes,1,'omitnan')';
respSign = 1;
% if most responses are negative, multiply by -1; this means that the unit
% is mostly suppressed by the visual stimuli
if sum(meanAmps) < 0
    meanAmps = -meanAmps;
    respSign = -1;
end
% ignore negative responses
meanAmps(meanAmps < 0) = 0;
% normalize responses to sum of 1
meanAmps = meanAmps ./ sum(meanAmps);
% translate stimulus directions to radians
dirsRadian = deg2rad(stimDirections);

% create vectors representing average responses to each stimulus direction
vectors = meanAmps .* exp(1i .* dirsRadian);
% sum vectors
finalVector = sum(vectors);

% direction of final vector in degrees
direction = mod(rad2deg(angle(finalVector)), 360);
% length of final vector
length = abs(finalVector);