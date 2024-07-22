function [direction, length] = vectorAveraging(amplitudes, stimDirections)
%VECTORAVERAGING   Determine angle and length of average vector.

% INPUT
% amplitudes        [repetitions x stim], response amplitudes to multiple
%                   repetitions of the same stimuli
% stimDirections    [stim x 1], direction of stimuli in degrees

% OUTPUTS
% direction         double, directin of average vector in degrees
% length            double, length of average vector (input vectors
%                   normalized by sum)

% determine average response across repetitions
medAmps = median(amplitudes,1,'omitnan')';
% if most responses are negative, multiply by -1
if sum(medAmps) < 0
    medAmps = -medAmps;
end
% ignore negative responses
medAmps(medAmps < 0) = 0;
% normalize responses to sum of 1
medAmps = medAmps ./ sum(medAmps);
% translate stimulus directions to radians
dirsRadian = deg2rad(stimDirections);

% create vectors representing average responses to each stimulus direction
vectors = medAmps .* exp(1i .* dirsRadian);
% average vectors
meanVector = mean(vectors);

% direction of average vector in degrees
direction = rad2deg(angle(meanVector));
% length of average vector
length = abs(meanVector);