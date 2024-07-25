function distance = determineDistance(x, y)
%DETERMINEDISTANCE   Calculate all pairwise distances in brain.

% INPUTS
% x         [ROIs x 1], horizontal position in brain
% y         [ROIs x 1], vertical position in brain

% OUTPUTS
% distance  [pairs x 1], pairwise distances

distance = sqrt((x-x').^2 + (y-y').^2);
distance = distance(tril(true(size(distance)),-1));