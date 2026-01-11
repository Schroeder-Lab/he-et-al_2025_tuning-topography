function distance = determineDistance(x, y)
%DETERMINEDISTANCE   Calculate all pairwise distances in brain.

% INPUTS
% x         [ROIs x 1], horizontal position in imaging plane
% y         [ROIs x 1], vertical position in imaging plane

% OUTPUTS
% distance  [pairs x 1], pairwise distances

distance = sqrt((x-x').^2 + (y-y').^2);
distance = distance(tril(true(size(distance)),-1));