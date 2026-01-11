function [stimMatrix, x, y, z] = ...
    getCirclesStimMatrix(xPos, yPos, diameter, isWhite)
%GETCIRCLESSTIMMATRIX   Get stimulus of every frame in circle stimulus
%paradigm.

% INPUTS        
% xPos          [t], azimuth of presented circle in each stimulus frame
% yPos          [t], elevation of presented circle in each stimulus frame
% diameter      [t], diameter of presented circle in each stimulus frame
% isWhite       [t], gray value of presented circle in each stimulus frame,
%               0: black, 1: white

% OUTPUTS
% stimMatrix    [t x rows x columns x sizes], position, size and gray value
%               (-1: black, 1: white) of circle in each stimulus frame
% x             [rows], list of unique azimuth positions of circles
% y             [columns], list of unique elevation positions of circles
% z             [sizes], list of unique diameters of circles

% determine tested circle positions and sizes
x = unique(xPos); % left to right
y = flip(unique(yPos)); % top to bottom
z = unique(diameter);

% check whether all combinations of positions and sizes were presented; if
% not, set stimMatrix for those combinations to NaN
allCombis = [ ...
    reshape(repmat(x(:), 1, length(y), length(z)), [], 1), ...
    reshape(repmat(y(:)', length(x), 1, length(z)), [], 1), ...
    reshape(repmat(permute(z(:),[2 3 1]), ...
    length(x), length(y), 1), [], 1)];
notTested = setdiff(allCombis, unique([xPos, yPos, diameter], 'rows'), 'rows');

% for each time point (1st dim), make grid of all tested positions (2nd + 
% 3rd dim), and all tested sizes (4th dim); 1: white circle, -1: black
% circle
stimMatrix = zeros(length(xPos), length(y), length(x), length(z));
% column indices for each time
c = xPos' == x;
indC = repmat((1:length(x))', 1, length(xPos));
c = indC(c);
% row indices for each time
r = yPos' == y;
indR = repmat((1:length(y))', 1, length(xPos));
r = indR(r);
% diameter indices for each time
d = diameter' == z;
indD = repmat((1:length(z))', 1, length(xPos));
d = indD(d);
% get linear indices of circle positions and sizes for each time
ind = sub2ind(size(stimMatrix), (1:length(xPos))', r, c, d);
stimMatrix(ind(isWhite == 1)) = 1;
stimMatrix(ind(isWhite ~= 1)) = -1;

for j = 1:size(notTested,1)
    stimMatrix(:,notTested(2),notTested(1),notTested(3)) = NaN;
end