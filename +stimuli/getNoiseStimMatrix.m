function [stimMatrix, edges, gridX, gridY] = ...
    getNoiseStimMatrix(edges, frames, stimOrder)
%GETNOISESTIMMATRIX  Get stimulus of every presented frame; ignore
%ipsilateral part

% INPUTS
% edges         [left right top bottom], location of stimulus boundaries in
%               visual space, locations above horizon <0
% frames        [uniqueStimuli x pixelRows x pixelColumns], gray-value
%               of each pixel for each unique stimulus (-1: black, 0: gray,
%               1: white)
% stimOrder     [k], indices of frames in order of presentation

% OUTPUTS
% stimMatrix    [k x pixelRows x pixelColumns], gray-value
%               of each pixel for each presented stimulus frame (-1: black, 
%               0: gray, 1: white)
% edges         [left right top bottom], location of stimulus boundaries in
%               visual space, locations above horizon >0
% gridX         Horizontal (along azimuth) locations of pixel centres
% gridY         Vertical (along elevation) locations of pixel centres

edges = double([edges([1 2]), -edges([3 4])]);
gridW = diff(edges(1:2)) / size(frames,3);
gridH = -diff(edges(3:4)) / size(frames,2);
% ignore pixels in ipsilateral (right) hemifield
if edges(1) * edges(2) < 0
    % determine right edge of all pixel columns
    rightEdges = edges(1) + (1:size(frames,3)) .* gridW;
    validPix = find(rightEdges <= 0);
    frames = frames(:,:,validPix);
    edges(2) = rightEdges(validPix(end));
end
stimMatrix = frames(stimOrder,:,:);
gridX = linspace(edges(1) + 0.5 * gridW, edges(2) - 0.5 * gridW, ...
    size(stimMatrix, 3));
gridY = linspace(edges(3) - 0.5 * gridH, edges(4) + 0.5 * gridH, ...
    size(stimMatrix, 2));