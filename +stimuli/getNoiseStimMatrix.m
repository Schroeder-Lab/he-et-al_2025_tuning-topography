function [stimMatrix, edges, gridX, gridY] = ...
    getNoiseStimMatrix(edges, frames, stimOrder)

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