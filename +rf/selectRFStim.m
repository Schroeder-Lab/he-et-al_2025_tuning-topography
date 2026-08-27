function stimTypes = selectRFStim(results, minEV, minPeak, dist2edge)
%SELECTRFSTIM   For each unit, select stimulus (visual noise or cirlces) 
% resulting in higher explained variance.

% INPUTS
% results   {2}, entries for two stimulus types
%   .gaussPars      [ROIs x parameters], (amplitude, xCenter, xSTD,
%                   yCenter, ySTD, rotation, offset)
%   .peaks  [units], amplitude of RF peak (of best type) in STD of
%           noise
%   .EVs    [units], explained variance of trace prediction based on
%           Gaussian map (gaussMasks) * temporal weights
%           (timeWeights)
% minEV     double, minimum explained variance for significant RF
% minPeak   double, minimum peak for significant RF
% dist2edge double, minimum distance of RF centre to monitor edge

% OUTPUTS
% stimTypes [units], 1: if 1st stimulus is better, 2: if 2nd stimulus is
%           better, NaN: if no stimulus yielded significant RF

if nargin < 4
    dist2edge = 0;
end

n = 0;
for stimType = 1:2
    if isempty(results{stimType})
        continue
    end
    if n == 0
        n = length(results{stimType}.peaks);
        EVs = NaN(n, 2);
    end
    dt = results{stimType};
    edges = getResultEdges(dt);
    pos = results{stimType}.gaussPars(:,[2 4]);
    % exclude all RFs too close to monitor edge
    invalidRF = dt.EV < minEV | dt.peaks < minPeak | ...
        pos(:,1) < edges(1)+dist2edge | ...
        pos(:,1) > edges(2)-dist2edge | ...
        pos(:,2) > edges(3)-dist2edge | ...
        pos(:,2) < edges(4)+dist2edge;
    EVs(~invalidRF, stimType) = results{stimType}.EV(~invalidRF);
end
[~, stimTypes] = max(EVs, [], 2);
stimTypes(all(isnan(EVs), 2)) = NaN;
end

function edges = getResultEdges(result)
if isfield(result, 'edges')
    edges = result.edges;
    return
end
if ~isfield(result, 'x') || ~isfield(result, 'y')
    error('RF result must contain either edges or x/y coordinates.')
end
gridX = result.x;
gridY = result.y;
gridW = median(diff(gridX));
gridH = median(-diff(gridY));
% edges: [left right top bottom] (above horizon: >0)
edges = [gridX(1)-0.5*gridW gridX(end)+0.5*gridW ...
    gridY(1)+0.5*gridH gridY(end)-0.5*gridH];
end
