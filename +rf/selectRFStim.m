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

% OUTPUTS
% stimTypes [units], 1: if 1st stimulus is better, 2: if 2nd stimulus is
%           better, NaN: if no stimulus yielded significant RF

edges = [-95 0 51 -27]; % [left, right, top, bottom] where bottom
                        % and top are negative if above horizon

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