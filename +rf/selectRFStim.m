function stimTypes = selectRFStim(results, minEV, minPeak)
%SELECTRFSTIM   For each unit, select stimulus (visual noise or cirlces) 
% resulting in higher explained variance.

% INPUTS
% results   {2}, entries for two stimulus types
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

n = 0;
for stimType = 1:2
    if isempty(results{stimType})
        continue
    end
    if n == 0
        n = length(results{stimType}.peaks);
        EVs = NaN(n, 2);
    end
    good = results{stimType}.EV >= minEV & ...
        results{stimType}.peaks >= minPeak;
    EVs(good, stimType) = results{stimType}.EV(good);
end
[~, stimTypes] = max(EVs, [], 2);
stimTypes(all(isnan(EVs), 2)) = NaN;