function stimTypes = selectRFStim(results, minEV, minPeak)

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