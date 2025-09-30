function [units, stimTypes] = selectRFStim(results, minEV, minPeak)

units = [];
EVs = zeros(0, 2);
n = 0;
for stimType = 1:2
    if isempty(results{stimType})
        continue
    end
    good = find(results{stimType}.EV >= minEV & ...
        results{stimType}.peaks >= minPeak);
    [member, ind] = ismember(results{stimType}.units(good), units);
    EVs(ind(member), stimType) = ...
        results{stimType}.EV(good(member));
    units = [units; results{stimType}.units(good(~member))];
    EVs(n+1:length(units), stimType) = ...
        results{stimType}.EV(good(~member));
    n = n + sum(~member);
end
[~, stimTypes] = max(EVs, [], 2);