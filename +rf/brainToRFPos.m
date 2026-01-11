function rfPos = brainToRFPos(brain, a1, b1, c1, a2, k2)
%BRAINTORFPOS   Determine location of RF in visual space given unit's
%location in imaging plane and fitted model parameters.

% INPUTS
% brain     [2*ROIs], first half: horizontal position of units in FOV,
%           second half: vertical position of units in FOV
% a1        double, model parameter
% b1        double, model parameter
% c1        double, model parameter
% a2        double, model parameter
% k2        double, model parameter

% OUTPUTS
% rfPos     [2*ROIs], first half: RF azimuth of units,
%           second half: RF elevation of units

nUnits = floor(length(brain) / 2);
brainX = brain(1 : nUnits);
brainY = brain(nUnits+1 : 2*nUnits);
rfPos = NaN(size(brain));
% RF azimuth
rfPos(1:nUnits) = a1 + b1 .* brainX + c1 .* brainY;

% satisfy orthogonality constraint
if b1 == 0
    b2 = 1;
    c2 = 0;
elseif c1 == 0
    b2 = 0;
    c2 = 1;
else
    c2 = -b1 / c1;
    l = sqrt(1^2 + c2^2);
    b2 = 1/l;
    c2 = c2/l;
end
b2 = b2 * k2;
c2 = c2 * k2;
% RF elevation
rfPos(nUnits+1 : 2*nUnits) = a2 + b2 * brainX + c2 * brainY;