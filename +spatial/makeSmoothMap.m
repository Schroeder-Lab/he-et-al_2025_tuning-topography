function [x, y, preferences, consistencies, counts] = ...
    makeSmoothMap(rfPos, directions, binSize, radius)

valid = ~any(isnan([rfPos, directions]), 2);

x = (floor(min(rfPos(valid,1))/binSize) : ceil(max(rfPos(valid,1))/binSize)) .* binSize;
y = (floor(min(rfPos(valid,2))/binSize) : ceil(max(rfPos(valid,2))/binSize)) .* binSize;
[x,y] = meshgrid(x,y);

% transform preferred directions into vectors
[dirX, dirY] = pol2cart(deg2rad(directions(valid)), ones(sum(valid),1));

% calculate mean direction, variance, and counts for each grid point
preferences = NaN(size(x));
consistencies = NaN(size(x));
counts = NaN(size(x));
for k = 1:numel(x)
    % distance of grid point to RF of each neuron
    d = sqrt((x(k) - rfPos(valid,1)).^2 + (y(k) - rfPos(valid,2)).^2);
    % find units within radius
    within = d <= radius;
    % ignore gridpoint with fewer than 5 within units
    if sum(within) < 5
        continue
    end
    % polar coordinates of vector average of direction preferences
    [theta, rho] = cart2pol(mean(dirX(within)), mean(dirY(within)));
    preferences(k) = theta;
    consistencies(k) = rho;
    counts(k) = sum(within);
end

preferences = mod(rad2deg(preferences), 360);