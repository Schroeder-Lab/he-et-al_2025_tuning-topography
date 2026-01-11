function plotRFOutlines(rfGaussPars, EVs, peaks, minEV, minPeak, ...
    exUnits, edges_rf, exColors)
%PLOTRFOUTLINES   Plot contours of all fitted 2D Gaussians.

% INPUTS
% rfGaussPars   [units x 7], parameters of fitted 2D Gaussians
% EVs           [units], explained variances of fitted spatio-temporal RFs
% peaks         [units], peak heights of RFs
% minEV         double, minimum EV for significant RFs
% minPeak       double, minimum peak for significant RFs
% exUnits       [k], indices of example units
% edges_rf      [left rigth top bottom] of RF map to be plotted
% exColors      [k, 3], RGBs for example units

ellipse_x = linspace(-pi, pi, 100);

figure
hold on
h = zeros(size(rfGaussPars,1),1);
for iUnit = 1:size(rfGaussPars,1)
    if EVs(iUnit) < minEV || peaks(iUnit) < minPeak
        continue
    end
    % ellipse at 2 STD (x and y), not rotated, not shifted
    x = rfGaussPars(iUnit,3) * cos(ellipse_x);
    y = rfGaussPars(iUnit,5) * sin(ellipse_x);
    % rotate and shift ellipse
    x_rot = rfGaussPars(iUnit,2) + ...
        x .* cos(rfGaussPars(iUnit,6)) - ...
        y .* sin(rfGaussPars(iUnit,6));
    y_rot = rfGaussPars(iUnit,4) + ...
        x .* sin(rfGaussPars(iUnit,6)) + ...
        y .* cos(rfGaussPars(iUnit,6));
    h(iUnit) = plot(x_rot, y_rot, 'k');
end
for k = 1:length(exUnits)
    iUnit = exUnits(k);
    set(h(iUnit),'Color',exColors(k,:))
    uistack(h(iUnit),'top')
    set(h(iUnit),'LineWidth',2)
end
legend(h(exUnits))
axis image
axis(edges_rf([1 2 4 3]))
xlabel('Azimuth (visual degrees)')
ylabel('Elevation (visual degrees)')