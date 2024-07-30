function [parameters, gaussMap] = fit2dGaussRF(profile, doPlot)
%FIT2DGAUSSRF   Fit 2D Gaussian to 2D profile.

% INPUTS
% profile       [rows x columns], 2D profile to which Gauss is fitted
% doPlot        logical, make plot if true

% OUTPUTS
% parameters    [amplitude, x-center, width, y-center, height, rotation],
%               parameters of Gaussian fit
% gaussMap      [rows x columns], map of 2D Gaussian

% create mesh with x- and y-coordinates of profile
xcoords = 1:size(profile,2);
ycoords = 1:size(profile,1);
[x, y] = meshgrid(xcoords, ycoords);
xdata = cat(3, x, y);

% x0 is the initial guess
[maxY, maxX] = find(profile == max(profile(:)),1);
x0 = [1, maxX, mean(diff(xcoords))*5, maxY, mean(diff(ycoords))*5, 0];
% lower and upper bounds of parameters
lb = [0, min(xcoords), 0, min(ycoords), 0, -pi/4];
ub = [realmax('double'), max(xcoords), max(xcoords)^2, max(ycoords), ...
    max(ycoords)^2, pi/4];
options = optimoptions('lsqcurvefit', 'Display', 'off');
% run fitting procedure to retrieve parameters
parameters = lsqcurvefit(@rf.D2GaussFunctionRot, x0, xdata, ...
    profile, lb, ub, options);
% create 2D map of Gaussian using parameters
gaussMap = rf.D2GaussFunctionRot(parameters, xdata);

if doPlot
    figure;
    subplot(3,1,1);
    imagesc(xcoords, ycoords, profile);
    colormap hot;

    subplot(3,1,2);
    imagesc(xcoords, ycoords, gaussMap);

    subplot(3,1,3);
    contour(gaussMap,1,'Color', 'k', 'LineWidth', 2.0);
    set(gca, 'YDir', 'reverse');
end