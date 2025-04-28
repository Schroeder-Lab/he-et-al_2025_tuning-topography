function Figure05(folders)

%% Parameters
[ds_trans, os_long, os_lat] = algebra.getDsOsAxes();

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure05');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Plot cartoons

% Cartoon of longitudes for direction preferences
% average body axis
body = mean([ds_trans(1,:); (ds_trans(2,:)+[180 0]).*[1 -1]], 1);
gravitation = mean([ds_trans(3,:); (ds_trans(4,:)+[180 0]).*[1 -1]], 1);

figure
hold on
[X,Y,Z] = sphere(100);
surf(X, Y, Z, 'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'none')
algebra.plotLongitudes(body, 12, 'b')
algebra.plotLongitudes(gravitation, 12, 'c')
algebra.plotLatitudes([0 90], 1, 'k')
algebra.plotLatitudes([0 0], 1, 'k')
view([-62, 25])
io.saveFigure(gcf, fPlots, 'cartoon_directionVectors');

figure
hold on
[X,Y,Z] = sphere(100);
surf(X, Y, Z, 'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'none')
algebra.plotLongitudes(os_long(1,:), 12, 'b')
algebra.plotLongitudes(os_long(2,:), 12, 'c')
algebra.plotLatitudes([0 90], 1, 'k')
algebra.plotLatitudes([0 0], 1, 'k')
view([-62, 25])
io.saveFigure(gcf, fPlots, 'cartoon_orientationLongitudeVectors');

figure
hold on
[X,Y,Z] = sphere(100);
surf(X, Y, Z, 'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'none')
algebra.plotLatitudes(os_lat(1,:), 8, 'r')
algebra.plotLatitudes(os_lat(2,:), 8, 'y')
algebra.plotLatitudes([0 90], 1, 'k')
algebra.plotLatitudes([0 0], 1, 'k')
view([-62, 25])
io.saveFigure(gcf, fPlots, 'cartoon_orientationLatitudeVectors');

%% Polar histograms of preferences per visual patch
Figure05_polarHists(data, fPlots, sets)