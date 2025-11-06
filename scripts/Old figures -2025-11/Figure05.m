function Figure05(folders, glob)

%% Parameters
sets = {'boutons', 'neurons'};
retinotopyRF = [false true]; % true: use RF positions estimated from 
                             % retinotopic mapping;
                             % false: use RF positions from RF mapping

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure05');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Plot cartoons

[ds_trans, os_long, os_lat] = algebra.getDsOsAxes();
% Cartoon of longitudes for direction preferences
% average body axis
body = mean( ...
    [mod(ds_trans(1,1), 360), ds_trans(1,2); ...
     mod(ds_trans(2,1)+180, 360), -ds_trans(2,2)], 1);
gravitation = mean( ...
    [mod(ds_trans(3,1), 360), ds_trans(3,2); ...
     mod(ds_trans(4,1)+180, 360), -ds_trans(4,2)], 1);

figure('Position', glob.figPositionDefault)
hold on
[X,Y,Z] = sphere(100);
surf(X, Y, Z, 'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'none')
algebra.plotLongitudes(body, 12, 'c')
algebra.plotLongitudes(gravitation, 12, 'b')
algebra.plotLatitudes([0 90], 1, 'k')
algebra.plotLatitudes([0 0], 1, 'k')
view([118, 25])
io.saveFigure(gcf, fPlots, 'cartoon_directionVectors');

figure('Position', glob.figPositionDefault)
hold on
[X,Y,Z] = sphere(100);
surf(X, Y, Z, 'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'none')
algebra.plotLongitudes(os_long(1,:), 12, 'c')
algebra.plotLongitudes(os_long(2,:), 12, 'b')
algebra.plotLatitudes([0 90], 1, 'k')
algebra.plotLatitudes([0 0], 1, 'k')
view([118, 25])
io.saveFigure(gcf, fPlots, 'cartoon_orientationLongitudeVectors');

figure('Position', glob.figPositionDefault)
hold on
[X,Y,Z] = sphere(100);
surf(X, Y, Z, 'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'none')
algebra.plotLatitudes(os_lat(1,:), 8, 'm')
algebra.plotLatitudes(os_lat(2,:), 8, 'r')
algebra.plotLatitudes([0 90], 1, 'k')
algebra.plotLatitudes([0 0], 1, 'k')
view([118, 25])
io.saveFigure(gcf, fPlots, 'cartoon_orientationLatitudeVectors');

%% Load data: RF position, tuning preferences
% data: .rfPos, .oriPref, .OSI, .dirPref, .DSI, .set
data = Figure04_loadData(folders, sets, retinotopyRF);

%% Polar histograms of preferences per visual patch
Figure05_polarHists(data, fPlots, sets)