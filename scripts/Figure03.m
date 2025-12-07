function Figure03(folders, glob)

%% Parameters
sets = {'boutons', 'neurons'};
retinotopyRF = [false true]; % true: use RF positions estimated from 
                             % retinotopic mapping;
                             % false: use RF positions from RF mapping
% for cartoons
cols = {'k', 'b', 'r', 'm'};

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure03');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Plot cartoons

[ds_trans, os_long, os_lat] = algebra.getDsOsAxes();

% 3D cartoons of topographies
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

% 2D cartoons of topographies
pos_azimuth = -120:30:0;
pos_elevation = 60:-30:-60;

figure('Position', glob.figPositionDefault)
tiledlayout(length(pos_elevation), length(pos_azimuth), "TileSpacing", "compact")
for y = pos_elevation
    for x = pos_azimuth
        nexttile
        for d = 1:4
            [~,theta] = algebra.getTranslationDir(ds_trans(d,:), [x y]);
            polarplot([0 theta], [0 1], cols{d}, 'LineWidth', 1)
            hold on
        end
        set(gca,'ThetaDir','clockwise', 'RTick', [], 'ThetaTickLabel', {})
    end
end
sgtitle(sprintf('Direction vectors [%d - %d azim, %d - %d elev]', ...
    pos_azimuth([1 end]), pos_elevation([1 end])))
io.saveFigure(gcf, fPlots, 'cartoon_directionVectors_2D')

figure('Position', glob.figPositionDefault)
tiledlayout(length(pos_elevation), length(pos_azimuth), "TileSpacing", "compact")
for y = pos_elevation
    for x = pos_azimuth
        nexttile
        for d = 1:2
            [~,theta] = algebra.getTranslationDir(os_long(d,:), [x y]);
            polarplot([theta theta+pi], [1 1], cols{d}, 'LineWidth', 1)
            hold on
        end
        set(gca,'ThetaDir','clockwise', 'RTick', [], 'ThetaTickLabel', {})
    end
end
sgtitle(sprintf('Orientation longitudes [%d - %d azim, %d - %d elev]', ...
    pos_azimuth([1 end]), pos_elevation([1 end])))
io.saveFigure(gcf, fPlots, 'cartoon_orientationLongitudeVectors_2D')

figure('Position', glob.figPositionDefault)
tiledlayout(length(pos_elevation), length(pos_azimuth), "TileSpacing", "compact")
for y = pos_elevation
    for x = pos_azimuth
        nexttile
        for d = 1:2
            [~,theta] = algebra.getLatitudeOrientation(os_lat(d,:), [x y]);
            polarplot([theta theta+pi], [1 1], cols{d+2}, 'LineWidth', 1)
            hold on
        end
        set(gca,'ThetaDir','clockwise', 'RTick', [], 'ThetaTickLabel', {})
    end
end
sgtitle(sprintf('Orientation latitudes [%d - %d azim, %d - %d elev]', ...
    pos_azimuth([1 end]), pos_elevation([1 end])))
io.saveFigure(gcf, fPlots, 'cartoon_orientationLatitudeVectors_2D')

%% Load data: RF position, tuning preferences
% data: .rfPos, .oriPref, .OSI, .dirPref, .DSI, .set
data = Figures_loadData(folders, sets, retinotopyRF);

%% Polar histograms of preferences per visual patch
Figure03_polarHists(glob, fPlots, data, sets)