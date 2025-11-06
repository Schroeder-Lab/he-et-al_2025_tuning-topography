function Figure06S(folders, glob)

%% Parameters
maxP = 0.05; % p-value threshold

% to evaluate RFs
minEV = 0.01; % minimum explained variance
minPeak = 5; % minimum peak of RF (compared to noise)

%% Examples
exAnimal = 'FG010';
exDate = '2024-10-19';

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure06S');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Example LFP data to determine SC layer boundaries
Figure06S_LFP(folders, fPlots, exAnimal, exDate)

%% Load data: tuning preferences, RF position
% include all units within SC that are responsive to gratings
data = Figure06_loadData(folders, maxP, minEV, minPeak);

%% Plot direction vs orientation preference and selectivity
tuned = [cat(1, data.dirTuned), cat(1, data.oriTuned)];
dirPreferences = cat(1, data.dirPreferences);
oriPreferences = cat(1, data.oriPreferences);
dirSel = cat(1, data.dirSel);
oriSel = cat(1, data.oriSel);

% Direction vs orientation preference scatterplot
figure('Position', glob.figPositionDefault)
hold on
plot([0 180], [0 180], 'Color', [1 1 1].*0.5)
plot([180 360], [0 180], 'Color', [1 1 1].*0.5)
scatter(dirPreferences(all(tuned,2)), oriPreferences(all(tuned,2)), 15, ...
    'k', 'filled')
axis equal
xlim([-10 370])
ylim([-10 190])
set(gca, "Box", "off", "XTick", 0:90:360, "YTick", 0:90:180)
xlabel('Direction (deg)')
ylabel('Orientation (deg)')
title(sprintf('n = %d', sum(all(tuned,2))))
io.saveFigure(gcf, fPlots, 'tuning_preference_dirVsOriScatter');

% DS vs OS scatterplot
figure('Position', glob.figPositionDefault)
h = gscatter(dirSel(any(tuned,2)), oriSel(any(tuned,2)), ...
    tuned(any(tuned,2),1) + 2*tuned(any(tuned,2),2), ...
    repmat([0.4 0.8 0]', 1, 3), [], 15);
% hold on
% scatter(dirSel(indExamples), oriSel(indExamples), 40, ...
%     lines(length(indExamples)), "filled")
l = legend(h, 'DS', 'OS', 'DS & OS', "Location", "bestoutside");
l.Box = "off";
axis padded equal
mini = -0.05;
maxi = 1.05;
axis([mini maxi mini maxi])
set(gca, "Box", "off")
xlabel('Direction selectivity')
ylabel('Orientation selectivity')
io.saveFigure(gcf, fPlots, 'tuning_selectivity_dirVsOriScatter');

%% Compare preferred direction / orientation to expected value based on RF location
Figure06S_polarHistograms(glob, fPlots, data)