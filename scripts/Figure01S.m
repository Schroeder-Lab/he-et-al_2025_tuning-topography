function Figure01S(folders)

%% Parameters
maxP = 0.05;
minR2 = 0.02;
exp = {'gratingsDrifting', 'bars', 'gratingsStatic'};

%% Examples
% datasets with drifting gratings, static gratings, and bars
ex = {'SS044', '2015-05-15', [205 382 5 401 431 126 422]};
    %'SS047', '2015-12-03'; 'SS048', '2015-11-09'

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure01S');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Example calcium traces and tuning curves
f = fullfile(folders.data, 'neurons', ex{1,1}, ex{1,2});
krnlFits = cell(1,3);
P_ex = [];
R2_ex = [];
for k = 1:3
    krnlFits{k} = io.getStimResponseFits(f, exp{k});
    P_ex = [P_ex, krnlFits{k}.pValue];
    R2_ex = [R2_ex, krnlFits{k}.R2];
end
R2_ex(P_ex > maxP) = -100;
[R2_sorted, order] = sort(mean(R2_ex,2), "descend");

%% Preferences across stimulus paradigms
% Collect data
subjDirs = dir(fullfile(folders.data, 'neurons', 'SS*'));
R2 = [];
dirPreferences = [];
oriPreferences = [];
dirTuned = [];
oriTuned = [];
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folders.data, 'neurons', name, '2*'));
    for d = 1:length(dateDirs) %dates
        date = dateDirs(d).name;
        f = fullfile(folders.data, 'neurons', name, date);
        r2 = [];
        for k = 1:3
            % ignore session if stimulus was not presented
            if ~isfile(fullfile(f, sprintf('_ss_%s.intervals.npy', exp{k})))
                continue
            end

            krnlFits = io.getStimResponseFits(f, exp{k});
            [dirTuning, oriTuning] = io.getTuningResults(f, exp{k});
            numUnits = length(krnlFits.pValue);
            if isempty(r2)
                r2 = NaN(numUnits,3);
                dp = NaN(numUnits,2);
                op = NaN(numUnits,3);
                dt = NaN(numUnits,2);
                ot = NaN(numUnits,3);
            end

            r2(:,k) = krnlFits.R2;
            if ~isempty(dirTuning)
                dp(:,k) = dirTuning.preference;
                dt(:,k) = dirTuning.pValue < maxP;
            end
            op(:,k) = oriTuning.preference;
            ot(:,k) = oriTuning.pValue < maxP;
        end
        R2 = [R2; r2];
        dirPreferences = [dirPreferences; dp];
        oriPreferences = [oriPreferences; op];
        dirTuned = [dirTuned; dt];
        oriTuned = [oriTuned; ot];

        % if strcmp(name, ex{s,1}) && strcmp(date, ex{s,2})
        %     indExamples = NaN(length(ex{s,4}),1);
        %     unitsResponsive = find(unitsResponsive);
        %     n = length(dataset) - length(unitsResponsive);
        %     for k = 1:length(indExamples)
        %         indExamples(k) = n + find(ex{s,4}(k) == unitsResponsive);
        %     end
        % end
    end
end

% Scatterplot: preferred direction from drifting gratings vs bars
ind = all(R2(:,[1 2]) > minR2, 2) & all(dirTuned(:, [1 2]), 2);
diffs = dirPreferences(ind,1) - dirPreferences(ind,2);
diffs(diffs < -180) = 360 + diffs(diffs < -180);
diffs(diffs > 180) = diffs(diffs > 180) - 360;
m = mean(diffs);
p = signrank(diffs);
% [~,p] = ttest(diffs);
figure
hold on
plot([0 360], [0 360], 'Color', [1 1 1].*0.5)
scatter(dirPreferences(ind,1), dirPreferences(ind,2), 15, 'k', 'filled')
axis equal
xlim([-10 370])
ylim([-10 370])
set(gca, "Box", "off", "XTick", 0:90:360, "YTick", 0:90:360)
xlabel(exp{1})
ylabel(exp{2})
title(sprintf('Preferred directions (n=%d, <x-y>=%.1f, p=%.3f)', ...
    sum(ind), m, p))

% Scatterplot: preferred orientation from drifting gratings vs bars
ind = all(R2(:,[1 2]) > minR2, 2) & all(oriTuned(:, [1 2]), 2);
diffs = oriPreferences(ind,1) - oriPreferences(ind,2);
diffs(diffs < -90) = 180 + diffs(diffs < -90);
diffs(diffs > 90) = diffs(diffs > 90) - 180;
m = mean(diffs);
p = signrank(diffs);
% [~,p] = ttest(diffs);
figure
hold on
plot([0 360], [0 360], 'Color', [1 1 1].*0.5)
scatter(oriPreferences(ind,1), oriPreferences(ind,2), 15, 'k', 'filled')
axis equal
xlim([-10 190])
ylim([-10 190])
set(gca, "Box", "off", "XTick", 0:45:180, "YTick", 0:45:180)
xlabel(exp{1})
ylabel(exp{2})
title(sprintf('Preferred orientations (n=%d, <x-y>=%.1f, p=%.3f)', ...
    sum(ind), m, p))

% Scatterplot: preferred orientation from drifting vs static gratings
ind = all(R2(:,[1 3]) > minR2, 2) & all(oriTuned(:, [1 3]), 2);
diffs = oriPreferences(ind,1) - oriPreferences(ind,3);
diffs(diffs < -90) = 180 + diffs(diffs < -90);
diffs(diffs > 90) = diffs(diffs > 90) - 180;
m = mean(diffs);
p = signrank(diffs);
% [~,p] = ttest(diffs);
figure
hold on
plot([0 360], [0 360], 'Color', [1 1 1].*0.5)
scatter(oriPreferences(ind,1), oriPreferences(ind,3), 15, 'k', 'filled')
axis equal
xlim([-10 190])
ylim([-10 190])
set(gca, "Box", "off", "XTick", 0:45:180, "YTick", 0:45:180)
xlabel(exp{1})
ylabel(exp{3})
title(sprintf('Preferred orientations (n=%d), <x-y>=%.1f, p=%.3f', ...
    sum(ind), m, p))

% Scatterplot: preferred orientation from static gratings vs bars
ind = all(R2(:,[3 2]) > minR2, 2) & all(oriTuned(:, [3 2]), 2);
diffs = oriPreferences(ind,3) - oriPreferences(ind,2);
diffs(diffs < -90) = 180 + diffs(diffs < -90);
diffs(diffs > 90) = diffs(diffs > 90) - 180;
m = mean(diffs);
p = signrank(diffs);
% [~,p] = ttest(diffs);
figure
hold on
plot([0 360], [0 360], 'Color', [1 1 1].*0.5)
scatter(oriPreferences(ind,3), oriPreferences(ind,2), 15, 'k', 'filled')
axis equal
xlim([-10 190])
ylim([-10 190])
set(gca, "Box", "off", "XTick", 0:45:180, "YTick", 0:45:180)
xlabel(exp{3})
ylabel(exp{2})
title(sprintf('Preferred orientations (n=%d), <x-y>=%.1f, p=%.3f', ...
    sum(ind), m, p))