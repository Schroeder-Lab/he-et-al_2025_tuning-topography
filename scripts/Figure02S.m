function Figure02S(folders)

%% Parameters
sets = {'boutons', 'neurons'};
% for evaluation of receptive fields (significance/goodness)
minEV = 0.01;
minPeak = 5;
% for evaluation of tuning selectivity
maxP = 0.05;
% to relate RF position to tuning preference
minUnits = 20;
smoothN = [201 29];
smoothPars = [0.1 0.5];
nPerm = 1000;
% plotting
labels = {'left', 'top'; 'direction', 'orientation'};
stdLims = [0.7 1.4];

%% For all plots
fPlots = fullfile(folders.plots, 'Figures', 'Figure02S');
if ~isfolder(fPlots)
    mkdir(fPlots)
end

%% Plot distance to monitor edge (left/top) vs pref. dir./ori.
for s = 1:2 % boutons and neurons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    rfDist = {}; % (1) from left edge, (2) from top edge
    prefs = {};
    tuned = {};

    % Collect data
    for subj = 1:length(subjDirs) % animals
        name = subjDirs(subj).name;
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            f = fullfile(folders.data, sets{s}, name, date);
            % ignore session if visual noise stimuolus was not present
            if ~isfile(fullfile(f, '_ss_rf.gaussFitPars.npy')) || ...
                    ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy'))
                continue
            end

            % load data
            data = io.getNoiseRFFits(f);
            fitPars = data.gaussPars;
            rfPos = fitPars(:, [2 4]); % (azimuth, elevation) in visual degrees
            ev_rf = data.EV;
            rf_peaks = data.peaks;
            edges = data.edges; % [left right top bottom]
            if isfile(fullfile(f, '_ss_rf.outliers.npy'))
                outliers = readNPY(fullfile(f, '_ss_rf.outliers.npy'));
            else
                outliers = false(size(ev_rf));
            end
            [dirTuning, oriTuning] = ...
                io.getTuningResults(f, 'gratingsDrifting');

            % only consider significant RFs
            valid = ev_rf >= minEV & rf_peaks >= minPeak;
            valid = valid & ~outliers;
            if sum(valid & dirTuning.pValue < maxP) < minUnits && ...
                    sum(valid & oriTuning.pValue < maxP) < minUnits
                continue
            end

            % collect data
            rfDist{end+1,1} = rfPos(valid,:).*[1 -1] + [-edges(1) edges(3)];
            prefs{end+1,1} = dirTuning.preference(valid);
            tuned{end+1,1} = dirTuning.pValue(valid);
            prefs{end,2} = oriTuning.preference(valid);
            tuned{end,2} = oriTuning.pValue(valid);
        end
    end
    dst = cat(1, rfDist{:});

    % Plot RF distance to monitor edge vs tuning preference for each unit
    pr = [cat(1, prefs{:,1}), cat(1, prefs{:,2})];
    tn = [cat(1, tuned{:,1}), cat(1, tuned{:,2})] < maxP;
    figure('Position', [40 580 1850 410])
    tiledlayout(1, 4)
    for feat = 1:2
        for ed = 1:2 % left / top edge
            nexttile
            hold on
            scatter(dst(tn(:,feat),ed), pr(tn(:,feat),feat), 15, 'k', ...
                "filled")
            xlim([0 max(dst(:,ed))])
            ylim([-5 365]./feat)
            set(gca, "Box", "off", "YTick", (0:90:360)./feat)
            xlabel(sprintf('Distance to %s monitor edge (deg)', labels{1,ed}))
            ylabel(sprintf('Preferred %s (deg)', labels{2,feat}))
            title(sprintf('n = %d', sum(tn(:,feat))))
        end
    end
    sgtitle(sets{s})
    io.saveFigure(gcf, fPlots, ...
        sprintf('rf-to-monitorEdge_%s_scatter', sets{s}))

    % Plot circular STD as function of RF distance to monitor edge
    pr = [cat(1, prefs{:,1}), cat(1, prefs{:,2})];
    tn = [cat(1, tuned{:,1}), cat(1, tuned{:,2})] < maxP;
    figure('Position', [40 580 1850 410])
    tiledlayout(1, 4)
    for feat = 1:2
        for ed = 1:2 % left / top edge
            nexttile
            hold on
            [dst_sorted, order] = sort(dst(tn(:,feat),ed), "ascend");
            pr_sorted = pr(tn(:,feat),feat);
            pr_sorted = pr_sorted(order);
            dst_sm = smoothdata(dst_sorted, "movmean", smoothN(1));
            dst_sm = dst_sm(ceil(smoothN(1)/2) : end-floor(smoothN(1)/2));
            std_sm  = NaN(length(pr_sorted) - smoothN(1) + 1, 1);
            for k = 1:length(std_sm)
                std_sm(k) = circ_std(deg2rad( ...
                    pr_sorted(k:k+smoothN(1)-1) .* feat));
            end
            f1 = fit(dst_sm, std_sm, "smoothingspline", ...
                "SmoothingParam", smoothPars(1));

            % permutation test
            std_null = NaN(size(std_sm,1), nPerm);
            rng('default');
            for n = 1:nPerm
                order = randperm(length(pr_sorted));
                pr_perm = pr_sorted(order);
                for k = 1:size(std_null,1)
                    std_null(k,n) = circ_std(deg2rad( ...
                        pr_perm(k:k+smoothN(1)-1) .* feat));
                end
            end
            std_prctl = prctile(std_null, [2.5 50 97.5], 2);
            f2 = cell(1,3);
            for k=1:3
                f2{k} = fit(dst_sm, std_prctl(:,k), "smoothingspline", ...
                    "SmoothingParam", smoothPars(1));
            end

            x = floor(dst_sm(1)*10)/10 : 0.1 : ceil(dst_sm(end)*10)/10;
            fill([x flip(x)], [f2{1}(x); flip(f2{3}(x))], 'k', ...
                'FaceAlpha', 0.3, 'EdgeColor', 'none')
            plot(x, f2{2}(x), 'k', 'LineWidth', 2)
            plot(x, f1(x), 'r', 'LineWidth', 2)

            xlim([0 max(dst(:,ed))])
            ylim(stdLims)
            set(gca, "Box", "off")
            xlabel(sprintf('Distance to %s monitor edge (deg)', labels{1,ed}))
            ylabel(sprintf('STD across %ss', labels{2,feat}))
        end
    end
    sgtitle(sets{s})
    io.saveFigure(gcf, fPlots, ...
        sprintf('rf-to-monitorEdge_%s_STD', sets{s}))
end