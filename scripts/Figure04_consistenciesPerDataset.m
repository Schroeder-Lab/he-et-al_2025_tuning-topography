function Figure04_consistenciesPerDataset(data, fPlots, sets, ...
    retinotopyRF, measures)

%% Parameters
% create surrogate global maps of direction/orientation tuning
numPerm = 1000;
% plotting
gridDist = 2;
gridRadius = 5;
binEdges = 0:0.1:1;
bins = binEdges(2:end) - 0.05;
binsSmooth = linspace(bins(1), bins(end), 100);

%% Determine consistencies and null distributions for each dataset
consistency = struct('direction', cell(2,1), 'orientation', []);
strPrefs = {'dirPref', 'oriPref'};
for s = 1:2
    numRecs = max(data(s).set);
    consistency(s).direction.original = cell(numRecs,1);
    consistency(s).direction.null = cell(numRecs,1);
    consistency(s).orientation.original = cell(numRecs,1);
    consistency(s).orientation.null = cell(numRecs,1);
    for d = 1:numRecs
        indRec = data(s).set == d;
        for m = 1:2
            prefs = data(s).(strPrefs{m})(indRec);
            rfs = data(s).rfPos(indRec,:);
            [x,y,~,consistencies] = ...
                spatial.makeSmoothMap(rfs, prefs, gridDist, gridRadius);
            % only consider datasets where RF positions span a distance of
            % at least gridRadius (otherwhise the permuted null data will
            % give the same results as the original data)
            if isempty(x) || diff(x([1 end]))<gridRadius || ...
                    diff(y([1 end]))<gridRadius
                continue
            end
            consistency(s).(measures{m}).original{d} = consistencies;
            consistency(s).(measures{m}).null{d} = ...
                NaN([size(consistencies) numPerm]);

            valid = find(~any(isnan([rfs, prefs]), 2));
            rng('default');
            for p = 1:numPerm
                order = randperm(length(valid));
                [~, ~,~, cons] = spatial.makeSmoothMap(...
                    rfs(valid,:), prefs(valid(order)), ...
                    gridDist, gridRadius);
                consistency(s).(measures{m}).null{d}(:,:,p) = cons;
            end
        end
    end
end

%% Plot consistencies compared to null distribution, per dataset
for s = 1:2
    if retinotopyRF(s)
        str = 'retinoptopic';
    else
        str = 'measured';
    end
    for m = 1:2
        probs = NaN(length(bins), ...
            length(consistency(s).(measures{m}).original));
        probsNull = NaN(length(bins), ...
            length(consistency(s).(measures{m}).original));
        probsNullDistr = NaN(length(bins), ...
            length(consistency(s).(measures{m}).original), numPerm);
        for d = 1:length(consistency(s).(measures{m}).original)
            if isempty(consistency(s).(measures{m}).original{d})
                continue
            end
            probs(:,d) = histcounts( ...
                consistency(s).(measures{m}).original{d}, binEdges, ...
                "Normalization", "probability");
            probsNull(:,d) = histcounts( ...
                consistency(s).(measures{m}).null{d}, binEdges, ...
                "Normalization", "probability");
            for p = 1:numPerm
                probsNullDistr(:,d,p) = histcounts( ...
                consistency(s).(measures{m}).null{d}(:,:,p), binEdges, ...
                "Normalization", "probability");
            end
        end
        valid = ~any(isnan(probs),1);
        lineColors = turbo(length(valid));

        % consistencies of null data
        smNull = NaN(length(binsSmooth), length(valid));
        smNull(:,valid) = interp1(bins, probsNull(:,valid), ...
            binsSmooth, "pchip");
        smNull = smNull ./ sum(smNull,1);
        % consistencies of original data
        smOrig = NaN(length(binsSmooth), length(valid));
        smOrig(:,valid) = interp1(bins, probs(:,valid), ...
            binsSmooth, "pchip");
        smOrig = smOrig ./ sum(smOrig,1);
        figure
        hold on
        for v = 1:length(valid)
            if ~valid(v)
                continue
            end
            % consistencies of null data
            plot(binsSmooth, smNull(:,v), ':', "LineWidth", 1, ...
                "Color", lineColors(v,:));
            % consistencies of original data
            plot(binsSmooth, smOrig(:,v), "LineWidth", 1, ...
                "Color", lineColors(v,:));
        end
        h = [0 0];
        h(1) = plot(binsSmooth, mean(smNull,2,"omitnan"), 'k:', ...
            "LineWidth", 2);
        h(2) = plot(binsSmooth, mean(smOrig,2,"omitnan"), 'k', ...
            "LineWidth", 2);
        set(gca, "Box", "off")
        legend(h, 'null', 'original')
        xlabel(sprintf('Consistency of %s preferences', measures{m}))
        ylabel('Probability')
        title(sprintf('%s consistency (%sRFs) - %s', ...
            measures{m}, str, sets{s}))
        io.saveFigure(gcf, fPlots, sprintf(...
            'consistency_%s_%s_histograms-per-dataset_%sRFs', ...
            sets{s}, measures{m}, str))

        % mean consistencies compared to null distribution of mean
        mn = cellfun(@mean, consistency(s).(measures{m}).original, ...
            repmat({"all"},length(valid),1), ...
            repmat({"omitnan"},length(valid),1));
        nullM = cellfun(@mean, consistency(s).(measures{m}).null, ...
            repmat({[1 2]},length(valid),1), ...
            repmat({"omitnan"},length(valid),1), 'UniformOutput', false);
        nullM = cellfun(@squeeze, nullM, 'UniformOutput', false);
        nullM_flat = NaN(numPerm, length(valid));
        nullM_flat(:,valid) = cat(2, nullM{valid});
        figure
        hold on
        h = [0 0 0];
        for v = 1:length(valid)
            if ~valid(v)
                continue
            end
            inv = prctile(nullM_flat(:,v), [2.5 97.5]);
            h(1) = plot(inv, [1 1].*v, ...
                "Color", lineColors(v,:));
            h(2) = plot(prctile(nullM_flat(:,v), 50), v, 'o', ...
                "MarkerFaceColor", lineColors(v,:), "MarkerEdgeColor", "none");
            h(3) = plot(mn(v), v, 'v', ...
                "MarkerFaceColor", lineColors(v,:), "MarkerEdgeColor", "none");
            if mn(v)>inv(2) || mn(v)<inv(1)
                plot(1.1, v, 'k*')
            end
        end
        set(gca, "Box", "off")
        legend(h, 'null interval', 'null median', 'original', "Location", "bestoutside")
        xlim([0 1.15])
        ylim([find(valid,1)-1 find(valid,1,'last')+1])
        xlabel(sprintf('Mean consistency of %s preferences', measures{m}))
        ylabel('Dataset')
        title('Original relative to null')
        io.saveFigure(gcf, fPlots, sprintf(...
            'consistency_%s_%s_intervals-per-dataset_%sRFs', ...
            sets{s}, measures{m}, str))
    end
end