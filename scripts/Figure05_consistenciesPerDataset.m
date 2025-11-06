function Figure05_consistenciesPerDataset(glob, fPlots, data, sets, ...
    retinotopyRF, measures)

%% Parameters
% create surrogate global maps of direction/orientation tuning
numPerm = 1000;
% plotting
gridDist = 2;
gridRadius = 5;

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
                spatial.makeSmoothMap(rfs, prefs .* m, gridDist, gridRadius);
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
                    rfs(valid,:), prefs(valid(order)) .* m, ...
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
        valid = ~cellfun(@isempty, consistency(s).(measures{m}).original);
        lineColors = turbo(length(valid));

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
        figure('Position', glob.figPositionDefault)
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