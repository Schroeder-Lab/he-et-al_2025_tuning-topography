function Figure06_exampleRFs(folders, fPlots, glob, tuningData, ...
    minEV, minPeak, exUnits, exColors)

edges_rf = [-94 -47 41.25 -22.125];
RFtypes = {'ON', 'OFF', 'ON+OFF'};

f = fullfile(folders.data, 'ephys', tuningData.animal, ...
    tuningData.date);
spikeData = io.getEphysData(f);
results = io.getNoiseRFFits(f);
edges = results.edges;
gridW = diff(edges(1:2)) / size(results.maps,3);
gridH = -diff(edges(3:4)) / size(results.maps,2);

% example RF maps of single units
for iUnit = 1:length(exUnits)
    unit = find(spikeData.clusterIDs == exUnits(iUnit));
    [~, mxTime] = max(results.timeWeights(unit,:));
    % rf: [rows x columns x subfield]
    rfield = squeeze(results.maps(unit,:,:,mxTime,:));
    rf.plotRF(rfield, results.gaussPars(unit,:), ...
        results.bestSubFields(unit), edges, edges_rf, gridW, gridH)
    sgtitle(sprintf('Unit %d (EV: %.3f, peak/noise: %.1f, %s)', ...
        exUnits(iUnit), results.EV(unit), results.peaks(unit), ...
        RFtypes{results.bestSubFields(unit)}))

    io.saveFigure(gcf, fPlots, sprintf('example_RF_%s_%s_%03d', ...
        tuningData.animal, ...
        tuningData.date, exUnits(iUnit)))
end

% RF outlines of all units in example session
indAll = spikeData.clusterDepths(:,2) > 0 & results.EV >= minEV & ...
    results.peaks >= minPeak;
[~, indEx] = ismember(exUnits, spikeData.clusterIDs(indAll));
rf.plotRFOutlines(results.gaussPars(indAll,:), results.EV(indAll), ...
    results.peaks(indAll), minEV, minPeak, indEx, edges_rf, exColors)
title(sprintf('%s %s (n = %d)', tuningData.animal, ...
    tuningData.date, sum(indAll)))
set(gcf, 'Position', glob.figPositionDefault)
io.saveFigure(gcf, fPlots, sprintf('example_RFoutlines_%s_%s', ...
    tuningData.animal, tuningData.date))

% depth of units within SC
cols = zeros(sum(indAll), 3);
cols(indEx,:) = exColors;

figure('Position', glob.figPositionDefault)
scatter(zeros(sum(indAll), 1), spikeData.clusterDepths(indAll,1) - ...
    tuningData.SO_depth, 30, cols, 'filled', '<');
set(gca, "YDir", "reverse")
ylabel('Depth from SGS-SO border (in um)')
title(sprintf('%s %s (n = %d)', tuningData.animal, ...
    tuningData.date, sum(indAll)))
io.saveFigure(gcf, fPlots, sprintf('example_depths_%s_%s', ...
    tuningData.animal, tuningData.date))