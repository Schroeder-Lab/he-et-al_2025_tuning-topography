function main_calculateDS_OS_ephys(folders)
% Given the spike responses to gratings,
% use vector averaging to determine preferred direction and orientation 
% (direction of final vector) as well as direction and orientation 
% selectivity (length of final vector). Use permutation test to determine
% significance of direction and orientation selectivity.

%% Parameters
% stimulus response
stimulusDur = 2;
durationSlack = 0.01;
baselineTime = -0.5;

% testing tuning selectivity
numShuffles = 1000; % number of permutations for testing significance of 
                    % direction and orientation selectivity
maxP = 0.05; % p-value threshold for direction/orientation selectivity

% for plotting
cols = [1 0 0; 1 0.8 0.8];
layers = {'sSC', 'dSC'};

%% Plot direction and orientation tuning curves
subjDirs = dir(fullfile(folders.data, 'ephys'));
subjDirs = subjDirs(~startsWith({subjDirs.name}, '.') & [subjDirs.isdir]);
for subj = 1:length(subjDirs) % animals
    name = subjDirs(subj).name;
    dateDirs = dir(fullfile(folders.data, 'ephys', name, '2*'));
    for dt = 1:length(dateDirs) %dates
        date = dateDirs(dt).name;
        f = fullfile(folders.data, 'ephys', name, date);
            
        % ignore session if stimulus was not presented
        if ~isfile(fullfile(f, '_ss_gratingsDrifting.intervals.npy'))
            continue
        end
        fPlots = fullfile(folders.plots, '02_TuningCurves', ...
            'gratingsDrifting', 'ephys', name, date);
        if ~isfolder(fPlots)
            mkdir(fPlots)
        end

        % load data
        spikeData = io.getEphysData(f);
        stimData = io.getGratingInfo(f, 'gratingsDrifting');

        % determine tuning related values
        stimDirs = stimData.directions;
        stimDist = median(diff(stimDirs));
        stimReps = max(histcounts(stimData.ids, 0.5:length(stimDirs)+1));
        stimDurs = diff(stimData.times,1,2);
        stimDirsCirc = [stimDirs; 360];
        stimOrisCirc = stimDirs(1:length(stimDirs)/2+1);

        % ignore spike data before/after visual paradigm
        t_ind = spikeData.times >= stimData.times(1) - 2 & ...
            spikeData.times <= stimData.times(end) + 2;
        spikeData.times = spikeData.times(t_ind);
        spikeData.clusters = spikeData.clusters(t_ind);

        % get baseline-subtracted firing rates for each trial
        units = spikeData.clusterIDs;
        amplitudes = NaN(stimReps, length(stimDirs), length(units));
        p_compareToBaseline = NaN(length(units), 1);
        p_compareAcrossStimuli = NaN(length(units), 1);
        % initialize results structures
        dirTuning = struct('preference', ...
            repmat({NaN(1,1)}, length(units), 1), ...
            'selectivity', NaN, 'pValue', NaN, 'responseSign', NaN);
        oriTuning = struct('preference', ...
            repmat({NaN(1,1)}, length(units), 1), ...
            'selectivity', NaN, 'pValue', NaN, 'responseSign', NaN);
        % identify trials that are too short
        invalidTrials = find(stimDurs < stimulusDur - durationSlack);
        % identify units within the SC
        unitsSC = find(spikeData.clusterDepths(:,2) > 0);

        for iUnit = 1:length(unitsSC)
            unit = unitsSC(iUnit);
            % determine baseline-subtracted visual responses for each trial
            t = spikeData.times(spikeData.clusters == units(unit));
            [sp_aligned, trial] = events.alignData(t, ...
                stimData.times(:,1), [baselineTime stimulusDur]);
            if isempty(sp_aligned)
                continue
            end
            for stim = 1:length(stimDirs)
                stimTrials = find(stimData.ids == stim);
                for rep = 1:length(stimTrials)
                    tr = stimTrials(rep);
                    if ismember(tr, invalidTrials)
                        continue
                    end
                    ind_spikes = trial == tr;
                    amplitudes(rep, stim, unit) = ...
                        sum(sp_aligned(ind_spikes) >= 0) / ...
                        min(stimDurs(tr), stimulusDur) - ...
                        sum(sp_aligned(ind_spikes) < 0) / (-baselineTime);
                end
            end

            % Test whether unit is responsive
            amps = reshape(amplitudes(:,:,unit), [], 1);
            groups = reshape(repmat(1:length(stimDirs), stimReps, 1), ...
                [], 1);
            invalid = isnan(amps);
            amps(invalid) = [];
            groups(invalid) = [];
            T = table(amps, groups, ...
                'VariableNames', {'amps', 'stim'});
            % Use a linear model (fixed effects) (not a mixed-effects model
            % because we do not treat stimuli as random effects). We make
            % the intercept equal to the grand mean (effects coding), so we
            % can test whether the "overall mean" is different from zero,
            % i.e. visual responses are different from baseline.
            mdl = fitlm(T, 'amps ~ 1 + stim', ...
                'CategoricalVars', 'stim', 'DummyVarCoding', 'effects');
            % Test whether "grand mean" is zero, which would mean that 
            % overall the visual responses are not different from baseline.
            coeffs = mdl.Coefficients;
            p_compareToBaseline(unit) = coeffs{'(Intercept)', 'pValue'};
            % Now test whether responses to stimuli are different from each
            % other (same as ANOVA).
            p_compareAcrossStimuli(unit) = coefTest(mdl);
            
            % disregard non-responsive units from further analysis
            if (p_compareToBaseline(unit) > maxP && ...
                    p_compareAcrossStimuli(unit) > maxP) 
                continue
            end

            % determine direction/orientation preference,
            % selectivity, and signficance (response sign: -1 if
            % most responses are negative)
            [drct, ori] = tuning.determineTuning(amplitudes(:,:,unit), ...
                stimDirs, numShuffles);
            dirTuning(unit) = drct;
            oriTuning(unit) = ori;

            % plot results
            % replicate responses for 0 degrees at 360 or 180 degrees
            ampsDir = amplitudes(:, [1:end,1], unit) .* drct.responseSign;
            ampsOri = [amplitudes(:,1:length(stimDirs)/2, unit); ...
                amplitudes(:,length(stimDirs)/2+1:end, unit)] .* ...
                ori.responseSign;
            ampsOri = ampsOri(:, [1:end 1]);

            figure('Position', [100 380 1180 420]);
            tiledlayout(1,2)
            % direction tuning curve
            maxi = max(ampsDir,[],'all');
            mini = min(ampsDir,[],'all');
            nexttile
            c = cols(1,:);
            if drct.pValue > maxP
                c = cols(2,:);
            end
            plot(stimDirsCirc' + randn(size(ampsDir)) .* (stimDist/10), ...
                ampsDir, '.', 'MarkerSize', 10, ...
                'Color', [1 1 1].*0.5);
            hold on
            plot(stimDirsCirc, mean(ampsDir,1, 'omitnan'), ...
                'k.-', 'MarkerSize', 30, 'LineWidth', 1);
            plot(drct.preference, maxi, 'v', 'MarkerSize', 8, ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', 'none');
            xlim([-10 370])
            ylim([min([0 1.1*mini]) 1.1*maxi])
            xlabel('Direction (deg)')
            title(sprintf('DS = %.3f, p = %.3f', drct.selectivity, ...
                drct.pValue))
            set(gca, 'box', 'off', 'XTick', stimDirsCirc(1:3:end))

            % orientation tuning curve
            maxi = max(ampsOri,[],'all');
            mini = min(ampsOri,[],'all');
            nexttile 
            c = cols(1,:);
            h = [0 0 0];
            if ori.pValue > maxP
                c = cols(2,:);
            end
            hs = plot(stimOrisCirc' + randn(size(ampsOri)) .* (stimDist/10), ...
                ampsOri, '.', 'MarkerSize', 10, ...
                'Color', [1 1 1].*0.5);
            h(1) = hs(1);
            hold on
            h(2) = plot(stimOrisCirc, mean(ampsOri,1, 'omitnan'), ...
                'k.-', 'MarkerSize', 30, 'LineWidth', 1);
            h(3) = plot(ori.preference, maxi, 'v', 'MarkerSize', 8, ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', 'none');
            legend(h, 'single trials','median','preference', ...
                'Location', 'bestoutside')
            set(gca, 'box', 'off', 'XTick', stimOrisCirc(1:2:end))
            xlim([-10 190])
            ylim([min([0 1.1*mini]) 1.1*maxi])
            xlabel('Orientation (deg)')
            title(sprintf('OS = %.3f, p = %.3f', ori.selectivity, ...
                ori.pValue))
            sgtitle(sprintf('Cluster %d in %s (%d um from surface)', ...
                units(unit), layers{spikeData.clusterDepths(unit,2)}, ...
                round(spikeData.clusterDepths(unit,1))))

            % save and close plot
            saveas(gcf, fullfile(fPlots, sprintf('Unit%03d.jpg', unit)))
            close(gcf)
        end
        % save data
        writeNPY(amplitudes, ...
            fullfile(f, '_ss_gratingsDriftingTrials.amplitudes.npy'))
        writeNPY(p_compareToBaseline, ...
            fullfile(f, '_ss_gratingsDriftingResponsive.p_grandMean.npy'))
        writeNPY(p_compareAcrossStimuli, ...
            fullfile(f, '_ss_gratingsDriftingResponsive.p_acrossStimuli.npy'))
        io.writeTuningResults(dirTuning, oriTuning, f, 'gratingsDrifting');
    end
end