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

% colors for plotting
cols = [1 0 0; 1 0.8 0.8];
phaseCols = lines(6);

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
        stimDurs = diff(stimTimes,1,2);
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
        invalidTrials = find(stimDurs < stimulusDur - durationSlack);
        for iUnit = 1:length(units)
            t = spikeData.times(spikeData.clusters == units(iUnit));
            [sp_aligned, trial] = events.alignData(t, ...
                stimData.times(:,1), [baselineTime stimulusDur]);
            for stim = 1:length(stimDirs)
                stimTrials = find(stimData.ids == stim);
                for rep = 1:length(stimTrials)
                    tr = stimTrials(rep);
                    if ismember(tr, invalidTrials)
                        continue
                    end
                    ind_spikes = trial == tr;
                    amplitudes(rep, stim, iUnit) = ...
                        sum(sp_aligned(ind_spikes) >= 0) / ...
                        min(stimDur(tr), stimulusDur) - ...
                        sum(sp_aligned(ind_spikes) , 0) / (-baselineTime);
                end
            end
            % test whether unit is responsive
            amps = amplitudes(:);
            groups = reshape(repmat(1:length(stimDirs), stimReps), [], 1);
            invalid = isnan(amps);
            amps(invalid) = [];
            groups(invalid) = [];
            T = table(amps, categorical(groups), ...
                'VaraibleNames', {'amps', 'stim'});
            lme = fitlme(T, 'amps ~ 1 + (1|stim)', 'FitMethod', 'REML');

            % CONTINUE HERE
        end


        % initialize results structures
        dirTuning = struct('preference', ...
            repmat({NaN(1,1)}, size(kernel,2), 1), ...
            'selectivity', NaN, 'pValue', NaN, 'responseSign', NaN);
        oriTuning = struct('preference', ...
            repmat({NaN(1,1)}, size(kernel,2), 1), ...
            'selectivity', NaN, 'pValue', NaN, 'responseSign', NaN);
        indResponsive = find(p_kernel < maxP & r2_kernel > minR2)';
        % loop over all ROIs
        for iCell = indResponsive
            % determine direction/orientation preference,
            % selectivity, and signficance (response sign: -1 if
            % most responses are negative)
            [drct, ori] = tuning.determineTuning(amplitudes(:,:,iCell), ...
                stimDirs, numShuffles);
            dirTuning(iCell) = drct;
            oriTuning(iCell) = ori;
            maxi = max(amplitudes(:,:,iCell),[],'all');
            mini = min(amplitudes(:,:,iCell),[],'all');

            % plot results
            % replicate responses at 0 degrees to plot at 360 or
            % 180 degrees
            switch type
                case {'gratingsDrifting', 'bars'}
                    ampsDir = amplitudes(:, [1:end,1], iCell) .* drct.responseSign;
                    ampsOri = [amplitudes(:,1:length(stimDirs)/2, iCell); ...
                        amplitudes(:,length(stimDirs)/2+1:end, iCell)] .* ...
                        ori.responseSign;
                    ampsOri = ampsOri(:, [1:end 1]);
                    figure('Position', [100 380 1740 420]);
                    tiledlayout(1,3)
                case 'gratingsStatic'
                    ampsOri = amplitudes(:, [1:end,find(stimDirs==0)'], ...
                        iCell) .* ori.responseSign;
                    figure('Position', [100 380 1200 420]);
                    tiledlayout(1,2)
            end
            nexttile % response kernel
            plot(time_kernel, kernel(:,iCell), 'k', 'LineWidth', 2)
            hold on
            plot([1 1].*stimDur, [min([0 min(kernel(:,iCell))]) 1], ...
                'k:')
            xlim(time_kernel([1 end]))
            title(sprintf('Kernel (R2 = %.2f, p = %.3f)', ...
                r2_kernel(iCell), p_kernel(iCell)))
            xlabel('Time (s)')
            set(gca, 'box', 'off')

            % direction tuning curve (not if grating was static)
            if any(strcmp(type, stimTypes(1:2)))
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
            end

            nexttile % orientation tuning curve
            c = cols(1,:);
            h = [0 0 0];
            if ori.pValue > maxP
                c = cols(2,:);
            end
            if strcmp(type, stimTypes{3}) % for static gratings
                hold on
                % color code responses for difference phases
                for ph = 1:length(phases)
                    indSt = stimPhasesCirc == phases(ph);
                    hs = plot(stimOrisCirc(indSt)' + ...
                        randn(size(ampsOri(:,indSt))) .* (stimDist/10), ...
                        ampsOri(:,indSt), '.', 'MarkerSize', 10, ...
                        'Color', phaseCols(ph,:));
                end
                h(1) = hs(1);
                % average responses for same orientation across
                % phases
                meanAmps = NaN(length(orientations),1);
                for st = 1:length(orientations)
                    indSt = stimOrisCirc == orientations(st);
                    meanAmps(st) = mean(ampsOri(:,indSt), 'all', 'omitnan');
                end
                h(2) = plot(orientations, meanAmps, ...
                    'k.-', 'MarkerSize', 30, 'LineWidth', 1);
                h(3) = plot(ori.preference, maxi, 'v', 'MarkerSize', 8, ...
                    'MarkerFaceColor', c, 'MarkerEdgeColor', 'none');
                legend(h, 'single trials per phase','median','preference', ...
                    'Location', 'bestoutside')
                set(gca, 'box', 'off', 'XTick', orientations)
            else
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
            end
            xlim([-10 190])
            ylim([min([0 1.1*mini]) 1.1*maxi])
            xlabel('Orientation (deg)')
            title(sprintf('OS = %.3f, p = %.3f', ori.selectivity, ...
                ori.pValue))

            % save and close plot
            saveas(gcf, fullfile(fPlots, sprintf('Unit%03d.jpg', iCell)))
            close(gcf)
        end
        % save data
        io.writeTuningResults(dirTuning, oriTuning, f, type);
    end
end