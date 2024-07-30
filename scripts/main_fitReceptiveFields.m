% Given the calcium responses to the visual noise stimulus, fit spatial
% receptive fields (ON and OFF fields).

%% Folders
getFolders;

%% Parameters
% for correcting baseline drifts of calcium traces at start of experiments
decayWin = 20; % in s, window to test whether baseline is higher than normal
decayThresh = 1.5; % in std, threshold for decay
correctWin = 150; % in s, window to fit exponential

% for receptive field estimates
% used for fitting 2 RFs (ON and OFF simultaneously), and fitting running
% kernels and RFs simultaneously
lambdasStim = 0.002; % 1e-4;
RFlimits = [0.2 0.4];
RFtypes = {'ON', 'OFF', 'ON+OFF'};

% % for evaluation of receptive fields (significance/goodness)
minEV = 0.005;
minPeak = 3.5;

% colormaps
[cm_ON, cm_OFF] = colmaps.getRFMaps;
cms = cat(3, cm_ON, cm_OFF);

titles = {'ON field','OFF field'};

%% Add paths
addpath(genpath(fullfile(folders.tools, 'npy-matlab')))
addpath(fullfile(folders.repo))

%% Fit RFs and get cross-validated explained variance
sets = {'neurons', 'boutons'};
for s = 1:2 % neurons and boutons
    subjDirs = dir(fullfile(folders.data, sets{s}, 'SS*'));
    for subj = 1:length(subjDirs)
        name = subjDirs(subj).name;
        fprintf('%s\n', name)
        dateDirs = dir(fullfile(folders.data, sets{s}, name, '2*'));
        for dt = 1:length(dateDirs) %dates
            date = dateDirs(dt).name;
            fprintf('  %s\n', date)
            f = fullfile(folders.data, sets{s}, name, date);

            if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
                continue
            end

            fPlots = fullfile(folders.plots, 'ReceptiveFields', ...
                sets{s}, name, date);
            if ~isfolder(fPlots)
                mkdir(fPlots)
            end

            % load data
            data = io.getCalciumData(f);
            time_tr = data.time;
            tr = data.traces;
            planes = data.planes;
            cellIDs = data.ids;
            delays = data.delays;
            data = io.getVisNoiseInfo(f);
            time_stim = data.times;
            stimMaps = data.frames;
            stimSeq = data.stimOrder;
            stimPos = data.edges; % [left right top bottom]
            % flip sign of stimulus borders along y-axis -> positive numbers
            % are above horizon/monitor centre
            stimPos(3:4) = -stimPos(3:4);

            % if stimulus covers both hemifields, only
            % consider right hemifield (affected datasets all recorded in left
            % hemisphere)
            squW = diff(stimPos(1:2)) / size(stimMaps,3);
            squH = -diff(stimPos(3:4)) / size(stimMaps,2);
            if stimPos(1) * stimPos(2) < 0
                % determine left edge of all pixel columns
                leftEdges = stimPos(1) + (0:size(stimMaps,3)-1) .* squW;
                validPix = leftEdges >= 0;
                stimMaps = stimMaps(:,:,validPix);
                stimPos(1) = leftEdges(find(validPix,1));
            end
            [x0, y0] = meshgrid(linspace(stimPos(1), stimPos(2), 100), ...
                flip(linspace(stimPos(4), stimPos(3), 100)));

            % ignore data outside visual noise stimulation
            t_ind = time_tr > time_stim(1) - 10 & ...
                time_tr < time_stim(end) + 10;
            tr = tr(t_ind,:);
            time_tr = time_tr(t_ind);
            % interpolate calcium traces to align all to same time
            [tr, time_tr] = traces.alignPlaneTraces(tr, time_tr, delays, planes);

            % remove strong baseline decay at start of experiment in cells that
            % show it
            tr = traces.removeInitialDecay(tr, time_tr, decayWin, decayThresh);

            % more stimulus information
            stimFrames = stimMaps(stimSeq,:,:);
            stimFrameDur = median(diff(time_stim));
            RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
                ceil(RFlimits(2) / stimFrameDur);
            clear stimMaps stimSeq

            % map RF
            rFields = ...
                rf.getReceptiveField( ...
                tr, time_tr, stimFrames, time_stim, ...
                RFtimesInFrames, lambdasStim);

            % fit Gaussian
            fitPars = NaN(size(tr,2), 6); % each row: [amplitude, xCenter,
            % xStd, yCenter, yStd, rotation]
            fitGaussians = NaN(size(tr,2), size(rFields,1), size(rFields,2), 2);
            fitWeights = NaN(size(tr,2), length(RFtimesInFrames));
            peak_from_noise = NaN(size(tr,2), 1);
            subFields = cell(size(tr,2), 1);
            [~, t_pred] = rf.makeStimToeplitz(stimFrames, ...
                time_stim, RFtimesInFrames);
            predictions = NaN(length(t_pred), size(tr,2));
            EVs = NaN(size(tr,2), 1);
            for iUnit = 1:size(fitPars,1)
                rfield = permute(rFields(:,:,:,:,iUnit), [1 2 4 3]); % [rows x cols x ON/OFF x t]
                [mx,mxTime] = max(max(reshape(abs(rfield), [], size(rfield,4)), [], 1));
                rf_T = squeeze(rfield(:,:,:,mxTime));

                % find whether Gaussian is best fit to only ON, only OFF, or
                % ON-OFF subfields
                [fitRFs, RFsigns, MSEs] = rf.findRFGaussianMask(rf_T);
                [~, bestSubField] = min(MSEs);
                subFields{iUnit} = RFtypes{bestSubField};
                fitGaussians(iUnit,:,:,:) = fitRFs(:,:,:,bestSubField);

                if bestSubField < 3
                    rf_sub = rf_T(:,:,bestSubField) .* RFsigns(bestSubField);
                else
                    rf_sub = (rf_T(:,:,1) .* RFsigns(1) + ...
                        rf_T(:,:,2) .* RFsigns(2)) ./ 2;
                end

                % interpolate RF so that pixels are square with edge length of
                % 1 visual degree
                rf_visDeg = rf.interpolateRFtoVisualDegrees(rf_sub, stimPos);
                % fit Gaussian
                fitPars(iUnit,:) = rf.fit2dGaussRF(rf_visDeg, false);

                % compare peak of fitted RF to amplitude of noise in RF map;
                % reproduce fitted Gaussian RF
                fitRF = rf.sampleGaussianAtPixelPos(rf_sub, stimPos, fitPars(iUnit,:));

                % subtract Gaussian from original RF map
                noise = rf_sub - fitRF;
                % distance of peak from noise
                peak_from_noise(iUnit) = (fitPars(iUnit,1) - ...
                    mean(noise(:))) / std(noise(:));

                % predict response from RF
                weights = reshape(fitRFs(:,:,:,bestSubField), [], 1) \ ...
                    reshape(rfield, [], size(rfield,4));
                spatTempMask = reshape(fitRFs(:,:,:,bestSubField), [], 1) * weights;
                spatTempMask = permute(reshape(spatTempMask, size(fitRFs,1), ...
                    size(fitRFs,2), 2, length(weights)), [1 2 4 3]); % [rows x cols x t x ON/OFF]
                [predictions(:, iUnit), EVs(iUnit)] = rf.predictFromRF(tr(:,iUnit), time_tr, stimFrames, ...
                    time_stim, RFtimesInFrames, spatTempMask);
                fitWeights(iUnit,:) = weights;

                % transform RF position to absolute values (relative to screen)
                fitPars(iUnit,2) = stimPos(1) + squW/2 + fitPars(iUnit,2);
                fitPars(iUnit,4) = stimPos(3) - squH/2 - fitPars(iUnit,4);
                % mirror RF orientation to account for flipped y-axis direction
                % (top is positive)
                fitPars(iUnit,6) = -fitPars(iUnit,6);

                % plot RF
                if EVs(iUnit) >= minEV && peak_from_noise(iUnit) >= minPeak
                    figure('Position', [75 195 1470 475])
                    for sf = 1:2
                        subplot(1,2,sf)
                        imagesc([stimPos(1)+squW/2 stimPos(2)-squW/2], ...
                            [stimPos(3)-squH/2 stimPos(4)+squH/2], ...
                            rf_T(:,:,sf),[-mx mx])
                        hold on
                        fitRF = rf.D2GaussFunctionRot(fitPars(iUnit,:), cat(3, x0, y0));
                        if any(strcmp(subFields{iUnit}, [RFtypes(sf), RFtypes(3)]))
                            contour(x0, y0, fitRF, [1 1] .* (fitPars(iUnit,1) * normpdf(1.5) / normpdf(0)), ...
                                'k', 'LineWidth', 2)
                        end
                        axis image
                        set(gca, 'box', 'off', 'XTick', round(stimPos(1:2)), ...
                            'YTick', [stimPos(4) stimPos(3)], 'YDir', 'normal', ...
                            'YTickLabel', [stimPos(4) stimPos(3)])
                        xlim(round([stimPos(1) stimPos(2)]))
                        colormap(gca, cms(:,:,sf))
                        title(titles{sf})
                        colorbar
                    end
                    sgtitle(sprintf('ROI %d (EV: %.3f, peak/noise: %.1f)', ...
                        iUnit, EVs(iUnit), peak_from_noise(iUnit)))

                    saveas(gcf, fullfile(fPlots, sprintf('Unit%03d.jpg', iUnit)));
                    close gcf
                end
            end

            % save results
            writeNPY(permute(rFields, [5 1 2 3 4]), fullfile(f, '_ss_rf.maps.npy'))
            writecell(subFields, fullfile(f, '_ss_rf.type.csv'))
            writeNPY(fitPars, fullfile(f, '_ss_rf.gaussFitPars.npy'))
            writeNPY(peak_from_noise, fullfile(f, '_ss_rf.peak.npy'))
            writeNPY(fitGaussians, fullfile(f, '_ss_rf.gaussMask.npy'))
            writeNPY(fitWeights, fullfile(f, '_ss_rf.gaussTimeWeights.npy'))
            writeNPY(EVs, fullfile(f, '_ss_rf.explVars.npy'))
            writeNPY(predictions, fullfile(f, '_ss_rfPrediction.traces.npy'))
            writeNPY(t_pred, fullfile(f, '_ss_rfPrediction.timestamps.npy'))
            writeNPY(stimPos, fullfile(f, '_ss_rfDescr.edges.npy'))
            writeNPY(RFtimesInFrames * stimFrameDur, fullfile(f, '_ss_rfDescr.timestamps.npy'))
        end
    end
end

%% Plot example receptive fields for paper
% name = 'SS048';
% date = '2015-09-26';
% ROIs = [29 31 71 137];
% stimAzimuth = -110;
% stimElevation = 15;
% stimSigma = 9;
% RFstd = 1.5;
% RFheight = normpdf(RFstd) / normpdf(0);
% folder = fullfile(folderBase, name, date);
% 
% stimPos = readNPY(fullfile(folder, '_ss_sparseNoiseArea.edges.npy')); % [left right top bottom]
% rFields = readNPY(fullfile(folder, '_ss_rf.maps.npy'));
% fitPars = readNPY(fullfile(folder, '_ss_rf.gaussFitPars.npy'));
% subFields = readcell(fullfile(folder, '_ss_rf.type.csv'));
% 
% stimPos(3:4) = -stimPos(3:4);
% squW = diff(stimPos(1:2)) / size(rFields,3);
% squH = -diff(stimPos(3:4)) / size(rFields,2);
% [x0, y0] = meshgrid(linspace(stimPos(1), stimPos(2), 100), ...
%     flip(linspace(stimPos(4), stimPos(3), 100)));
% for ex = 1:length(ROIs)
%     rfield = squeeze(rFields(ROIs(ex),:,:,:,:)); % [rows x cols x t x ON/OFF]
%     rfield(:,:,:,2) = -rfield(:,:,:,2);
%     [mx,mxTime] = max(max(reshape(permute(abs(rfield), [1 2 4 3]), ...
%         [], size(rfield,3)), [], 1));
%     fitRF = whiteNoise.D2GaussFunctionRot(fitPars(ROIs(ex),:), cat(3, x0, y0));
%     figure('Position', [75 195 1470 475])
%     for f = 1:2
%         subplot(1,2,f)
%         imagesc([stimPos(1)+squW/2 stimPos(2)-squW/2], ...
%             [stimPos(3)-squH/2 stimPos(4)+squH/2], ...
%             rfield(:,:,mxTime,f),[-mx mx])
%         hold on
%         if any(strcmp(subFields{ROIs(ex)}, [RFtypes(f), RFtypes(3)]))
%             contour(x0, y0, fitRF, [1 1] .* (fitPars(ROIs(ex),1) * RFheight), ...
%                 'k', 'LineWidth', 2)
%         end
%         rectangle('Position', [stimAzimuth-stimSigma ...
%             stimElevation-stimSigma 2*stimSigma 2*stimSigma], ...
%             'Curvature', 1, 'LineStyle', '--', 'LineWidth', 2)
%         axis image
%         set(gca, 'box', 'off', 'XTick', round(stimPos(1:2)), ...
%             'YTick', [stimPos(4) 0 stimPos(3)], 'YDir', 'normal', ...
%             'YTickLabel', [stimPos(4) 0 stimPos(3)])
%         xlim(round([stimPos(1) stimPos(2)]))
%         colormap(gca, cms(:,:,f))
%         sgtitle(sprintf('ROI %d', ROIs(ex)))
%         colorbar
%     end
% 
%     filename = fullfile('C:\Users\Sylvia\OneDrive - University of Sussex\Paper\Plots\Fig 1', ...
%         sprintf('RF%03d.eps', ROIs(ex)));
%     myprint
%     close gcf
% end