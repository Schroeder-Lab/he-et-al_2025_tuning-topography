function [predictions, EVs] = makeAllNoisePredictions(zTraces, toeplitz, ...
    fitGaussians, fitWeights, sizeTuning)

% INPUTS
% zTraces           [t x ROIs]; z-scored calcium traces of units, sampled
%                   at stimulus presentation times
% toeplitz          [t x pixels]; noise stimulus
% fitGaussians      [ROIs x rows x columns x ON/OFF], fitted Gaussians
%                   sampled at same loctions as mapped RFs
% fitWeights        [ROIs x RFframes], fitted weights of 2D Gaussians to
%                   describe temporal dimension of spatio-temporal RFs
% (sizeTuning)      [ROIs x sizes], size tuning curve, only in response to
%                   circle stimuli

% OUTPUTS
% predictions       [t x ROIs], predicted calcium traces based on fitted
%                   spatio-temporal RFs and stimulus sequence
% EVs               [ROIs], explained variance of fitted spatio-temporal RF

numUnits = size(zTraces, 2);
predictions = NaN(size(toeplitz, 1), numUnits);
EVs = NaN(numUnits, 1);

for iUnit = 1:numUnits
    % generate spatio-temporal RF from fitted Gaussian and temporal weights
    spatTempMask = reshape(fitGaussians(iUnit,:,:,:), [], 1) * ...
        fitWeights(iUnit,:); % [pix x t]

    if nargin < 5 % visual noise
        % spatTempMask: [rows x cols x t x ON/OFF]
        spatTempMask = permute(reshape(spatTempMask, size(fitGaussians,2), ...
            size(fitGaussians,3), 2, size(fitWeights,2)), [1 2 4 3]);
        spatTempMask(:,:,:,2) = -spatTempMask(:,:,:,2);
    else % circle paradigm
        % repeat spatio-temporal RF for each diameter, weighted according
        % to size tuning
        tuning = sizeTuning(iUnit,:);
        if all(isnan(tuning))
            continue
        end
        [~, peakSize] = max(abs(tuning), [], 2, 'omitnan');
        sgn = sign(tuning(peakSize));
        if sgn == 0 || isnan(sgn)
            sgn = 1;
        end
        tuning = tuning .* sgn;
        peakTuning = max(tuning, [], 2, 'omitnan');
        if peakTuning == 0 || isnan(peakTuning)
            continue
        end
        spatTempMask = reshape(spatTempMask, [], 1) * ...
            (tuning ./ peakTuning); % [(pix * t) x diameters]
        % spatTempMas: [rows x cols x diameters x t x ON/OFF]
        spatTempMask = permute(reshape( ...
            spatTempMask, size(fitGaussians,2), size(fitGaussians,3), ...
            2, size(fitWeights,2), size(sizeTuning,2)), [1 2 5 4 3]);
        spatTempMask(:,:,:,:,2) = -spatTempMask(:,:,:,:,2);
    end

    % predict calcium trace based on generated spatio-temporal
    % RF
    [predictions(:, iUnit), EVs(iUnit)] = ...
        rf.predictFromRF(zTraces(:,iUnit), toeplitz, spatTempMask);
end
