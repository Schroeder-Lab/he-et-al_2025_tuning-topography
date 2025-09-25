function [rfGaussPars, fitGaussians, fitWeights, peakNoiseRatio, ...
    bestSubFields, subFieldSigns, predictions, EVs] = ...
    fitAllRFs(rFields, rfBins, gridX, gridY, zTraces, toeplitz)

numUnits = size(zTraces, 2);
% parameters of fitted Gaussian:
% [amplitude, xCenter, xStd, yCenter, yStd, rotation]
rfGaussPars = NaN(numUnits, 7);
% fitted 2D Gaussian map (one for ON and OFF)
fitGaussians = NaN(numUnits, size(rFields,1), size(rFields,2), 2);
fitWeights = NaN(numUnits, length(rfBins));
peakNoiseRatio = NaN(numUnits, 1);
bestSubFields = NaN(numUnits, 1);
subFieldSigns = NaN(numUnits, 2);
predictions = NaN(size(toeplitz, 1), numUnits);
EVs = NaN(numUnits, 1);
for iUnit = 1:numUnits
    % rfield: [rows x cols x t x ON/OFF]
    rfield = rFields(:,:,:,:,iUnit);
    % continue if RF is invalid (all NaNs)
    if all(isnan(rfield), "all")
        continue
    end
    % invert polarity of OFF field so that positive values
    % mean: unit is driven by black square
    rfield(:,:,:,2) = -rfield(:,:,:,2);
    % find best subfield (combination): find whether Gaussian
    % is best fit to only ON, only OFF, or ON-OFF subfields,
    % and what optimal sign of each subfield is
    % average across time
    rf_tmp = squeeze(mean(rfield,3));
    % fitRFs: Gaussian masks with best sign (pos or neg)
    [fitRFs, RFsigns, MSEs] = rf.findRFGaussianMask(rf_tmp);
    [~, bestSubField] = min(MSEs);
    bestSubFields(iUnit) = bestSubField;
    subFieldSigns(iUnit,:) = RFsigns;
    fitGaussians(iUnit,:,:,:) = fitRFs(:,:,:,bestSubField);

    if bestSubField < 3
        rf_sub = rf_tmp(:,:,bestSubField) .* RFsigns(bestSubField);
    else
        rf_sub = (rf_tmp(:,:,1) .* RFsigns(1) + ...
            rf_tmp(:,:,2) .* RFsigns(2)) ./ 2;
    end

    % fit Gaussian
    [rfGaussPars(iUnit,:), rf_gauss] = rf.fit2dGaussRF(...
        rf_sub, false, gridX, gridY);
    % mirror RF orientation to account for flipped y-axis direction
    % (top is positive)
    rfGaussPars(iUnit,6) = -rfGaussPars(iUnit,6);

    % subtract Gaussian from original RF map
    noise = rf_sub - rf_gauss;
    % distance of peak from noise
    peakNoiseRatio(iUnit) = rfGaussPars(iUnit,1) / std(noise(:));

    % predict response from RF
    % amplitudes (weights) of spatial RF across time span of RF
    weights = reshape(fitRFs(:,:,:,bestSubField), [], 1) \ ...
        reshape(permute(rfield, [1 2 4 3]), [], size(rfield,3));
    % generate spatio-temporal RF from fitted Gaussian and
    % temporal weights
    spatTempMask = reshape(fitRFs(:,:,:,bestSubField), [], 1) * ...
        weights; % [pix x t]
    % spatTempMas: [rows x cols x t x ON/OFF]
    spatTempMask = permute(reshape(spatTempMask, size(fitRFs,1), ...
        size(fitRFs,2), 2, length(weights)), [1 2 4 3]);
    spatTempMask(:,:,:,2) = -spatTempMask(:,:,:,2);
    % predict calcium trace based on generated spatio-temporal
    % RF
    [predictions(:, iUnit), EVs(iUnit)] = ...
        rf.predictFromRF(zTraces(:,iUnit), toeplitz, ...
        spatTempMask);
    fitWeights(iUnit,:) = weights;
end