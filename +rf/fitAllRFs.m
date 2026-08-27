function [rfGaussPars, fitGaussians, fitWeights, peakNoiseRatio, ...
    bestSubFields, subFieldSigns, sizeTuning] = ...
    fitAllRFs(rFields, rfBins, gridX, gridY, zTraces)
%FITALLRFS   Fit 2D Gaussian to each mapped RF.

% INPUTS
% rFields           [rows x columns x time x ON/OFF x units], mapped RF for
%                   each unit and subfield
% rfBins            1 x RFframes]; time of spatio-temporal RF relative
%                   to response in number of sitmulus frames
% gridX             Horizontal (along azimuth) locations of pixel centres
% gridY             Vertical (along elevation) locations of pixel centres
% zTraces           [t x ROIs]; z-scored calcium traces of units, sampled
%                   at stimulus presentation times

% OUTPUTS
% rfGaussPars       [ROIs x 7], parameters of fitted Gaussians: amplitude, 
%                   xCenter, xStd, yCenter, yStd, rotation, offset
% fitGaussians      [ROIs x rows x columns x ON/OFF], fitted Gaussians
%                   sampled at same loctions as mapped RFs
% fitWeights        [ROIs x RFframes], fitted weights of 2D Gaussians to
%                   describe temporal dimension of spatio-temporal RFs
% peakNoiseRatio    [ROIs], peak height of RF in SDs of residuals between
%                   mapped RF and fitted Gaussian
% bestSubFields     [ROIs], 1: ON, 2: OFF, 3: ON+OFF
% subFieldSigns     [ROIs x 2], 1: enhanced response, -1: suppressed
%                   response
% sizeTuning        [ROIs x sizes], size tuning curve, only in response to
%                   circle stimuli, otherwise NaN

thresh_diam = 0.75;

numUnits = size(zTraces, 2);
% parameters of fitted Gaussian:
% [amplitude, xCenter, xStd, yCenter, yStd, rotation, offset]
rfGaussPars = NaN(numUnits, 7);
% fitted 2D Gaussian map (one for ON and OFF)
fitGaussians = NaN(numUnits, size(rFields,1), size(rFields,2), 2);
fitWeights = NaN(numUnits, length(rfBins));
peakNoiseRatio = NaN(numUnits, 1);
bestSubFields = NaN(numUnits, 1);
subFieldSigns = NaN(numUnits, 2);

if ndims(rFields) == 5 % visual noise
    sizeTuning = NaN;
elseif ndims(rFields) == 6 % circle paradigm
    sizeTuning = NaN(numUnits, size(rFields,3));
end
for iUnit = 1:numUnits
    if ndims(rFields) == 5 % visual noise
        % rfield: [rows x cols x t x ON/OFF]
        rfield = rFields(:,:,:,:,iUnit);
        % invert polarity of OFF field so that positive values
        % mean: unit is driven by black square
        rfield(:,:,:,2) = -rfield(:,:,:,2);
    elseif ndims(rFields) == 6 % circle paradigm
        % continue with good circle sizes (average across sizes that drive 
        % unit well)
        % rfield: [rows x cols x diameter x t x ON/OFF]
        rfield = rFields(:,:,:,:,:,iUnit);
        % invert polarity of OFF field so that positive values
        % mean: unit is driven by black square
        rfield(:,:,:,:,2) = -rfield(:,:,:,:,2);
        % average across time
        rf_tmp = squeeze(mean(rfield,4,"omitnan"));
        % most driven RF "pixel", regardless of whether response is
        % enhanced or suppressed
        [~,mpx] = max(abs(rf_tmp), [], "all", "omitnan");
        % translate "pixel" to indices in RF dimensions
        [mr,mc,md,ms] = ind2sub(size(rf_tmp), mpx);
        % get size tuning for optimal RF location and subfield
        mx = squeeze(rf_tmp(mr,mc,:,ms))';
        % sign of strongest response (enhanced or suppressed)
        sgn = sign(mx(md));
        sizeTuning(iUnit,:) = mx;
        % determine circle sizes that drive or suppress cell to
        % >thresh_diam of the sign-aligned peak response
        mx_signed = mx .* sgn;
        gd = mx_signed ./ max(mx_signed, [], 2, 'omitnan') > thresh_diam;
        % average across well-driving circle sizes
        rfield = squeeze(mean(rfield(:,:,gd,:,:),3));
    end
    % skip if RF is invalid (all NaNs)
    if all(isnan(rfield), "all")
        continue
    end
    % find best subfield (combination): find whether Gaussian
    % is best fit to only ON, only OFF, or ON-OFF subfields,
    % and optimal sign of each subfield.
    % average across time
    rf_tmp = squeeze(mean(rfield,3)); % [rows x cols x ON/OFF]
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

    % amplitudes (weights) of spatial RF across time span of RF
    weights = reshape(fitRFs(:,:,:,bestSubField), [], 1) \ ...
        reshape(permute(rfield, [1 2 4 3]), [], size(rfield,3));
    fitWeights(iUnit,:) = weights;
end
