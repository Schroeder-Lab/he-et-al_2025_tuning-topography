function [fitRFs, RFsigns, MSEs] =  findRFGaussianMask(rfield)
%FINDGAUSSIANMASK   Fit Gaussian mask to ON field only (* 1 and (-1)), 
% OFF field only (* 1 and (-1)), and average of ON and OFF fields (using 
% best signs determined before); determine errors of fits.

% INPUTS
% rfield    [rows x columns x ON/OFF], response-triggered ON and OFF receptive
%           field

% OUTPUTS
% fitRFs    [rows x columns x subfield x ON/OFF/ON-OFF], Gaussian masks fitted
%           to each sign of each subfield (including mean of ON and OFF),
%           for ON RF, fitRFs(:,:,2,1) is zero, and accordingly for OFF RF
% RFsigns   [3], 1 or -1, optimal sign of each subfield, if -1 responses
%           are suppressed 
% MSEs      [3], mean-squared-error of each Gaussian fit compared to
%           original RF

fitRFs = NaN(size(rfield,1), size(rfield,2), size(rfield,3), 3);
RFsigns = NaN(1, 3);
MSEs = Inf(1, 3);

for sub = 1:3
    if sub < 3 % start with separate ON and OFF fields
        subRF = rfield(:,:,sub);
    else % for ON-OFF RF
        % take average of original ON and OFF subfield, multiplied by
        % previously found optimal sign (make suppressed responses go
        % positive)
        subRF = (rfield(:,:,1) .* RFsigns(1) + rfield(:,:,2) .* RFsigns(2)) ./ 2;
    end
    for sign = [1 -1]
        if sub == 3 && sign == -1 % we already found the optimal signs for ON and OFF fields
            continue
        end
        % fit Gaussian mask to current subfield multiplied by current sign
        [~, fitRF] = rf.fit2dGaussRF(subRF .* sign, false);
        mask = zeros(size(rfield));
        if sub < 3
            % remember fitted masks multiplied by current sign (to match
            % original RF)
            mask(:,:,sub) = fitRF * sign;
        else
            mask(:,:,1) = fitRF * RFsigns(1);
            mask(:,:,2) = fitRF * RFsigns(2);
        end
        % determine MSE
        mse = sum((rfield(:) - mask(:)) .^ 2) / numel(rfield);
        if mse < MSEs(sub) % for sign resulting in smaller MSE
            % remember fitted Gaussian mask, MSE, and optimal sign
            fitRFs(:,:,:,sub) = mask;
            MSEs(sub) = mse;
            RFsigns(sub) = sign;
        end
    end
end