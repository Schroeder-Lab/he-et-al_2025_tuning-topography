function [fitRFs, RFsigns, MSEs] =  findRFGaussianMask(rfield)

% fit Gaussian mask to ON field only, OFF field only, or average of ON and
% OFF fields; determine which fit explains most of the RF structure

fitRFs = NaN(size(rfield,1), size(rfield,2), size(rfield,3), 3);
RFsigns = NaN(1, 3);
MSEs = Inf(1, 3);

% start with separate ON and OFF fields
for sub = 1:3
    if sub < 3
        subRF = rfield(:,:,sub);
    else
        subRF = (rfield(:,:,1) .* RFsigns(1) + rfield(:,:,2) .* RFsigns(2)) ./ 2;
    end
    for sign = [1 -1]
        if sub == 3 && sign == -1 % we already found the optimal signs for ON and OFF fields
            continue
        end
        [~, fitRF] = rf.fit2dGaussRF(subRF .* sign, false);
        mask = zeros(size(rfield));
        if sub < 3
            mask(:,:,sub) = fitRF * sign;
        else
            mask(:,:,1) = fitRF * RFsigns(1);
            mask(:,:,2) = fitRF * RFsigns(2);
        end
        mse = sum((rfield(:) - mask(:)) .^ 2) / numel(rfield);
        if mse < MSEs(sub)
            fitRFs(:,:,:,sub) = mask;
            MSEs(sub) = mse;
            RFsigns(sub) = sign;
        end
    end
end