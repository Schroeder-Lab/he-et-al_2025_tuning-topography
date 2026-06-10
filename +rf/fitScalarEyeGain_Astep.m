function [sHat, out] = fitScalarEyeGain_Astep( ...
    toeplitz, traces, eyePos, rfPars, bestSubFields, subFieldSigns, ...
    gridX, gridY, rfBins, timeWeights, opts)
% fitScalarEyeGain_Astep
%
% Optimizes A = s*I given fixed Gaussian RFs.
%
% INPUTS
%   stim   : [ny x nx x nFrames] stimulus movie.
%            Use signed contrast: gray = 0, white = +1, black = -1.
%
%   resp   : [nFrames x nCells] neural responses aligned to stimulus frames.
%
%   pupil  : [nFrames x 2] pupil position in video coordinates.
%            Should be aligned to stimulus frames.
%
%   rf     : [nCells x 7]
%            Columns:
%              1 amplitude
%              2 xCenter
%              3 xStd
%              4 yCenter
%              5 yStd
%              6 rotation
%              7 offset
%
%   xDeg   : [1 x nx] or [ny x nx] x coordinates in visual degrees.
%   yDeg   : [ny x 1] or [ny x nx] y coordinates in visual degrees.
%
%   lags   : [1 x nLags] stimulus-response lags in frames.
%            response(t) is predicted from stim(t-lag).
%            Use lags = 0 if you have no temporal kernel.
%
%   timeW  : [nCells x nLags] temporal weights.
%            Use ones(nCells,1) if lags = 0.
%
%   opts   : optional struct.
%
% OUTPUTS
%   sHat   : fitted scalar gain. Units: visual deg per pupil unit.
%   out    : diagnostics and final prediction.

% -------------------------
% Defaults
% -------------------------

if nargin < 9 || isempty(rfBins)
    rfBins = 0;
end

if nargin < 10 || isempty(timeWeights)
    timeWeights = ones(size(rfPars,1), numel(rfBins));
end

if nargin < 11 || isempty(opts)
    opts = struct();
end

opts = setDefault(opts, 'sBounds', [-5 5]);          % deg per pupil unit
opts = setDefault(opts, 'nGrid', 41);                % coarse grid before fminbnd
opts = setDefault(opts, 'TolX', 1e-3);
opts = setDefault(opts, 'MaxIter', 80);
opts = setDefault(opts, 'Display', 'iter');

opts = setDefault(opts, 'thetaInDegrees', false);    % true if rf(:,6) is degrees
opts = setDefault(opts, 'centerPupil', false);
opts = setDefault(opts, 'pupilSign', [-1 -1]);        % e.g. [1 -1] if video y is inverted
opts = setDefault(opts, 'shiftSign', 1);             % try +1 and -1 if sign convention unclear

opts = setDefault(opts, 'useSpatialOffset', false); % usually false
opts = setDefault(opts, 'fitResponseIntercept', false);
opts = setDefault(opts, 'fitCellGain', true);       % profile out gain + intercept per cell
opts = setDefault(opts, 'normalizeCells', false);    % prevents high-variance cells dominating

% -------------------------
% Dimensions and checks
% -------------------------

[nFrames, nPix] = size(toeplitz);
[nRespFrames, nCells] = size(traces);

assert(nRespFrames == nFrames, 'resp must be [nFrames x nCells].');
assert(size(eyePos,1) == nFrames && size(eyePos,2) == 2, ...
    'pupil must be [nFrames x 2].');
assert(size(rfPars,1) == nCells && size(rfPars,2) >= 7, ...
    'rf must be [nCells x 7].');

rfBins = rfBins(:)';
nLags = numel(rfBins);
nPix = nPix / nLags;

if isvector(timeWeights)
    timeWeights = timeWeights(:)';
end

if size(timeWeights,1) == 1
    timeWeights = repmat(timeWeights, nCells, 1);
end

assert(size(timeWeights,1) == nCells && size(timeWeights,2) == nLags, ...
    'timeW must be [nCells x nLags] or [1 x nLags].');

% -------------------------
% Prepare coordinate grid
% -------------------------

if isvector(gridX) && isvector(gridY)
    [X, Y] = meshgrid(gridX(:)', gridY(:));
else
    X = gridX;
    Y = gridY;
end

assert(numel(X) == nPix && numel(Y) == nPix, ...
    'xDeg/yDeg must match stimulus spatial dimensions.');

xPix = X(:);
yPix = Y(:);

% -------------------------
% Prepare RF parameters
% -------------------------

amp  = rfPars(:,1);
x0   = rfPars(:,2);
sx   = rfPars(:,3);
y0   = rfPars(:,4);
sy   = rfPars(:,5);
th   = rfPars(:,6);
rfOffset = rfPars(:,7);

if opts.thetaInDegrees
    th = deg2rad(th);
end

sx = max(sx, eps);
sy = max(sy, eps);

ct = cos(th);
st = sin(th);

goodCell = isfinite(amp) & isfinite(x0) & isfinite(y0) & ...
           isfinite(sx) & isfinite(sy) & isfinite(th) & ...
           sx > 0 & sy > 0;

% Replace invalid cells with harmless values; they will be excluded.
amp(~goodCell) = 0;
x0(~goodCell) = 0;
y0(~goodCell) = 0;
sx(~goodCell) = 1;
sy(~goodCell) = 1;
ct(~goodCell) = 1;
st(~goodCell) = 0;
rfOffset(~goodCell) = 0;
timeWeights(~goodCell,:) = 0;

% -------------------------
% Prepare pupil trace
% -------------------------

if opts.centerPupil
    eyePos = eyePos - median(eyePos, 1, 'omitnan');
end

eyePos(:,1) = eyePos(:,1) .* opts.pupilSign(1);
eyePos(:,2) = eyePos(:,2) .* opts.pupilSign(2);

pupilOK = all(isfinite(eyePos), 2);

% -------------------------
% Valid response samples
% -------------------------

validMask = true(nFrames, nCells) & isfinite(traces);
validMask(:, ~goodCell) = false;

% -------------------------
% Sparse representation of stimulus movie
% This is much faster for sparse noise.
% -------------------------

nzIdx = cell(nFrames, 1);
nzVal = cell(nFrames, 1);

for k = 1:nFrames
    idx = find(isfinite(toeplitz(k,:)) & toeplitz(k,:) ~= 0);
    nzIdx{k} = idx;
    nzVal{k} = toeplitz(k,idx);
end

% -------------------------
% Objective function
% -------------------------

objFun = @(s) objectiveForS(s);

% Coarse grid first. This is useful because sparse stimuli can make
% the objective mildly irregular.
sGrid = linspace(opts.sBounds(1), opts.sBounds(2), opts.nGrid);
sseGrid = nan(size(sGrid));

for j = 1:numel(sGrid)
    sseGrid(j) = objFun(sGrid(j));
end

[~, bestIdx] = min(sseGrid);

if bestIdx == 1
    lo = sGrid(1);
    hi = sGrid(2);
elseif bestIdx == numel(sGrid)
    lo = sGrid(end-1);
    hi = sGrid(end);
else
    lo = sGrid(bestIdx-1);
    hi = sGrid(bestIdx+1);
end

% Refine with bounded scalar optimization.
optimOpts = optimset( ...
    'TolX', opts.TolX, ...
    'MaxIter', opts.MaxIter, ...
    'Display', opts.Display);

[sHat, fval, exitflag, optimOutput] = fminbnd(objFun, lo, hi, optimOpts);

[predHat, predInfo] = predictForS(sHat);

out = struct();
out.sHat = sHat;
out.sse = fval;
out.pred = predHat;
out.predInfo = predInfo;
out.sGrid = sGrid;
out.sseGrid = sseGrid;
out.bracket = [lo hi];
out.exitflag = exitflag;
out.optimOutput = optimOutput;
out.validMask = validMask;
out.opts = opts;

% ============================================================
% Nested function: objective
% ============================================================

    function mse = objectiveForS(s)

        [pred, ~] = predictForS(s);

        resid = pred - traces;
        resid(~validMask) = NaN;

        if opts.normalizeCells
            cellScale = std(traces, 0, 1, 'omitnan');
            cellScale(~isfinite(cellScale) | cellScale <= eps) = 1;
            resid = resid ./ cellScale;
        end

        r = resid(validMask);
        r = r(isfinite(r));

        if isempty(r)
            mse = Inf;
        else
            mse = mean(r.^2);
        end
    end

% ============================================================
% Nested function: prediction for a given scalar s
% ============================================================

    function [pred, info] = predictForS(s)

        pred = zeros(nFrames, nCells);

        % Loop over stimulus frames.
        for frame = 1:nFrames

            if ~pupilOK(frame)
                continue
            end

            index = nzIdx{frame};
            vals = nzVal{frame};

            if isempty(index)
                continue
            end

            % to index into the RFs, separate the time delayed frames in 
            % each row of the toeplitz matrix
            lagInds = ceil(index / nPix);
            index = mod(index - 1, nPix) + 1;
            vals = vals .* timeWeights(lagInds);

            % Coordinates and values of active stimulus pixels/squares.
            xp = xPix(index)';       % [1 x nActive]
            yp = yPix(index)';       % [1 x nActive]
            vals = vals(:);        % [nActive x 1]

            % Shift RF centers according to A = s*I.
            dxEye = opts.shiftSign * s * eyePos(frame,1);
            dyEye = opts.shiftSign * s * eyePos(frame,2);

            cx = x0 + dxEye;       % [nCells x 1]
            cy = y0 + dyEye;       % [nCells x 1]

            % Difference between active stimulus locations and shifted RF centers.
            dx = xp - cx;          % [nCells x nActive]
            dy = yp - cy;          % [nCells x nActive]

            % Rotate coordinates into RF's principal axes.
            xr =  ct .* dx + st .* dy; % [nCells x nActive]
            yr = -st .* dx + ct .* dy; % [nCells x nActive]

            % Evaluate Gaussian RF at active stimulus locations.
            G = exp(-0.5 .* ((xr ./ sx).^2 + (yr ./ sy).^2)); % [nCells x nActive]

            % Spatial RF value.
            F = amp .* G;          % [nCells x nActive]

            % Usually I would leave this false.
            % The Gaussian offset from fitting RF images is often not a useful
            % response baseline in a stimulus convolution model.
            if opts.useSpatialOffset
                F = F + rfOffset;
            end

            % Dot product between sparse stimulus and shifted RF.
            ind = bestSubFields == 1; % ON only
            v_on = vals;
            v_on(v_on < 0) = 0;
            pred(frame, ind) = subFieldSigns(ind,1) .* F(ind,:) * v_on;

            ind = bestSubFields == 2; % OFF only
            v_off = vals;
            v_off(v_off > 0) = 0;
            pred(frame, ind) = subFieldSigns(ind,2) .* F(ind,:) * -v_off;

            ind = bestSubFields == 3; % ON & OFF
            pred(frame, ind) = ...
                subFieldSigns(ind,1) .* F(ind,:) * v_on + ...
                subFieldSigns(ind,2) .* F(ind,:) * -v_off;
        end

        info = struct();
        info.cellGain = ones(1, nCells);
        info.cellIntercept = zeros(1, nCells);

        % Optional: analytically fit gain and intercept per cell for this s.
        % This profiles out amplitude/baseline mismatch and makes the A-step
        % focus more on timing/spatial alignment.
        if opts.fitCellGain

            for iCell = 1:nCells
                m = validMask(:,iCell) & ...
                    isfinite(traces(:,iCell)) & ...
                    isfinite(pred(:,iCell));

                if nnz(m) < 3
                    pred(:,iCell) = NaN;
                    continue
                end

                Xfit = [pred(m,iCell), ones(nnz(m),1)];
                beta = Xfit \ traces(m,iCell);

                info.cellGain(iCell) = beta(1);
                info.cellIntercept(iCell) = beta(2);

                pred(:,iCell) = beta(1) .* pred(:,iCell) + beta(2);
            end

        elseif opts.fitResponseIntercept

            for iCell = 1:nCells
                m = validMask(:,iCell) & ...
                    isfinite(traces(:,iCell)) & ...
                    isfinite(pred(:,iCell));

                if nnz(m) < 1
                    pred(:,iCell) = NaN;
                    continue
                end

                b = mean(traces(m,iCell) - pred(m,iCell), 'omitnan');
                info.cellIntercept(iCell) = b;

                pred(:,iCell) = pred(:,iCell) + b;
            end
        end
    end
end

% ============================================================
% Helper: set default options
% ============================================================

function opts = setDefault(opts, name, value)
if ~isfield(opts, name) || isempty(opts.(name))
    opts.(name) = value;
end
end