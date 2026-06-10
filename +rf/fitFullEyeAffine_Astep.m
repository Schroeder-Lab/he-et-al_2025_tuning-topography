function [AHat, out] = fitFullEyeAffine_Astep( ...
    toeplitz, traces, eyePos, rfPars, bestSubFields, subFieldSigns, ...
    gridX, gridY, rfBins, timeWeights, opts)
% fitFullEyeAffine_Astep
%
% Optimizes full A matrix given fixed Gaussian RFs:
%
%   [dx; dy] = A * [eyeX; eyeY]
%
% where
%
%   A = [a11 a12;
%        a21 a22]
%
% The shifted RF center for each stimulus/response frame is:
%
%   cx = x0 + shiftSign * dx
%   cy = y0 + shiftSign * dy
%
% INPUTS
%   toeplitz       : [nFrames x nPix*nLags] stimulus design matrix
%   traces         : [nFrames x nCells] neural responses
%   eyePos         : [nFrames x 2] eye position / pupil position
%   rfPars         : [nCells x 7]
%                    columns:
%                       1 amplitude
%                       2 xCenter
%                       3 xStd
%                       4 yCenter
%                       5 yStd
%                       6 rotation
%                       7 offset
%   bestSubFields  : [nCells x 1]
%                    1 = ON only
%                    2 = OFF only
%                    3 = ON and OFF
%   subFieldSigns  : [nCells x 2]
%                    column 1: ON sign
%                    column 2: OFF sign
%   gridX, gridY   : stimulus grid coordinates
%   rfBins         : [1 x nLags] lag offsets in frames
%   timeWeights    : [nCells x nLags], [1 x nLags], or empty
%   opts           : optional struct
%
% OUTPUTS
%   AHat           : fitted 2 x 2 matrix
%   out            : diagnostics and prediction

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

% A parameterization:
% theta = [a11; a12; a21; a22]
if ~isfield(opts, 'A0') || isempty(opts.A0)
    opts.A0 = zeros(2,2);
end

if ~isfield(opts, 'ABounds') || isempty(opts.ABounds)
    if isfield(opts, 'sBounds') && ~isempty(opts.sBounds)
        opts.ABounds = repmat(opts.sBounds(:)', 4, 1);
    else
        opts.ABounds = repmat([-5 5], 4, 1);
    end
end

opts = setDefault(opts, 'linearA', []);
opts = setDefault(opts, 'linearb', []);

opts = setDefault(opts, 'TolX', 1e-4);
opts = setDefault(opts, 'TolFun', 1e-6);
opts = setDefault(opts, 'MaxIter', 150);
opts = setDefault(opts, 'MaxFunEvals', 1000);
opts = setDefault(opts, 'Display', 'iter');

% Optional multistart. Usually start with nStarts = 1 in the alternating loop.
opts = setDefault(opts, 'nStarts', 1);
opts = setDefault(opts, 'startScale', 0.25);
opts = setDefault(opts, 'randomSeed', 1);

opts = setDefault(opts, 'thetaInDegrees', false);
opts = setDefault(opts, 'centerPupil', false);
opts = setDefault(opts, 'pupilSign', [1 1]);
opts = setDefault(opts, 'shiftSign', 1);

% Important for toeplitz design matrices:
% If row t contains stimulus frames t-rfBins(q), then the eye position
% for lag q should be eyePos(t-rfBins(q),:), not eyePos(t,:).
opts = setDefault(opts, 'useLagSpecificEyePos', true);
opts = setDefault(opts, 'lagSign', -1);  % eyeFrame = frame + lagSign * rfBins(q)

opts = setDefault(opts, 'useSpatialOffset', false);
opts = setDefault(opts, 'fitResponseIntercept', false);
opts = setDefault(opts, 'fitCellGain', true);
opts = setDefault(opts, 'normalizeCells', false);

% Optional weak regularization of A.
% Useful if full A becomes unstable or off-diagonal terms become implausibly large.
opts = setDefault(opts, 'lambdaA', 0);
opts = setDefault(opts, 'APrior', opts.A0);

% Penalty used only for fminsearch if it tries to leave the bounds.
opts = setDefault(opts, 'boundPenaltyScale', 1e6);

% -------------------------
% Dimensions and checks
% -------------------------

[nFrames, nCols] = size(toeplitz);
[nRespFrames, nCells] = size(traces);

assert(nRespFrames == nFrames, 'traces must be [nFrames x nCells].');
assert(size(eyePos,1) == nFrames && size(eyePos,2) == 2, ...
    'eyePos must be [nFrames x 2].');
assert(size(rfPars,1) == nCells && size(rfPars,2) >= 7, ...
    'rfPars must be [nCells x 7].');

bestSubFields = bestSubFields(:);
assert(numel(bestSubFields) == nCells, ...
    'bestSubFields must have one entry per cell.');
assert(size(subFieldSigns,1) == nCells && size(subFieldSigns,2) >= 2, ...
    'subFieldSigns must be [nCells x 2].');

rfBins = rfBins(:)';
nLags = numel(rfBins);

nPix = nCols / nLags;
assert(abs(nPix - round(nPix)) < 1e-9, ...
    'Number of toeplitz columns must be nPix * nLags.');
nPix = round(nPix);

% -------------------------
% Prepare time weights
% -------------------------

if isempty(timeWeights)
    timeWeights = ones(nCells, nLags);
elseif isvector(timeWeights)
    if numel(timeWeights) == nLags
        timeWeights = repmat(timeWeights(:)', nCells, 1);
    elseif numel(timeWeights) == nCells && nLags == 1
        timeWeights = timeWeights(:);
    else
        error(['Vector timeWeights must have length nLags, or length nCells ' ...
               'when nLags == 1.']);
    end
elseif size(timeWeights,1) == 1 && size(timeWeights,2) == nLags
    timeWeights = repmat(timeWeights, nCells, 1);
end

assert(size(timeWeights,1) == nCells && size(timeWeights,2) == nLags, ...
    'timeWeights must be [nCells x nLags] or [1 x nLags].');

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
    'gridX/gridY must match the number of spatial stimulus pixels.');

xPix = X(:);
yPix = Y(:);

% -------------------------
% Prepare RF parameters
% -------------------------

amp      = rfPars(:,1);
x0       = rfPars(:,2);
sx       = rfPars(:,3);
y0       = rfPars(:,4);
sy       = rfPars(:,5);
th       = rfPars(:,6);
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

% Replace invalid cells with harmless values; exclude them from scoring.
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
% Prepare eye trace
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

if opts.useLagSpecificEyePos
    validT = true(nFrames, 1);

    for q = 1:nLags
        lagOffset = round(rfBins(q));
        eyeFrame = (1:nFrames)' + opts.lagSign * lagOffset;

        ok = eyeFrame >= 1 & eyeFrame <= nFrames;
        thisValid = false(nFrames, 1);
        thisValid(ok) = pupilOK(eyeFrame(ok));

        validT = validT & thisValid;
    end
else
    validT = pupilOK;
end

validMask = repmat(validT, 1, nCells) & isfinite(traces);
validMask(:, ~goodCell) = false;

if opts.normalizeCells
    cellScale = std(traces, 0, 1, 'omitnan');
    cellScale(~isfinite(cellScale) | cellScale <= eps) = 1;
else
    cellScale = ones(1, nCells);
end

% -------------------------
% Sparse representation of stimulus design matrix
% -------------------------

nzIdx = cell(nFrames, 1);
nzVal = cell(nFrames, 1);

for frame = 1:nFrames
    [~, idx, vals] = find(toeplitz(frame,:));

    idx = idx(:);
    vals = full(vals(:));

    keep = isfinite(vals) & vals ~= 0;
    nzIdx{frame} = idx(keep);
    nzVal{frame} = vals(keep);
end

% -------------------------
% Optimize A
% -------------------------

theta0 = packA(opts.A0);
thetaPrior = packA(opts.APrior);

lb = opts.ABounds(:,1);
ub = opts.ABounds(:,2);

assert(numel(theta0) == 4, 'opts.A0 must be 2 x 2.');
assert(all(size(opts.ABounds) == [4 2]), ...
    'opts.ABounds must be [4 x 2], in order [a11; a12; a21; a22].');

theta0 = min(max(theta0, lb), ub);
thetaStarts = makeThetaStarts(theta0, lb, ub, opts);

bestF = Inf;
bestTheta = theta0;
bestExitflag = NaN;
bestOptimOutput = struct();
startResults = struct([]);

for iStart = 1:size(thetaStarts, 2)
    thisTheta0 = thetaStarts(:, iStart);

    optimOpts = optimset( ...
        'TolX', opts.TolX, ...
        'TolFun', opts.TolFun, ...
        'MaxIter', opts.MaxIter, ...
        'MaxFunEvals', opts.MaxFunEvals, ...
        'Display', opts.Display);

    [thetaFit, fval, exitflag, optimOutput] = fmincon( ...
        @objectiveForTheta, thisTheta0, ...
        opts.linearA, opts.linearb, [], [], lb, ub, [], optimOpts);

    startResults(iStart).theta0 = thisTheta0;
    startResults(iStart).thetaFit = thetaFit;
    startResults(iStart).A = unpackA(thetaFit);
    startResults(iStart).fval = fval;
    startResults(iStart).exitflag = exitflag;
    startResults(iStart).optimOutput = optimOutput;

    if fval < bestF
        bestF = fval;
        bestTheta = thetaFit;
        bestExitflag = exitflag;
        bestOptimOutput = optimOutput;
    end
end

AHat = unpackA(bestTheta);

[predHat, predInfo] = predictForA(AHat);

out = struct();
out.sse = bestF;
out.pred = predHat;
out.predInfo = predInfo;
out.exitflag = bestExitflag;
out.optimOutput = bestOptimOutput;
out.startResults = startResults;
out.validMask = validMask;
out.opts = opts;

% ============================================================
% Nested function: bounded objective for fminsearch fallback
% ============================================================

    function val = objectiveWithBounds(theta)

        theta = theta(:);

        below = max(lb - theta, 0);
        above = max(theta - ub, 0);
        outOfBounds = below + above;

        if any(outOfBounds > 0)
            range = ub - lb;
            range(~isfinite(range) | range <= 0) = 1;

            val = opts.boundPenaltyScale * ...
                  (1 + sum((outOfBounds ./ range).^2));
            return
        end

        val = objectiveForTheta(theta);
    end

% ============================================================
% Nested function: objective
% ============================================================

    function mse = objectiveForTheta(theta)

        A = unpackA(theta);
        [pred, ~] = predictForA(A);

        resid = pred - traces;
        resid(~validMask) = NaN;

        if opts.normalizeCells
            resid = resid ./ cellScale;
        end

        r = resid(validMask);
        r = r(isfinite(r));

        if isempty(r)
            mse = Inf;
        else
            mse = mean(r.^2);
        end

        % Optional regularization on A.
        if opts.lambdaA > 0 && isfinite(mse)
            dTheta = theta(:) - thetaPrior(:);
            mse = mse + opts.lambdaA * mean(dTheta.^2);
        end
    end

% ============================================================
% Nested function: prediction for a given full A
% ============================================================

    function [pred, info] = predictForA(A)

        pred = zeros(nFrames, nCells);

        % Loop over response/design rows.
        for k = 1:nFrames

            idxAll = nzIdx{k};
            valsAll = nzVal{k};

            if isempty(idxAll)
                continue
            end

            lagIndsAll = ceil(idxAll / nPix);
            pixIndsAll = mod(idxAll - 1, nPix) + 1;

            % Process one lag at a time because each lag can have a different
            % eye position if toeplitz row frame contains stim(frame-rfBins(q)).
            for l = 1:nLags

                useQ = lagIndsAll == l;

                if ~any(useQ)
                    continue
                end

                if opts.useLagSpecificEyePos
                    lagOffset = round(rfBins(l));
                    eyeFrame = k + opts.lagSign * lagOffset;
                else
                    eyeFrame = k;
                end

                if eyeFrame < 1 || eyeFrame > nFrames || ~pupilOK(eyeFrame)
                    continue
                end

                index = pixIndsAll(useQ);
                vals = valsAll(useQ);

                % Coordinates and values of active stimulus pixels/squares.
                xp = xPix(index)';     % [1 x nActive]
                yp = yPix(index)';     % [1 x nActive]
                vals = vals(:);        % [nActive x 1]

                % Full A mapping:
                %
                %   [dxEye; dyEye] = A * [eyeX; eyeY]
                %
                eyeVec = eyePos(eyeFrame, :)';
                dEye = opts.shiftSign * (A * eyeVec);

                dxEye = dEye(1);
                dyEye = dEye(2);

                cx = x0 + dxEye;       % [nCells x 1]
                cy = y0 + dyEye;       % [nCells x 1]

                % Difference between active stimulus locations and shifted RF centers.
                dx = xp - cx;          % [nCells x nActive]
                dy = yp - cy;          % [nCells x nActive]

                % Rotate coordinates into RF's principal axes.
                xr =  ct .* dx + st .* dy;
                yr = -st .* dx + ct .* dy;

                % Evaluate Gaussian RF at active stimulus locations.
                G = exp(-0.5 .* ((xr ./ sx).^2 + (yr ./ sy).^2));

                % Spatial RF value.
                F = amp .* G;

                if opts.useSpatialOffset
                    F = F + rfOffset;
                end

                % Apply temporal weight for this lag, cell by cell.
                F = F .* timeWeights(:, l);

                % Split stimulus into ON and OFF components.
                v_on = vals;
                v_on(v_on < 0) = 0;

                v_off = vals;
                v_off(v_off > 0) = 0;
                v_off = -v_off;

                driveOn = F * v_on;      % [nCells x 1]
                driveOff = F * v_off;    % [nCells x 1]

                ind = bestSubFields == 1; % ON only
                if any(ind)
                    pred(k, ind) = pred(k, ind) + ...
                        (subFieldSigns(ind,1) .* driveOn(ind))';
                end

                ind = bestSubFields == 2; % OFF only
                if any(ind)
                    pred(k, ind) = pred(k, ind) + ...
                        (subFieldSigns(ind,2) .* driveOff(ind))';
                end

                ind = bestSubFields == 3; % ON and OFF
                if any(ind)
                    pred(k, ind) = pred(k, ind) + ...
                        (subFieldSigns(ind,1) .* driveOn(ind) + ...
                         subFieldSigns(ind,2) .* driveOff(ind))';
                end
            end
        end

        info = struct();
        info.cellGain = ones(1, nCells);
        info.cellIntercept = zeros(1, nCells);

        % Optional: analytically fit gain and intercept per cell for this A.
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

% ============================================================
% Helper: pack/unpack A
% ============================================================

function theta = packA(A)
theta = [A(1,1); A(1,2); A(2,1); A(2,2)];
end

function A = unpackA(theta)
theta = theta(:);
A = [theta(1), theta(2); ...
     theta(3), theta(4)];
end

% ============================================================
% Helper: multistart initialization
% ============================================================

function thetaStarts = makeThetaStarts(theta0, lb, ub, opts)

nStarts = max(1, round(opts.nStarts));
thetaStarts = nan(numel(theta0), nStarts);
thetaStarts(:,1) = theta0(:);

if nStarts == 1
    return
end

rng(opts.randomSeed);

range = ub - lb;
range(~isfinite(range) | range <= 0) = 1;

for k = 2:nStarts
    theta = theta0(:) + opts.startScale .* range(:) .* randn(numel(theta0), 1);
    theta = min(max(theta, lb), ub);
    thetaStarts(:,k) = theta;
end
end