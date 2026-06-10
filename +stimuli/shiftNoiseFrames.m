function stim_shifted = shiftNoiseFrames(stim, stimSize, ...
    gridY, gridX, stimShifts)

% find largest shift in x and y in units of pixels
maxi = max(abs(stimShifts), [], 1, "omitnan");
gridWidth = median(diff(gridX));
gridHeight = -median(diff(gridY));
maxShifts = ceil(maxi ./ [gridWidth, gridHeight]);

if all(maxShifts == 0)
    stim_shifted = stim;
    return
end

stim = reshape(stim, [size(stim, 1), stimSize]); % [frames, y, x]

% pad stimulus so shift is always inside matrix
stim = padarray(stim, [0 maxShifts(2), maxShifts(1)], 0, 'both'); % Pad the stimulus

% Initialize the output array
stim_shifted = zeros(size(stim)); 

% pad pixel grid, then make meshgrid
if maxShifts(1) > 0
    gridX = [(-maxShifts(1):-1) .* gridWidth + gridX(1), gridX, ...
        (1:maxShifts(1)) .* gridWidth + gridX(end)];
end
if maxShifts(2) > 0
    gridY = [(maxShifts(2):-1:1) .* gridHeight + gridY(1), gridY, ...
        (-1:-1:-maxShifts(2)) .* gridHeight + gridY(end)];
end
[gridX, gridY] = meshgrid(gridX, gridY);

% shift each frame, use linear interpolation
for i = 1:size(stim, 1)
    stim_shifted(i,:,:) = interp2(gridX, gridY, squeeze(stim(i,:,:)), ...
        gridX + stimShifts(i,1), gridY + stimShifts(i,2));
end

% cut padded edges off
stim_shifted = stim_shifted(:, maxShifts(2)+1:end-maxShifts(2), ...
    maxShifts(1)+1:end-maxShifts(1), :);

stim_shifted = reshape(stim_shifted, size(stim_shifted, 1), []); % [frames, pixels]