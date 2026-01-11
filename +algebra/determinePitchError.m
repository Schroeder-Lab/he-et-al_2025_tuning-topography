function [error, errDir, errOri] = determinePitchError(pitchAngle, ...
    dirPrefs, oriPrefs, rfPos)
%DETERMINEPITCHERROR   Determine difference between measured and predicted 
%direction and orientation preferences, given a specific pitch of the 
%animal's head.

% INPUTS
% pitchAngle    double, pitch of head (zero if lambda and bregma at same
%               height)
% dirPrefs      [units], preferred directions
% oriPrefs      [units], preferred orientations
% rfPos         [units, 2], azimuth and elevation of RF locations

% OUTPUTS
% error         double, MSE across all units for direction and orientation
%               preferences
% errDir        [units], differences between measured and predicted
%               direction preference
% errOri        [units], differences between measured and predicted
%               orientation preference

[ds_trans, os_long, os_lat] = algebra.getDsOsAxes(pitchAngle);

% direction
valid_dir = find(~isnan(dirPrefs) & ~any(isnan(rfPos),2));
transVectors = NaN(length(valid_dir), 4); % 4 direction vectors of DSGCs at RF position
for k = 1:length(valid_dir)
    for v = 1:4
        [~,~, predDir] = algebra.getTranslationDir(ds_trans(v,:), ...
            rfPos(valid_dir(k),:));
        transVectors(k,v) = predDir;
    end
end
err = abs(dirPrefs(valid_dir) - transVectors);
ind = err > 180;
err(ind) = 360 - err(ind);
errDir = min(err, [], 2);

% orientation
% tranform angles from direction of motion (axial motion) to
% orientation of grating
oriPrefs = mod(oriPrefs + 90, 180);
valid_ori = find(~isnan(oriPrefs) & ~any(isnan(rfPos),2));
longVectors = NaN(length(valid_ori), 2); % 2 longitudinal vectors of OSGCs at RF position
latVectors = NaN(length(valid_ori), 2); % 2 latitudinal vectors of OSGCs at RF position
for k = 1:length(valid_ori)
    for v = 1:2
        [~,~, predOri] = algebra.getTranslationDir(os_long(v,:), ...
            rfPos(valid_ori(k),:));
        longVectors(k,v) = mod(predOri, 180);
        [~,~, predOri] = algebra.getLatitudeOrientation(os_lat(v,:), ...
            rfPos(valid_ori(k),:));
        latVectors(k,v) = predOri;
    end
end
err = abs(oriPrefs(valid_ori) - [longVectors, latVectors]);
ind = err > 90;
err(ind) = 180 - err(ind);
errOri = min(err, [], 2);

% total error (MSE)
error = (sum(errDir .^ 2) + sum(errOri .^ 2)) / ...
    (length(errDir) + length(errOri));