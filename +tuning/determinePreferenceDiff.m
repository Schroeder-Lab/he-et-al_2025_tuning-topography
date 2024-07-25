function differences = determinePreferenceDiff(angles, type)

differences = [];
diff = abs(angles - angles');
diff = diff(tril(true(size(diff)),-1));
if strcmp(type, 'dir')
    ind = diff > 180;
    diff(ind) = 360 - diff(ind);
elseif  strcmp(type, 'ori')
    ind = diff > 90;
    diff(ind) = 180 - diff(ind);
else
    return
end
differences = diff;