function angles = continuousAngles(angles, type)

switch type
    case 'dir'
        limit = 360;
    case 'ori'
        limit = 180;
end

% m = median(angles);
% ind = angles - m > limit/2;
% angles(ind) = angles(ind) - limit;
% ind = angles - m < -limit/2;
% angles(ind) = angles(ind) + limit;

while any(abs(diff(angles)) > limit/2)
    for k = 1:length(angles)-1
        d = diff(angles(k:k+1));
        if d > limit/2
            angles(k+1) = angles(k+1) - limit;
        elseif d < -limit/2
            angles(k+1) = angles(k+1) + limit;
        end
    end
end

m = mean(angles);
angles = angles - sign(m) * floor(abs(m) / limit) * limit;