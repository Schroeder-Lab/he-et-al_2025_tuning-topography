function plotRF(rf, rfGaussPars, bestSubField, edges, edges_rf, gridW, ...
    gridH)
%PLOTRF   Plot mapped ON and OFF subfields with contour of fitted 2D
%Gaussian.

% INPUTS
% rf            [rows x columns x subfield], mapped RF
% rfGaussPars   [1 x 7], parameters of fitted 2D Gaussian
% bestSubField  integer, 1: 'ON', 2: 'OFF', or 3: 'ON+OFF'
% edges         [left rigth top bottom] of RF map
% edges_rf      [left rigth top bottom] of RF map to be plotted
% gridW         double, sample spacing along azimuth
% gridH         double, sample spacing along elevation

ellipse_x = linspace(-pi, pi, 100);
[cm_ON, cm_OFF] = colmaps.getRFMaps;
cms = cat(3, cm_ON, cm_OFF);
titles = {'ON field','OFF field'};

rf(:,:,2) = -rf(:,:,2);
mx = max(abs(rf),[],"all");

figure('Position', [75 195 1470 475])
for sf = 1:2
    subplot(1,2,sf)
    % STA
    imagesc([edges(1)+gridW/2 edges(2)-gridW/2], ...
        [edges(3)-gridH/2 edges(4)+gridH/2], ...
        rf(:,:,sf),[-mx mx])
    hold on
    if ismember(bestSubField, [sf 3])
        line = '-';
    else
        line = ':';
    end
    % ellipse at 1 STD (x and y), not rotated, not shifted
    x = rfGaussPars(3) * cos(ellipse_x);
    y = rfGaussPars(5) * sin(ellipse_x);
    % rotate and shift ellipse
    x_rot = rfGaussPars(2) + ...
        x .* cos(rfGaussPars(6)) - ...
        y .* sin(rfGaussPars(6));
    y_rot = rfGaussPars(4) + ...
        x .* sin(rfGaussPars(6)) + ...
        y .* cos(rfGaussPars(6));
    n = x_rot < edges(1) | x_rot > edges(2) | ...
        y_rot > edges(3) | y_rot < edges(4);
    x_rot(n) = NaN;
    y_rot(n) = NaN;
    plot(x_rot, y_rot, ['k' line], 'LineWidth', 2)
    axis image
    set(gca, 'box', 'off', 'YDir', 'normal')
    xlim(edges_rf([1 2]))
    ylim(edges_rf([4 3]))
    colormap(gca, cms(:,:,sf))
    title(titles{sf})
    colorbar
end