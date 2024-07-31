function colormap = getGCaMPMap()
%GETGCAMPMAP   Return color map used for plotting mean imaging frame.

% OUTPUTS
% colormap    [300 x 3], black-green color gradient

green = [.5 1 .5];

grad = linspace(0,1,300)';

colormap = green .* grad;