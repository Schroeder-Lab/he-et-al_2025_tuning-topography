function [colormapON, colormapOFF] = getRFMaps()
%GETRFMAPS   Return color maps used for plotting receptive fields (ON and
%OFF separately).

% OUTPUTS
% colormapON    [100 x 3], grey-white-pink color gradient
% colormapOFF   [100 x 3], grey-white-blue color gradient

red = [1 0 .5];
blue = [0 .5 1];
black = [1 1 1].*0.5;

grad = linspace(0,1,100)';

reds = red.*flip(grad) + [1 1 1].*grad;
blacks = black.*flip(grad) + [1 1 1].*grad;
colormapON = [blacks; flip(reds(1:end-1,:),1)];

blues = blue.*flip(grad) + [1 1 1].*grad;
colormapOFF = [blacks; flip(blues(1:end-1,:),1)];