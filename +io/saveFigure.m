function saveFigure(fig, folder, filename)
%SAVEFIGURE   Save Matlab figure in matlab, jpeg, and epsc format, then
%close figure.

% INPUTS
% fig       figure handle
% folder    path to folder to save figure
% filename  filename for saving, without suffix

fMatlab = fullfile(folder, 'matlab');
if ~isfolder(fMatlab)
    mkdir(fMatlab)
end
fJpeg = fullfile(folder, 'jpeg');
if ~isfolder(fJpeg)
    mkdir(fJpeg)
end
savefig(fig, fullfile(fMatlab, filename), 'compact')
saveas(fig, fullfile(fJpeg, filename), 'jpeg')
% fig.Renderer = 'painters';
% saveas(fig, fullfile(folder, filename), 'epsc')
exportgraphics(fig, fullfile(folder, [filename '.eps']), "ContentType", "vector")
% saveas(fig, fullfile(folder, filename), 'pdf')

close(fig)