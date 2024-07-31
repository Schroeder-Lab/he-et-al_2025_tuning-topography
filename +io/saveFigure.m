function saveFigure(fig, folder, filename)

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
fig.Renderer = 'painters';
saveas(fig, fullfile(folder, filename), 'epsc')

close(fig)