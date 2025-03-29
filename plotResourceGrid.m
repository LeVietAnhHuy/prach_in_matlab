function plotResourceGrid(grid)

    cmap = parula(64);
    axRG = axes;
    image(100*abs(grid));
    axis(axRG,'xy');
    colormap(axRG, cmap);
    title(axRG,sprintf('PRACH Resource Grid (Size [%s])',strjoin(string(size(grid)),' ')));
    xlabel(axRG,'Symbols'); ylabel(axRG,'Subcarriers in RE');
end