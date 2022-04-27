function PlotCellSizeAndDivision(lineage)
    plot(lineage.timepoints, lineage.cellSize, 'k-*', 'LineWidth', 2);
    hold on;
    for dividingTime = find(lineage.division)
        plot([lineage.timepoints(dividingTime) lineage.timepoints(dividingTime)], ...
            [0 max(lineage.cellSize) * 1.2], ...
            'Color', [0.85 0.85 0.85], 'LineStyle', '--');
    end
    set(gcf, 'Position',  [100, 100, 1500, 800]);
end
