function PlotFluorescenceAndLeakageInOnePanel(lineage)

    box on;
    
    set(gcf, 'Position',  [100, 100, 900, 250]);
    set(gca, 'Fontname', 'Helvetica', 'FontSize', 20);

    ax1 = gca; % current axes

    hold on;
    for dividingTime = find(lineage.division)
        plot([lineage.timepoints(dividingTime) lineage.timepoints(dividingTime)], ...
            [max(lineage.fluorescence) * [0.85 1]], ...
            'Color', [0.85 0.85 0.85], 'LineStyle', '-', 'LineWidth', 2);
    end

    yyaxis left;
    hold on;
    ax1.XColor = 'r';
    ax1.YColor = 'r';
    plot(lineage.timepoints, lineage.fluorescence, 'r-', 'LineWidth', 3);
    if isfield(lineage, 'backgroundFluorescence')
        plot(lineage.timepoints, lineage.inferredBackgroundFluorescence, 'g.');
    end
    yLimPosition = ylim;
    ylabel('fluorescence (a.u.)');
    xLimit = xlim;
    xlim([min(0, xLimit(1)), xLimit(2)]);
    ylim([0 max(lineage.fluorescence)]);
    
    ax1 = gca; % current axes
    yyaxis right;
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    plot(lineage.timepoints, lineage.inferredLeakage, 'k-', 'LineWidth', 2);
    ylim(yLimPosition);
    ylabel('promoter activity (a.u.)');
    xlim([min(0, xLimit(1)), xLimit(2)]);
    ylim([0 max(lineage.fluorescence)]);
    %set(gca, 'YScale', 'log');

    title('fluorescence and promoter activity');
    xlabel('time (hr)');

end
