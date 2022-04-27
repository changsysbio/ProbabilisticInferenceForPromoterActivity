function PlotRealFluorescenceAndLeakage(lineage)

    figure;
    subplot(2, 1, 1);
    hold on;
    plot(lineage.timepoints, lineage.fluorescence, 'r.');
    if isfield(lineage, 'inferredRealFluorescence')
        plot(lineage.timepoints, lineage.inferredRealFluorescence, 'b.', ...
            lineage.timepoints, lineage.inferredExpectedFluorescence, 'g.');
    end
    for dividingTime = find(lineage.division)
        plot([lineage.timepoints(dividingTime) lineage.timepoints(dividingTime)], ...
            [min(lineage.fluorescence) - 50 max(lineage.fluorescence) + 50], ...
            'Color', [0.85 0.85 0.85], 'LineStyle', '--');
    end
    title('fluorescence (red) and inferred fluorescence (blue)');
    ylabel('fluorescence (a.u.)');
    subplot(2, 1, 2);
    plot(lineage.timepoints, lineage.inferredLeakage, 'k-');
    index1 = lineage.inferredLeakage > 0;
    hold on;
    plot(lineage.timepoints(index1), lineage.inferredLeakage(index1), 'k+');

    if isfield(lineage, 'realLeakage')
        hold on;
        index2 = find(lineage.realLeakage > 0);
        plot(lineage.timepoints(index2), lineage.realLeakage(index2), 'r+');
    end

    title(['inferred leakage'], 'Interpreter', 'none');
    hold on;
    count = 1;
    index = find(lineage.division);
    for dividingTime = index'
        plot([lineage.timepoints(dividingTime) lineage.timepoints(dividingTime)], ...
            [min(lineage.fluorescence) - 50 max(lineage.fluorescence) + 50], ...
            'Color', [0.85 0.85 0.85], 'LineStyle', '--');
        text(lineage.timepoints(dividingTime), 300, num2str(count));
        count = count + 1;
    end
    ylabel('fluorescence (a.u.)');
    xlabel('time (hr)');
end
