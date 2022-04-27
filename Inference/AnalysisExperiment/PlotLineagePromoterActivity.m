function PlotLineagePromoterActivity(lineage, methodParameter, prefix)

    %
    % @description: plot the promoter activity by lineage.
    %

    if strcmp(methodParameter.method, 'probabilistic')
        for (i = 1 : methodParameter.repeat)
            PlotRealFluorescenceAndLeakage(lineage);
            inferenceMovie(i) = getframe(gcf);
            close;
        end

        v = VideoWriter([prefix '_InferenceMovie.avi'], 'Motion JPEG AVI');
        v.Quality = 95;
        open(v);
        writeVideo(v, inferenceMovie);
        close(v);

        PlotRealFluorescenceAndLeakage(lineage);

        saveas(gcf, [prefix '_ExperimentAnalysis.eps'], 'epsc');
        saveas(gcf, [prefix '_ExperimentAnalysis.fig'], 'fig');
        close;
    end
    
    PlotFluorescenceAndLeakageInOnePanel(lineage);

    saveas(gcf, [prefix '_ExperimentAnalysis2.eps'], 'epsc');
    saveas(gcf, [prefix '_ExperimentAnalysis2.fig'], 'fig');
    close;

    PlotCellSizeAndDivision(lineage);
    saveas(gcf, [prefix '_CellSizeAndDivision.eps'], 'epsc');
    saveas(gcf, [prefix '_CellSizeAndDivision.fig'], 'fig');
    close;

end
