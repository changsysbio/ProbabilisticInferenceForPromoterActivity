function information = AnalysisExperiment(lineages, methodParameter, ...
    autoFluorescence)

    %
    % @description: analysis experiment and extract information.
    %

    information = Information();

    for i = 1 : length(lineages)

        if strcmp(methodParameter.method, 'probabilistic')
            lineage = ProbabilisticLeakageInference(lineages(i), ...
                autoFluorescence, methodParameter);
        elseif strcmp(methodParameter.method, 'threshold')
            lineage = ThresholdLeakageInference(lineages(i), ...
                methodParameter.threshold, autoFluorescence, ...
                methodParameter.timeInterval);
        else
            assert(0);
        end

        PlotLineagePromoterActivity(lineage, methodParameter, lineage.prefix);

        lineageInformation = Information(lineage);
        information = Merge(information, lineageInformation);

        save([lineage.prefix '.mat'], 'lineage', 'methodParameter', ...
            'autoFluorescence', 'lineageInformation');
    end
    
    information = OverallStatistics(information);

end
