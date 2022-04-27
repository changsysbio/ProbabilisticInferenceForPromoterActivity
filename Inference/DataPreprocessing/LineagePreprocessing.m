function lineage = LineagePreprocessing(lineage, methodParameter)

    %
    % @description: Data preprocessing: determine cell division, apply signal
    % filter, etc.
    %

    if ~ isfield(lineage, 'division') || isempty(lineage.division)
        lineage.division = DetermineDivision(lineage, methodParameter.partition);
    end

    if strcmp(methodParameter.filter, 'filter')
    
        lineage.fluorescenceBeforeSignalFilter = lineage.fluorescence;
        lineage.fluorescence = SignalFilter(lineage.division, ...
            lineage.fluorescence);
            
        if isfield(lineage, 'realLeakage')
            lineage.realLeakageBeforeSignalFilter = lineage.realLeakage;
            lineage.realLeakage = SignalFilter(lineage.division, ...
                lineage.realLeakage);
        end
        
        lineage.fluorescence_derivative = FluorescenceDerivative(...
            lineage.fluorescence, 1.0);
            
    elseif strcmp(methodParameter.filter, 'no-filter')
        lineage.fluorescence_derivative = FluorescenceDerivative(...
            lineage.fluorescence, 1.0);
        true;
    else
        assert(0);
    end

    lineage.cellCycle = CellCycleInformation(lineage.division, lineage.timeInterval);
end
