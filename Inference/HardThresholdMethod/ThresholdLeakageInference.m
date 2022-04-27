function lineage = ThresholdLeakageInference(lineage, threshold, ...
    autoFluorescence, timeInterval)

    %
    % @description: The 'hard threshold' method.
    %

    lineage.inferredLeakage = zeros(size(lineage.fluorescence_derivative));
    for i = 1 : length(lineage.fluorescence_derivative)
        if lineage.fluorescence_derivative(i) > (threshold * timeInterval / 60)
            lineage.inferredLeakage(i) = lineage.fluorescence_derivative(i);
        end
    end

    lineage.inferredBackgroundFluorescence = FirstOrderStatistics(...
        autoFluorescence.sizeFluorescenceCoefficients, lineage.cellSize);
end

