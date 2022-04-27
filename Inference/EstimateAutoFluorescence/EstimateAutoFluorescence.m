function autoFluorescence = EstimateAutoFluorescence(controlData)

    %
    % @description: estimate coefficients for auto-fluorescence.
    %

    autoFluorescence = FluorescenceNonnegativeDerivative(controlData);

    properties = {'cellSize', 'fluorescence'};
    controlAttribute = CombineAttribute(controlData, properties);

    [autoFluorescence.sizeFluorescenceCoefficients ...
        autoFluorescence.sizeDependentVarianceCoefficients] ...
        = EstimateFirstSecondOrderCoefficients(controlAttribute.cellSize, ...
        controlAttribute.fluorescence);

end
