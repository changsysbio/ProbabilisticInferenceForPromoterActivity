function [background nonnegativeDerivative cellSizeWithNonnegativeDerivative] = ...
    FluorescenceNonnegativeDerivative(info)

    properties = {'cellSize', 'fluorescence_derivative'};
    attribute = CombineAttribute(info, properties);

    nonnegativeIndex = find(attribute.fluorescence_derivative >= 0);
    nonnegativeDerivative = attribute.fluorescence_derivative(nonnegativeIndex);
    cellSizeWithNonnegativeDerivative = attribute.cellSize(nonnegativeIndex);

    background.meanDerivative = mean(nonnegativeDerivative);
    background.stdDerivative = std(nonnegativeDerivative);

    [background.sizeDerivativeFirstOrderCoefficients ...
    background.sizeDerivativeSecondOrderCoefficients] = ...
        EstimateFirstSecondOrderCoefficients(...
            attribute.cellSize(nonnegativeIndex), nonnegativeDerivative);
end
