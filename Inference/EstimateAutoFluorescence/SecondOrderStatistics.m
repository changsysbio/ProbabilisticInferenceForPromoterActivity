function statistics = SecondOrderStatistics(coefficients, independentVariable)

    %
    % @description: calculate the variance based on given coefficients and
    % independent variable.
    %

    statistics = coefficients(1) * independentVariable.^2 + ...
        coefficients(2) * independentVariable +coefficients(3);
end
