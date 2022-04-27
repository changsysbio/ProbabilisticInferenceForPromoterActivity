function statistics = FirstOrderStatistics(coefficients, independentVariable)

    %
    % @description: calculate the mean based on given coefficients and
    % independent variable.
    %

   statistics = independentVariable * coefficients(2) + coefficients(1);
end
