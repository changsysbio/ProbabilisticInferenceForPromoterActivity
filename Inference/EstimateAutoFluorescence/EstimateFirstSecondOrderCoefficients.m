function [firstOrderCoefficients secondOrderCoefficients] = ...
    EstimateFirstSecondOrderCoefficients(x, y)

    %
    % @description: estimate the first order (mean) and the second order
    % (variance) dependence of y on x.
    %

    x = Preprocessing(x);
    y = Preprocessing(y);

    mdl = fitlm(x, y);
    firstOrderCoefficients = mdl.Coefficients.Estimate;
    predictedY = FirstOrderStatistics(firstOrderCoefficients, x);

    b = (y - predictedY).^2;
    A = [x.^2, x, ones(size(x))];
    secondOrderCoefficients = pinv(A) * b;

end

function x = Preprocessing(x)
    if size(x, 1) == 1
        x = x';
    else
        assert(size(x, 2) == 1);
    end
end
