function derivative = FluorescenceDerivative(fluorescence, timeInterval)

    % this ensures the last derivative is zero
    derivative = zeros(size(fluorescence));

    for t = 1 : (length(fluorescence) - 1)
        derivative(t) = (fluorescence(t + 1) - fluorescence(t)) / timeInterval;
    end
end
