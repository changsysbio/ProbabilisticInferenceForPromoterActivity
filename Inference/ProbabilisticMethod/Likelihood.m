classdef Likelihood

    %
    % @description: The probabilistic algorithm.
    %

    properties
        leakage;
        realFluorescence;
        likelihood;
    end

    properties
        repeat;
        backgroundFluorescence;
        backgroundVariance;
        backgroundDerivative;
        backgroundDerivativeVariance;
        cellSize;
        fluorescence;
        division;
        background;
        sigma;
    end

    % public functions
    methods

        % constructor
        function obj = Likelihood(lineage, background, methodParameter)

            obj.leakage = zeros(length(lineage.timepoints), methodParameter.repeat);
            obj.realFluorescence = zeros(length(lineage.timepoints), ...
                methodParameter.repeat);
            obj.likelihood = zeros(1, methodParameter.repeat);

            obj.repeat = methodParameter.repeat;
            obj.fluorescence = lineage.fluorescence;
            obj.division = lineage.division;

            obj.background = background;

            obj.backgroundFluorescence = FirstOrderStatistics(...
                background.sizeFluorescenceCoefficients, lineage.cellSize);

            obj.backgroundVariance = SecondOrderStatistics(...
                background.sizeDependentVarianceCoefficients, lineage.cellSize);

            obj.backgroundDerivative = FirstOrderStatistics(...
                background.sizeDerivativeFirstOrderCoefficients, lineage.cellSize);

            obj.backgroundDerivativeVariance = SecondOrderStatistics(...
                background.sizeDerivativeSecondOrderCoefficients, lineage.cellSize);

            obj.cellSize = lineage.cellSize;
            obj.sigma = methodParameter.sigma;

        end

        function likelihood = UpdateParameters(likelihood)

            % initialize leakage, allowing negative value
            likelihood = UpdateParametersWithThreshold(likelihood, -100.0);

            % using previous optimal result as starting point to find optimal
            % result with constraints leakage >= 0
            likelihood = UpdateParametersWithThreshold(likelihood, 0.0);
            
            % Retaining the same burst size for an ON interval, but leakage in
            % the peak would be 'diffused' to nearby timepoints without leakage
            % if likelihood could be improved.
            likelihood = UpdateParametersByDiffusion(likelihood);
        end

        function obj = UpdateParametersWithThreshold(obj, threshold)

            for (likelihoodIndex = 1 : obj.repeat)
                for (j = 1:5)
                    for ( timeIndex = [1 (randperm(length(obj.fluorescence) - 1) + 1)] )

                        obj.leakage(:, likelihoodIndex) = ...
                            Likelihood.UpdateLeakage(obj, likelihoodIndex, ...
                                timeIndex, obj.sigma, threshold);
                        obj.realFluorescence(:, likelihoodIndex) = ...
                            Likelihood.UpdateRealFluorescence(obj, likelihoodIndex);

                    end
                end
                obj.likelihood(likelihoodIndex) = Likelihood.SquareError(obj, ...
                    likelihoodIndex);
            end
        end
        
        function obj = UpdateParametersByDiffusion(obj)
        
            for (likelihoodIndex = 1 : obj.repeat)
                
                onIndices = CountIntervals(...
                    obj.leakage(:, likelihoodIndex) > 0, 0);
                
                for burstIndex = 1 : max(onIndices)
                    obj = Likelihood.DiffusionUpdate(obj, likelihoodIndex, ...
                        onIndices, burstIndex, 'left');
                    obj = Likelihood.DiffusionUpdate(obj, likelihoodIndex, ...
                        onIndices, burstIndex, 'right');
                end
            end
        end
        
    end

    % private functions
    methods(Static)

        function obj = DiffusionUpdate(obj, likelihoodIndex, onIndices, ...
            burstIndex, choice)
        
            assert(strcmp(choice, 'left') || strcmp(choice, 'right'));
            
            if strcmp(choice, 'left')
                pulseIndex = find(onIndices == burstIndex, 1);
            elseif strcmp(choice, 'right')
                pulseIndex = find(onIndices == burstIndex, 1, 'last');
            end

            likelihood = Likelihood.SquareError(obj, likelihoodIndex);
            while true
                        
                deltaSignal = obj.fluorescence - obj.backgroundFluorescence
                    - obj.realFluorescence(:, likelihoodIndex);
                    
                if strcmp(choice, 'left')
                    [updatedLeakage pulseIndex] = Likelihood.DiffusionToLeft(...
                        obj.leakage(:, likelihoodIndex), deltaSignal, ...
                        pulseIndex, obj.background.meanDerivative + ...
                        obj.background.stdDerivative, obj.division);
                elseif strcmp(choice, 'right')
                    [updatedLeakage pulseIndex] = Likelihood.DiffusionToRight(...
                        obj.leakage(:, likelihoodIndex), deltaSignal, ...
                        pulseIndex, obj.background.meanDerivative + ...
                        obj.background.stdDerivative, obj.division);
                end
                
                updatedObj = obj;
                updatedObj.leakage(:, likelihoodIndex) = updatedLeakage;
                updatedObj.realFluorescence(:, likelihoodIndex) = ...
                    Likelihood.UpdateRealFluorescence(updatedObj, ...
                    likelihoodIndex);
                updatedLikelihood = Likelihood.SquareError(updatedObj, ...
                    likelihoodIndex);
                        
                if updatedLikelihood > likelihood
                    obj = updatedObj;
                    likelihood = updatedLikelihood;
                    
                    clearvars updatedObj, updatedLeakage;
                else
                    break;
                end
                        
            end
        end
        
        function [leakage index] = DiffusionToLeft(leakage, delta, index, threshold, division)
        
            if index == 1
                return;
            end
            
            if leakage(index - 1) > 0
                return;
            end
            
            diffusedLeakage = max(0.0, min(leakage(index), delta(index - 1)));
            if diffusedLeakage < threshold * (1 - division(index - 1))
                return;
            end
            
            leakage(index - 1) = diffusedLeakage;
            leakage(index) = leakage(index) - diffusedLeakage;
            index = index - 1;
        end

        function [leakage index] = DiffusionToRight(leakage, delta, index, threshold, division)
        
            if index == length(leakage)
                return;
            end
            
            if leakage(index + 1) > 0
                return;
            end
            
            diffusedLeakage = max(0.0, min(leakage(index), delta(index)));
            if diffusedLeakage < threshold * (1 - division(index))
                return;
            end
            
            leakage(index + 1) = diffusedLeakage;
            leakage(index) = leakage(index) - diffusedLeakage;
            index = index + 1;
        end
        
        function likelihood = SquareError(obj, likelihoodIndex)

            expectedFluorescence = obj.realFluorescence(:, likelihoodIndex) + ...
                obj.backgroundFluorescence;
            likelihood = - sum( (obj.fluorescence - expectedFluorescence).^2 ...
                ./ obj.backgroundVariance);

        end

        function realFluorescence = UpdateRealFluorescence(obj, likelihoodIndex)

            realFluorescence = zeros(size(obj.fluorescence));

            realFluorescence(1) = obj.leakage(1, likelihoodIndex);
            for t = 2 : length(realFluorescence)
                realFluorescence(t) = (1 - obj.division(t - 1)) ...
                    * realFluorescence(t - 1) + obj.leakage(t, likelihoodIndex);
            end
        end

        function effectiveLeakage = UpdateEffectiveLeakage(leakage, ...
            division, timeIndex)

            if timeIndex > 1
                effectiveLeakage(1 : timeIndex - 1) = zeros(1, timeIndex - 1);
            end

            effectiveLeakage(timeIndex) = leakage(timeIndex);
            for t = (timeIndex + 1) : length(leakage)
                effectiveLeakage(t) = (1 - division(t - 1)) ...
                    * effectiveLeakage(t - 1) + leakage(t);
            end
        end

        function leakage = UpdateLeakage(obj, likelihoodIndex, timeIndex, ...
            sigma, threshold)

            assert(threshold <= 0.0);

            leakage = obj.leakage(:, likelihoodIndex);
            divisionCoefficient = Likelihood.DivisionCoefficient(...
                obj.division, timeIndex);
            effectiveLeakage = Likelihood.UpdateEffectiveLeakage(...
                leakage, obj.division, timeIndex);

            numerator = 0.0;
            denominator = 0.0;
            for tau = timeIndex : length(obj.fluorescence)
                numerator = numerator + divisionCoefficient(tau) * ...
                    (obj.fluorescence(tau) - effectiveLeakage(tau) - ...
                    obj.backgroundFluorescence(tau)) / ...
                    obj.backgroundVariance(tau);
                denominator = denominator + divisionCoefficient(tau)^2 ...
                    / obj.backgroundVariance(tau);
            end

            if (timeIndex == 1)
                expectedLeakage = max(threshold, numerator / denominator);
            else
                assert(timeIndex > 1);
                expectedLeakage = max(threshold, numerator / denominator - ...
                    obj.realFluorescence(timeIndex - 1, likelihoodIndex) * ...
                    (1 - obj.division(timeIndex - 1)));
            end

            if Likelihood.PassThreshold(sigma, obj, ...
                timeIndex, expectedLeakage, likelihoodIndex)
                leakage(timeIndex) = expectedLeakage;
            elseif threshold >= 0.0
                if strcmp(sigma.nonPassMethod, 'zero')
                    leakage(timeIndex) = 0.0;
                elseif strcmp(sigma.nonPassMethod, 'nonzero')
                    leakage(timeIndex) = max(0.0, leakage(timeIndex));
                else
                    assert(0);
                end
            end

        end

        function status = PassThreshold(sigma, obj, ...
            timeIndex, expectedLeakage, likelihoodIndex)

            if ~ isempty(sigma.signal)

                statusSignalFrameOfReference = ...
                    ((timeIndex == 1 && ...
                    obj.fluorescence(timeIndex) ...
                    > obj.backgroundFluorescence(timeIndex) ...
                    + sigma.signal * sqrt(obj.backgroundVariance(timeIndex))) ...
                    || ...
                    (timeIndex > 1 && ...
                    obj.fluorescence(timeIndex) ...
                    > obj.backgroundFluorescence(timeIndex) ...
                    + obj.realFluorescence(timeIndex - 1, likelihoodIndex) ...
                    + sigma.signal * sqrt(obj.backgroundVariance(timeIndex))));

            end

            if ~ isempty(sigma.derivative)

                if timeIndex > 1

                    derivative = obj.fluorescence(timeIndex) - ...
                        obj.fluorescence(timeIndex - 1) * ...
                        (1 - obj.division(timeIndex - 1));

                    statusFluorescenceConstantDerivative = (...
                        derivative > obj.background.meanDerivative ...
                        + sigma.derivative * obj.background.stdDerivative);
                    statusFluorescenceSizeDependentDerivative = (derivative > ...
                        obj.backgroundDerivative(timeIndex) ...
                        + sigma.derivative * ...
                            sqrt(obj.backgroundDerivativeVariance(timeIndex)));
                else
                    assert(timeIndex == 1);
                    statusFluorescenceConstantDerivative = true;
                    statusFluorescenceSizeDependentDerivative = true;
                end
            end

            if strcmp(sigma.method, 'signal-frame-of-reference')
                status = statusSignalFrameOfReference;
                return;
            elseif strcmp(sigma.method, 'fluorescence-constant-derivative')
                status = statusFluorescenceConstantDerivative;
                return;
            elseif strcmp(sigma.method, 'fluorescence-size-dependent-derivative')
                status = statusFluorescenceSizeDependentDerivative;
                return;
            elseif strcmp(sigma.method, ...
                'signal-frame-of-reference-fluorescence-constant-derivative')
                status = (statusSignalFrameOfReference && ...
                    statusFluorescenceConstantDerivative);
                return;
            elseif strcmp(sigma.method, ...
                'signal-frame-of-reference-fluorescence-size-dependent-derivative')
                status = (statusSignalFrameOfReference && ...
                    statusFluorescenceSizeDependentDerivative);
                return;
            else
                assert(false);
            end
        end

        function divisionCoefficient = DivisionCoefficient(division, timeIndex)

            divisionCoefficient = zeros(size(division));
            for timepoint = 1 : length(division)

                if timepoint < timeIndex
                    divisionCoefficient(timepoint) = 0;
                    continue;
                elseif timepoint == timeIndex
                    % the starting timepoint should be included.
                    divisionCoefficient(timepoint) = 1;
                    continue;
                end
                assert(timepoint > timeIndex);
                divisionCoefficient(timepoint) = ...
                    divisionCoefficient(timepoint - 1) ...
                    * (1 - division(timepoint - 1));
            end
        end
    end
end
