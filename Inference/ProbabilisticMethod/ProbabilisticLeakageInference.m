function lineage = ProbabilisticLeakageInference(lineage, background, ...
    methodParameter)

    %
    % @description: maximum likelihood solution for promoter activity
    %

    likelihood = Likelihood(lineage, background, methodParameter);
    likelihood = UpdateParameters(likelihood);

    [maxLikelihood, maxLikelihoodIndex] = max(likelihood.likelihood(end, :));

    lineage.likelihood = maxLikelihood;
    lineage.inferredLeakage = likelihood.leakage(:, maxLikelihoodIndex);
    lineage.inferredRealFluorescence = likelihood.realFluorescence(...
        :, maxLikelihoodIndex);
    lineage.inferredBackgroundFluorescence = likelihood.backgroundFluorescence;
    lineage.inferredExpectedFluorescence = lineage.inferredRealFluorescence ...
        + lineage.inferredBackgroundFluorescence;
        
end
