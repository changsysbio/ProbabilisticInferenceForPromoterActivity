function Inference(experimentStrain, method, choice)

    %
    % @description: The entry point of the program.
    %
    % @parameter experimentStrain: the name of the experimental strain.
    %
    % @parameter method: 'probabilistic' or 'threshold', where the former is
    % our probabilistic pulse detection algorithm, and the latter is the
    % 'hard threshold' method.
    %

    methodParameter = MethodParameter(experimentStrain, method);

    controlData = LineageSpecificData(methodParameter, ...
        methodParameter.controlStrain);
    autoFluorescence = EstimateAutoFluorescence(controlData);
    experimentData = LineageSpecificData(methodParameter, ...
        methodParameter.experimentStrain);

    if exist('choice') && strcmp(choice, 'do-not-analyze-control')
        % do nothing
    else
        controlInformation = AnalysisExperiment(...
            controlData.lineageSpecific, methodParameter, autoFluorescence);
    end
    experimentInformation = AnalysisExperiment(experimentData.lineageSpecific, ...
        methodParameter, autoFluorescence);

    save([methodParameter.experimentStrain '.mat']);
end
