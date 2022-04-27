function InferenceWithInput(experimentStrain, method, autoFluorescence)

    %
    % @description: Idential to @function Inference, except @parameter
    % autoFluorescence is input by the user.
    %
    % @parameter autoFluorescence: coefficients of auto-fluorescence estimated
    % from a background strain without any fluorescence reporter.
    %

    methodParameter = MethodParameter(experimentStrain, method);

    experimentData = LineageSpecificData(methodParameter, ...
        methodParameter.experimentStrain);

    experimentInformation = AnalysisExperiment(experimentData.lineageSpecific, ...
        methodParameter, autoFluorescence);

    save([methodParameter.experimentStrain '.mat']);
end

