function methodParameter = MethodParameter(experimentStrain, method)

    %
    % @description: It specifies the details of the inference method.
    %

    methodParameter.timeInterval = 5.0;
    methodParameter.controlStrain = 'CL-41alpha';
    methodParameter.experimentStrain = experimentStrain;
    methodParameter.method = method;
    for i = 1 : 16
        FOVs(i) = string(num2str(i, '%03.f'));
    end
    methodParameter.FOVs = FOVs;

    if strcmp(method, 'threshold')

        methodParameter.threshold = 1250;
        methodParameter.filter = 'filter';
        methodParameter.partition = 'none';
        methodParameter.gapThreshold = 0;
        methodParameter.outputFolder = ['threshold=' ...
            num2str(methodParameter.threshold)];

    elseif strcmp(method, 'probabilistic')

        sigma.signal = 2.5;
        sigma.derivative = 1.5;
        sigma.method = 'signal-frame-of-reference-fluorescence-constant-derivative';
        sigma.nonPassMethod = 'zero';

        methodParameter.sigma = sigma;
        methodParameter.filter = 'filter';
        methodParameter.partition = 'partition-by-fluorescence';
        methodParameter.repeat = 10;
        methodParameter.gapThreshold = 3;
        methodParameter.outputFolder = [sigma.method ...
            '---signal=' num2str(sigma.signal) ...
            '_derivative=' num2str(sigma.derivative) ...
            '_partition=' methodParameter.partition ...
            '_repeat=' num2str(methodParameter.repeat)];

    else
        assert(0);
    end

    if ~exist(methodParameter.outputFolder, 'dir')
        mkdir(methodParameter.outputFolder);
    end

    outputFilename = [methodParameter.outputFolder '/method_parameter.mat'];
    if  ~ exist(outputFilename, 'file')
        save(outputFilename, 'methodParameter');
    end

end
