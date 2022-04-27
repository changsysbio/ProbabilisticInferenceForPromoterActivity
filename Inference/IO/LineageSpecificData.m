function data = LineageSpecificData(methodParameter, strain)

    %
    % @description: read data and organize by lineage. Data preprocessing is
    % also applied (determine cell division, apply signal filter, etc.). 
    %

    logFilename = [strain '_log.txt'];
    logFile = fopen(logFilename, 'w');

    data.lineageSpecific = [];
    for FOV = methodParameter.FOVs

        FOV = char(FOV);
        subFolder = [FOV '/' FOV '_' strain '/'];
        prefix = [subFolder FOV '_' strain];
        channelToDeleteFilename = [prefix '.ctd'];
        dataFilename = [prefix '.tsv'];

        if ~ exist(subFolder, 'dir') | ...
            ~ exist(channelToDeleteFilename, 'file') | ...
            ~ exist(dataFilename, 'file')
            continue;
        else
            fprintf(logFile, 'Included, %s : \n', subFolder);
        end

        % read channels-to-delete file
        channelsToDelete = ReadChannelToDeleteFile(channelToDeleteFilename);
        fprintf(logFile, 'exclude channels: ');
        fprintf(logFile, ' %d,', channelsToDelete);
        fprintf(logFile, '\n');

        % read lineage data
        data = ReadLineageFile(data, dataFilename, channelsToDelete, ...
            methodParameter, prefix);

    end

    for i = 1 : length(data.lineageSpecific)
        data.lineageSpecific(i) = LineagePreprocessing(...
            data.lineageSpecific(i), methodParameter);
    end
end
