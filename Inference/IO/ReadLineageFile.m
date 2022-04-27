function data = ReadLineageFile(data, inputFilename, channelsToDelete, ...
    methodParameter, prefix)

    %
    % @description: read lineage data.
    %

    rawData = tdfread(inputFilename);

    channels = setdiff(unique(rawData.channel_in_multipoint'), channelsToDelete);
    for channel = channels
    
        indices = find(rawData.channel_in_multipoint == channel);
        
        data.lineageSpecific(end + 1).timepoints = rawData.cell_age(indices);
        TimepointCheck(data.lineageSpecific(end).timepoints);
        
        data.lineageSpecific(end).fluorescence = ...
            rawData.fluorescence_total_0(indices);
        data.lineageSpecific(end).length = rawData.length(indices);
        data.lineageSpecific(end).channel_width = rawData.channel_width(indices);
        data.lineageSpecific(end).cellSize = ...
            data.lineageSpecific(end).length .* ...
            data.lineageSpecific(end).channel_width;
            
        data.lineageSpecific(end).prefix = [prefix '_' num2str(channel)];
        data.lineageSpecific(end).timeInterval = methodParameter.timeInterval;
        data.lineageSpecific(end).gapThreshold = methodParameter.gapThreshold;
  
        data.lineageSpecific(end).division = [];
        data.lineageSpecific(end).fluorescenceBeforeSignalFilter = [];
        data.lineageSpecific(end).fluorescence_derivative = [];
        data.lineageSpecific(end).cellCycle = [];
    end

end

function TimepointCheck(timepoints)

    % timepoints increase monotonically
    assert(sum(diff(timepoints) <= 0) == 0);
    
    % and it inreases with a same dwell time
    assert(length(unique(round(diff(timepoints), 3))) == 1);
end
