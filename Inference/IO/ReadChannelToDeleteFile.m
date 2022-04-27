function channelsToDelete = ReadChannelToDeleteFile(channelToDeleteFilename)

    %
    % @description: read info about which channels are not include in analysis
    %
    
    channelToDeleteFile = fopen(channelToDeleteFilename, 'r');
    channelsToDelete = fscanf(channelToDeleteFile, '%d');
    channelsToDelete = channelsToDelete';

end
