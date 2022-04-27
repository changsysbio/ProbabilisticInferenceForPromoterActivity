function statistics = PromoterActivityStatistics(promoterActivity, gapThreshold, normalizedCellCycle)

    %
    % @description: count the statistics from the inferred promoter activity.
    %
    % @parameter gapThreshold: a threshold for merging neighboring ON intervals
    % that are too close to each other.
    %
    % @concept definition
    %
    % An OFF interval is defined as the duration between two bursts. Thus, the
    % interval between the starting timepoint, ending timepoint or a burst would
    % not be counted.
    %
    % Similarly, if an ON interval includes the starting or ending timepoint,
    % it will also not be counted.
    %
    % Since the inference for extremely short OFF intervals are not reliable,
    % which is likely to be a part of an ON intervals, we will merge ON
    % intervals if their gaps are smaller or equal to gapThreshold. An OFF
    % interval, if is smaller or equal to gapThreshold, will not be counted.
    % The burst size of the pulse for a merged ON interval is also merged.
    %

    assert(gapThreshold >= 0);

    statistics.off = [];
    statistics.on = [];
    statistics.amplitude = [];
    statistics.burstSize = [];
    statistics.onIndices = [];
    statistics.offIndices = [];
    
    if exist('normalizedCellCycle')
       statistics.offByCellCycle = [];
    end
    
    if isempty(promoterActivity)
        return;
    end

    [paggt beginIndex endIndex] = PromoterActivityGivenGapThreshold(...
        promoterActivity, gapThreshold);
    if exist('normalizedCellCycle')
        normalizedCellCycle = normalizedCellCycle(beginIndex : endIndex);
    end
    if isempty(paggt)
        return;
    end

    [onIndices, offIndices] = CountIntervals(paggt > 0, gapThreshold);
    statistics.onIndices = [zeros(1, max(1, beginIndex - 1)) onIndices ...
        zeros(1, length(promoterActivity) - endIndex)];
    statistics.offIndices = [zeros(1, max(1, beginIndex - 1)) offIndices ...
        zeros(1, length(promoterActivity) - endIndex)];
    
    for i = 1 : max(offIndices)
    
        index = find(offIndices == i);
        statistics.off(end + 1) = length(index);
        
        if exist('normalizedCellCycle')
            statistics.offByCellCycle(end + 1) = ...
                CalculateOffByCellCycle(normalizedCellCycle(index));
        end
    end

    for i = 1 : max(onIndices)
        index = find(onIndices == i);
        statistics.on(end + 1) = length(index);
        statistics.burstSize(end + 1) = sum(paggt(index));
    end

    statistics.amplitude = paggt(find(paggt));
    
    if exist('normalizedCellCycle')
        statistics.burstFrequency = length(statistics.burstSize) / ...
            CalculateOffByCellCycle(...
                normalizedCellCycle(normalizedCellCycle>=0));
    end
end

function offByCellCycle = CalculateOffByCellCycle(offTimepoints)

    offByCellCycle = 0;
        
    for j = 1 : (length(offTimepoints) - 1)
            
        % For the first/last cell cycle (offTimepoints == -1),
        % we cannot determine its exact cell cycle duration,
        % thus are not included in this analysis
        if offTimepoints(j + 1) < 0
            offByCellCycle = NaN;
            break;
        end
                
        if offTimepoints(j + 1) > offTimepoints(j)
            offByCellCycle = offByCellCycle + ...
                offTimepoints(j + 1) - offTimepoints(j);
        end
    end
end

function [paggt beginIndex endIndex] = PromoterActivityGivenGapThreshold(...
    promoterActivity, gapThreshold)

    % index all ON intervals given gapThreshold
    onIndices = CountIntervals(promoterActivity > 0, gapThreshold);

    % remove the OFF intervals that are too close to the start time point
    beginIndex = find(onIndices, 1);
    endIndex = find(onIndices, 1, 'last');
    if isempty(beginIndex) || isempty(endIndex)
        paggt = [];
        return;
    end

    % remove the first ON interval if it is too close to the start time point
    if beginIndex <= 1 + gapThreshold && onIndices(beginIndex) > 0
        beginIndex = find(onIndices == onIndices(beginIndex), 1, 'last') + 1;
    end

    % remove the last ON interval if it is too close to the end time point
    if endIndex >= length(onIndices) - gapThreshold && onIndices(endIndex) > 0
        endIndex = find(onIndices == onIndices(endIndex), 1) - 1;
    end

    paggt = promoterActivity(beginIndex : endIndex);

end

function PromoterActivityStatisticsTest

    promoterActivity = [];
    gapThreshold = 0;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(isempty(statistics.off));
    assert(isempty(statistics.on));
    assert(isempty(statistics.amplitude));
    assert(isempty(statistics.burstSize));

    promoterActivity = [0 0 0 0 0 0];
    gapThreshold = 0;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(isempty(statistics.off));
    assert(isempty(statistics.on));
    assert(isempty(statistics.amplitude));
    assert(isempty(statistics.burstSize));

    promoterActivity = [0 1 2 0 0 1 3 0.1 0 0 0 1 1];
    gapThreshold = 0;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(~ isempty(statistics.off));
    assert(sum(statistics.off ~= [2 3]) == 0);
    assert(~ isempty(statistics.on));
    assert(sum(statistics.on ~= [2 3]) == 0);
    assert(~ isempty(statistics.amplitude));
    assert(sum(statistics.amplitude ~= [1 2 1 3 0.1000]) == 0);
    assert(~ isempty(statistics.burstSize));
    assert(sum(statistics.burstSize ~= [3 4.1000]) == 0);

    promoterActivity = [0 1 2 0 0 1 3 0.1 0 0 0 1 1];
    gapThreshold = 1;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(~ isempty(statistics.off));
    assert(sum(statistics.off ~= [2 3]) == 0);
    assert(~ isempty(statistics.on));
    assert(sum(statistics.on ~= [3]) == 0);
    assert(~ isempty(statistics.amplitude));
    assert(sum(statistics.amplitude ~= [1 3 0.1000]) == 0);
    assert(~ isempty(statistics.burstSize));
    assert(sum(statistics.burstSize ~= [4.1000]) == 0);

    promoterActivity = [0 1 2 0 0 1 3 0.1 0 0 0 1 1];
    gapThreshold = 2;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(~ isempty(statistics.off));
    assert(sum(statistics.off ~= [3]) == 0);
    assert(isempty(statistics.on));
    assert(isempty(statistics.amplitude));
    assert(isempty(statistics.burstSize));

    promoterActivity = [0 0 0 1 2 0 0 1 3 0.1 0 0 0 1 1];
    gapThreshold = 2;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(~ isempty(statistics.off));
    assert(sum(statistics.off ~= [3]) == 0);
    assert(~ isempty(statistics.on));
    assert(sum(statistics.on ~= [7]) == 0);
    assert(~ isempty(statistics.amplitude));
    assert(sum(statistics.amplitude ~= [1 2 1 3 0.1000]) == 0);
    assert(~ isempty(statistics.burstSize));
    assert(sum(statistics.burstSize ~= [7.1000]) == 0);

    promoterActivity = [0 0 0 1 2 0 0 0 1 3 0.1 0 0 0 1 1];
    gapThreshold = 2;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(~ isempty(statistics.off));
    assert(sum(statistics.off ~= [3 3]) == 0);
    assert(~ isempty(statistics.on));
    assert(sum(statistics.on ~= [2 3]) == 0);
    assert(~ isempty(statistics.amplitude));
    assert(sum(statistics.amplitude ~= [1 2 1 3 0.1000]) == 0);
    assert(~ isempty(statistics.burstSize));
    assert(sum(statistics.burstSize ~= [3 4.1000]) == 0);

    promoterActivity = [0 0 0 0 1 2 0 0 0 1 3 0.1 0 0 0 0 1 1];
    gapThreshold = 3;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(~ isempty(statistics.off));
    assert(sum(statistics.off ~= [4]) == 0);
    assert(~ isempty(statistics.on));
    assert(sum(statistics.on ~= [8]) == 0);
    assert(~ isempty(statistics.amplitude));
    assert(sum(statistics.amplitude ~= [1 2 1 3 0.1000]) == 0);
    assert(~ isempty(statistics.burstSize));
    assert(sum(statistics.burstSize ~= [7.1000]) == 0);

    promoterActivity = [0 0 0 0 1 2 0 0 0 1 3 0.1 0 0 0 1 1];
    gapThreshold = 3;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(isempty(statistics.off));
    assert(isempty(statistics.on));
    assert(isempty(statistics.amplitude));
    assert(isempty(statistics.burstSize));

    promoterActivity = [0 0 0 0 1 2 0 0 0 0 1 3 0.1 0 0 0 1 1];
    gapThreshold = 3;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(~ isempty(statistics.off));
    assert(sum(statistics.off ~= [4]) == 0);
    assert(~ isempty(statistics.on));
    assert(sum(statistics.on ~= [2]) == 0);
    assert(~ isempty(statistics.amplitude));
    assert(sum(statistics.amplitude ~= [1 2]) == 0);
    assert(~ isempty(statistics.burstSize));
    assert(sum(statistics.burstSize ~= [3]) == 0);

    promoterActivity = [28.7917768    0    802.9886628    0    21.71285108    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    53.9615645    691.688902    20.54450301    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    35.60330754    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    34.91562202    555.3314841    0    0    317.0276155    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0];
    gapThreshold = 3;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(~ isempty(statistics.off));
    assert(sum(statistics.off ~= [92 72 108]) == 0);
    assert(~ isempty(statistics.on));
    assert(sum(statistics.on ~= [3 1 5]) == 0);
    assert(~ isempty(statistics.amplitude));
    assert(sum(statistics.amplitude ~= [53.9615645000000 691.688902000000 20.5445030100000 35.6033075400000 34.9156220200000 555.331484100000 317.027615500000]) == 0);
    assert(~ isempty(statistics.burstSize));
    assert(sum(statistics.burstSize ~= [766.194969510000    35.6033075400000    907.274721620000]) == 0);

    promoterActivity = [28.7917768    0    802.9886628    0    21.71285108    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    53.9615645    691.688902    20.54450301    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    35.60330754    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    34.91562202    555.3314841    0    0    317.0276155    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    1    0    0    0];
    gapThreshold = 3;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(~ isempty(statistics.off));
    assert(sum(statistics.off ~= [92 72 108 196]) == 0);
    assert(~ isempty(statistics.on));
    assert(sum(statistics.on ~= [3 1 5]) == 0);
    assert(~ isempty(statistics.amplitude));
    assert(sum(statistics.amplitude ~= [53.9615645000000 691.688902000000 20.5445030100000 35.6033075400000 34.9156220200000 555.331484100000 317.027615500000]) == 0);
    assert(~ isempty(statistics.burstSize));
    assert(sum(statistics.burstSize ~= [766.194969510000    35.6033075400000    907.274721620000]) == 0);

    promoterActivity = [28.7917768    0    802.9886628    0    21.71285108    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    53.9615645    691.688902    20.54450301    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    35.60330754    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    34.91562202    555.3314841    0    0    317.0276155    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    1    0    0    0];
    gapThreshold = 0;
    statistics = PromoterActivityStatistics(promoterActivity, gapThreshold);
    assert(~ isempty(statistics.off));
    assert(sum(statistics.off ~= [1 1 92 72 108 2 196 1]) == 0);
    assert(~ isempty(statistics.on));
    assert(sum(statistics.on ~= [1 1 3 1 2 1 1 1]) == 0);
    assert(~ isempty(statistics.amplitude));
    assert(sum(statistics.amplitude ~= [802.9886628 21.71285108 53.9615645 691.688902 20.54450301 35.60330754 34.91562202 555.3314841 317.0276155 1 1]) == 0);
    assert(~ isempty(statistics.burstSize));
    assert(sum(statistics.burstSize ~= [802.9886628 21.71285108 766.194969510000 35.60330754 590.24710612 317.0276155 1 1]) == 0);

end
