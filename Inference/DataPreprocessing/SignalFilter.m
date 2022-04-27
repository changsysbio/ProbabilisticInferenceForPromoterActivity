function fluorescenceAfterSignalFilter = SignalFilter(division, fluorescence)

    %
    % @description: filter applied on the total fluorescence per cell data
    % for data smoothing.
    %

    dividingTimepointIndices = find(division)';
    dividingIntervals = [1 dividingTimepointIndices length(fluorescence)];

    for i = 1 : length(dividingIntervals) - 1

        if i == 1
            left = 1;
        else
            left = dividingIntervals(i) + 1;
        end

        right = dividingIntervals(i + 1);

        intervalIndices = left:right;

        if length(fluorescence(intervalIndices)) >= 3
            fluorescenceAfterSignalFilter(intervalIndices) = ...
                sgolayfilt(fluorescence(intervalIndices), 1, 3);
        else
            fluorescenceAfterSignalFilter(intervalIndices) = ...
                fluorescence(intervalIndices);
        end

    end

    fluorescenceAfterSignalFilter = fluorescenceAfterSignalFilter';
    assert(length(fluorescence) == length(fluorescenceAfterSignalFilter));

end
