function [onIndices, offIndices] = CountIntervals(status, gapThreshold)

    currentStatus = status(1);
    onIndices = zeros(1, length(status));
    offIndices = zeros(1, length(status));

    if status(1) == 0
        currentOffIntervalIndex = 1;
        currentOnIntervalIndex = 0;
    else
        currentOffIntervalIndex = 0;
        currentOnIntervalIndex = 1;
    end
    i = 1;
    while i <= length(status)

        if currentStatus == status(i) && currentStatus == 0
            offIndices(i) = currentOffIntervalIndex;
        elseif currentStatus == status(i) && currentStatus == 1
            onIndices(i) = currentOnIntervalIndex;
        elseif currentStatus ~= status(i) && currentStatus == 0
            currentOnIntervalIndex = currentOnIntervalIndex + 1;
            onIndices(i) = currentOnIntervalIndex;
            currentStatus = ~currentStatus;
        elseif currentStatus ~= status(i) && currentStatus == 1

            hasGap = false;
            gapIndex = [];

            for j = 0 : min(gapThreshold, length(status) - i)
                if status(i + j) == 1
                    hasGap = true;
                    gapIndex = i + j;
                    break;
                end
            end

            if hasGap
                onIndices(i : gapIndex) = currentOnIntervalIndex;
                i = gapIndex;
            else
                currentOffIntervalIndex = currentOffIntervalIndex + 1;
                offIndices(i) = currentOffIntervalIndex;
                currentStatus = ~currentStatus;
            end

        else
            assert(false);
        end

        i = i + 1;

    end
end
