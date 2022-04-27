function division = DetermineDivision(lineage, partitionMethod)

    %
    % @description: determine the timepoints of cell division and estimate the
    % ratio of the fluorescence signal that does not go to the next timepoint.
    %

    neighboringLengthThreshold = 0.8;
    lengthThreshold = 2.2;

    division = zeros(length(lineage.length), 1);
    for i = 1 : (length(lineage.length) - 1)

        neighboringCheck = (lineage.length(i + 1) ...
            < neighboringLengthThreshold * lineage.length(i));
        lengthCheck = (lineage.length(i) > lengthThreshold);

        if (neighboringCheck & lengthCheck)

            if strcmp (partitionMethod, 'partition-by-fluorescence')
                division(i) = 1.0 - ...
                    lineage.fluorescence(i + 1) / lineage.fluorescence(i);
            elseif strcmp (partitionMethod, 'partition-by-cell-size')
                division(i) = 1.0 - lineage.cellSize(i + 1) / lineage.cellSize(i);
            elseif strcmp (partitionMethod, 'equal-partition') || ...
                strcmp (partitionMethod, 'none')
                division(i) = 0.5;
            else
                assert(0);
            end

            division(i) = max(0.30, division(i));
            division(i) = min(0.70, division(i));

        else
            division(i) = 0;
        end

        assert(division(i) <= 1.0 && division(i) >= 0.0);
    end
    division(end) = 0;

end
