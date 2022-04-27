function cellCycle = CellCycleInformation(division, timeInterval)

    %
    % @return [cellCycle.startingAnalysis cellCycle.endAnalysis]
    % defines the range for analysis. It only includes complete cell cycles.
    % for those not in complete cell cycles,in cellCycle.unnormalizedCellCycle
    % and cellCycle.normalizedCellCycle, the values are assigned as -1.
    %
    % @return cellCycle.unnormalizedCellCycle
    % For the unnormalized cell cycle index, now it would starts from 0.
    % e.g., 0, 5, 10, ..., 55 min, and total duration is 60 min.
    %
    % @return cellCycle.normalizedCellCycle
    % e.g., 0, 5/55, 10/55, ..., 55/55 (55 = 60 - 5).
    %
    % @test
    %
    % division = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];
    % cellCycle = CellCycleInformation(division, 5);
    %    startingAnalysis: 5
    %    endAnalysis: 16
    %    unnormalizedCellCycle: [-1 -1 -1 -1 0 5 10 15 20 25 30 35 40 45 50 55
    %        -1 -1 -1]
    %    normalizedCellCycle: [-1	-1	-1	-1	0	0.0909090909090909
    %        0.181818181818182 0.272727272727273	0.363636363636364
    %        0.454545454545455	0.545454545454545	0.636363636363636
    %        0.727272727272727	0.818181818181818	0.909090909090909
    %        1	-1	-1	-1]
    %    cellCycleLength: 60
    %
    % @test
    %
    % division = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0];
    % cellCycle = CellCycleInformation(division, 5);
    %     startingAnalysis: 5
    %     endAnalysis: 20
    %     unnormalizedCellCycle: [-1 -1 -1 -1 0 5 10 15 20 25 30 35 40 45 50 55
    %         0 5 10 15 -1 -1 -1]
    %     normalizedCellCycle: [1×23 double]
    %     cellCycleLength: [60 20]
    %
    
    dividing = find(division > 0);
    cellCycle.startingAnalysis = dividing(1) + 1;
    cellCycle.endAnalysis = dividing(end);

    cellCycle.unnormalizedCellCycle = ones(size(division)) * -1;
    cellCycle.normalizedCellCycle = ones(size(division)) * -1;
    cellCycle.cellCycleLength = [];

    for i = 1 : length(dividing) - 1
        cellCycle.cellCycleLength(i) = timeInterval * (dividing(i + 1) - dividing(i));

        for j = (dividing(i) + 1) : dividing(i + 1)
            cellCycle.unnormalizedCellCycle(j) = (j - dividing(i) - 1) * timeInterval;
            cellCycle.normalizedCellCycle(j) = cellCycle.unnormalizedCellCycle(j) / (cellCycle.cellCycleLength(i) - timeInterval);
        end
    end
    
end
