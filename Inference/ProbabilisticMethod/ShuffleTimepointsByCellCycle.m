function shuffledIndex = ShuffleTimepointsByCellCycle(unnormalizedCellCycle)

    cellCycleIndex = CellCycleToIndex(unnormalizedCellCycle);
    cellCycleIndex = cellCycleIndex';
    allCellCycleIndices = unique(cellCycleIndex);
    
    shuffledOrder = allCellCycleIndices(randperm(length(allCellCycleIndices)));
    
    shuffledIndex = [];
    for i = 1 : length(allCellCycleIndices)
        shuffledIndex = [shuffledIndex find(shuffledOrder(i) == cellCycleIndex)];
    end
        
end
