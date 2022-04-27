classdef Information

    %
    % @description: information to be extracted from experiments.
    %

    % The former organizes the information seperately by lineage, and the latter
    % gives the overall combined statistics.
    properties
        lineageSpecific = [];
        overallStatistics = [];
    end

    methods
        
        function obj = Information(lineage)
        
            if nargin == 0
                return;
            end

            statistics = PromoterActivityStatistics(lineage.inferredLeakage, ...
                lineage.gapThreshold);
                        
            obj.lineageSpecific.id = cellstr(lineage.prefix);
            obj.lineageSpecific.timeInterval = lineage.timeInterval;
            obj.lineageSpecific.on = statistics.on;
            obj.lineageSpecific.off = statistics.off;
            obj.lineageSpecific.amplitude = statistics.amplitude;
            obj.lineageSpecific.burstSize = statistics.burstSize;
            obj.lineageSpecific.leakage = lineage.inferredLeakage;
            obj.lineageSpecific.timepoints = lineage.timepoints;
            obj.lineageSpecific.cellSize = lineage.cellSize;
            obj.lineageSpecific.fluorescence = lineage.fluorescence;
            obj.lineageSpecific.fluorescenceBeforeSignalFilter = ...
                lineage.fluorescenceBeforeSignalFilter;
            obj.lineageSpecific.division = lineage.division;
            obj.lineageSpecific.cellCycleLength = ...
                lineage.cellCycle.cellCycleLength;

            if isfield(lineage, 'realFluorescence')
                obj.lineageSpecific.realFluorescence = lineage.inferredRealFluorescence;
            end
            
        end

        function obj = OverallStatistics(obj)

            properties = fields(obj.lineageSpecific);
            obj.overallStatistics = CombineAttribute(obj, properties);
        end

        function receiver = Merge(receiver, sender)
        
            properties = fields(sender);
            for i = 1 : length(properties)
                
                previous = getfield(receiver, properties{i});
                if size(previous, 1) > 1
                    previous = previous';
                end
                
                current = getfield(sender, properties{i});
                if size(current, 1) > 1
                    current = current';
                end
                
                receiver = setfield(receiver, properties{i}, [previous current]);
            end
            
        end

    end
end
