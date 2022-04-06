% BisectionEventData (subclass of event.EventData) Encapsulates data used for
%   bisection.

classdef (ConstructOnLoad) BisectionEventData < event.EventData
    properties (GetAccess='public', SetAccess='protected')
        bisectedEdges (:,1) logical
        bisecGroups (:,1) cell
    end
    
    properties (Dependent)
        nBisecGroups
    end
    
    methods (Access='public')
        function obj = BisectionEventData(bisectedEdges, bisecGroups)
            obj.bisectedEdges = bisectedEdges;
            obj.bisecGroups = bisecGroups;
        end
    end
    
    methods    
        function val = get.nBisecGroups(obj)
            val = numel(obj.bisecGroups);
        end
    end
end

