% BisectionEventData (subclass of event.EventData) Encapsulates data used for
%   bisection.

classdef (ConstructOnLoad) BisectionEventData < event.EventData
    properties (GetAccess='public', SetAccess='protected')
        bisectedEdges (:,1) logical
        bisecGroups (:,1) cell
        nRefinedElements (:,1) double
        nDescendants (:,1) double
        nInnerNodes (:,1) double
        nInnerEdges (:,1) double
    end
    
    properties (Dependent)
        nBisecGroups
    end
    
    methods (Access='public')
        function obj = BisectionEventData(bisectedEdges, bisecGroups)
            obj.bisectedEdges = bisectedEdges;
            obj.bisecGroups = bisecGroups;
            obj.nRefinedElements = cellfun(@(x) x.nMembers, bisecGroups);
            obj.nDescendants = cellfun(@(x) x.nDescendants, bisecGroups);
            obj.nInnerNodes = cellfun(@(x) x.nInnerNodes, bisecGroups);
            obj.nInnerEdges = cellfun(@(x) x.nInnerEdges, bisecGroups);
        end
    end
    
    methods    
        function val = get.nBisecGroups(obj)
            val = numel(obj.bisecGroups);
        end
    end
end

