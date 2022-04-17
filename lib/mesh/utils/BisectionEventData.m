% BisectionEventData (subclass of event.EventData) Encapsulates data used for
%   bisection.

classdef (ConstructOnLoad) BisectionEventData < event.EventData
    properties (GetAccess=public, SetAccess=protected)
        nElements (1,1) double
        bisectedEdges (:,1) logical
        bisection (:,1) cell
        nNonRefinedElements (1,1) double
        nRefinedElements (:,1) double
        nDescendants (:,1) double
        nInnerNodes (:,1) double
        nInnerEdges (:,1) double
    end
    
    properties (Access=protected)
        nonRefinedElements
        refinedElements
    end
    
    properties (Dependent)
        nBisecMethods
    end
    
    methods (Access=public)
        function obj = BisectionEventData(nElements, bisectedEdges, bisecMethods, refinedElements)
            arguments
                nElements
                bisectedEdges
            end
            arguments (Repeating)
                bisecMethods
                refinedElements
            end
            
            obj.nElements = nElements;
            obj.bisectedEdges = bisectedEdges;
            obj.bisection = bisecMethods;
            obj.refinedElements = refinedElements;
            
            obj.nRefinedElements = cellfun(@(x) numel(x), refinedElements);
            obj.nDescendants = cellfun(@(x) x.nDescendants, bisecMethods);
            obj.nInnerNodes = cellfun(@(x) size(x.innerNodes,2), bisecMethods);
            obj.nInnerEdges = cellfun(@(x) x.nInnerEdges, bisecMethods);
            
            elements = true(nElements, 1);
            elements(horzcat(obj.refinedElements{:})) = false;
            obj.nonRefinedElements = find(elements);
            obj.nNonRefinedElements = numel(obj.nonRefinedElements);
        end
        
        function children = getNChildrenPerElement(obj)
            children = accumarray(horzcat(obj.refinedElements{:})', ...
                repelem(obj.nDescendants, obj.nRefinedElements), ...
                [obj.nElements, 1], [], 1);
        end
        
        % 0 = non-refined, 1,... = refined according to bisection{k}
        function idx = getRefinedElementIdx(obj, k)
            assert((0 <= k) && (k <= obj.nBisecMethods))
            if k == 0
                idx = obj.nonRefinedElements;
            else
                idx = obj.refinedElements{k};
            end
        end
    end
    
    methods    
        function val = get.nBisecMethods(obj)
            val = numel(obj.bisection);
        end
    end
end

