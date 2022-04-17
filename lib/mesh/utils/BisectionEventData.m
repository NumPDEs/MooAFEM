% BisectionEventData (subclass of event.EventData) Encapsulates data used for
%   bisection.

classdef (ConstructOnLoad) BisectionEventData < event.EventData
    properties (GetAccess='public', SetAccess='protected')
        nElements (1,1) double
        bisectedEdges (:,1) logical
        bisecMethods (:,1) cell
        nonRefinedElements (:,1) double
        nNonRefinedElements (1,1) double
        refinedElements (:,1) cell
        nRefinedElements (:,1) double
        nDescendants (:,1) double
        nInnerNodes (:,1) double
        nInnerEdges (:,1) double
    end
    
    properties (Dependent)
        nBisecMethods
    end
    
    methods (Access='public')
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
            obj.bisecMethods = bisecMethods;
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
    end
    
    methods    
        function val = get.nBisecMethods(obj)
            val = numel(obj.bisecMethods);
        end
        
        function children = getNChildrenPerElement(obj)
            children = ones(obj.nElements, 1);
            for k = 1:obj.nBisecMethods
                children(obj.refinedElements{k}) = obj.nDescendants(k);
            end
        end
    end
end

