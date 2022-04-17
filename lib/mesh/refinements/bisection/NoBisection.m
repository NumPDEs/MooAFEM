% NoBisection (subclass of AbstracBisection) Encapsulate 'bisection' of an
%   element which is not bisected.

classdef NoBisection < AbstractBisection
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        innerNodes = []
        nInnerEdges = 0
        nDescendants = 1
    end
    
    %% methods
    methods (Access='public')
        function newElements = refineElement(~, oldNodes, ~, ~)
            newElements = oldNodes;
        end
        
        function newEdges = createNewEdges(~, ~, ~, ~)
            newEdges = [];
        end
        
        function newElement2Edges = refineEdgeConnectivity(~, oldEdges, ~, ~, ~)
            newElement2Edges = oldEdges;
        end
        
        function newOrientation = refineEdgeOrientation(obj, old)
            newOrientation = old;
        end
    end
end
