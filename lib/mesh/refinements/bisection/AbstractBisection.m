% AbstractBisection (value class) Abstract interface for bisection rules on 2D
%   triangles.

classdef AbstractBisection
    %% properties
    properties (Abstract, GetAccess='public', SetAccess='protected')
        innerNodes (3,:) double
        nInnerEdges (1,1) double
        nDescendants (1,1) double
    end
    
    %% methods
    methods (Access='protected')
        function [T, F] = getTrueFalseArrays(obj, n)
            T = true(1, n);
            F = false(1, n);
        end
    end
    
    methods (Abstract, Access='public')
        newElements = refineElement(obj, elements, edge2newNodes, innerNodes)
        newEdges = createNewEdges(obj, elements, edge2newNodes, innerNodes)
        newElement2Edges = refineEdgeConnectivity(obj, element2edges, left, right, innerEdges)
        newOrientation = refineEdgeOrientation(obj, oldOrientation)
    end
end
