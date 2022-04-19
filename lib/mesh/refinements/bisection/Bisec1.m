% Bisec1 (subclass of AbstracBisection) Encapsulate bisection of an element
%   where only the refinement edge (1st edge) is bisected.

classdef Bisec1 < AbstractBisection
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        innerNodes = []
        nInnerEdges = 1
        nDescendants = 2
    end
    
    %% methods
    methods (Access=public)
        function newElements = refineElement(~, oldNodes, edgeMidPts, ~)
            newElements = [oldNodes(3,:),   oldNodes(2,:); ...
                oldNodes(1,:),   oldNodes(3,:); ...
                edgeMidPts(1,:), edgeMidPts(1,:)];
        end
        
        function newEdges = createNewEdges(~, oldNodes, edgeMidPts, ~)
            newEdges = [edgeMidPts(1,:); oldNodes(3,:)];
        end
        
        function newElement2Edges = refineEdgeConnectivity(~, oldEdges, left, right, innerEdges)
            newElement2Edges = [oldEdges(3,:), oldEdges(2,:); ...
                left(1,:),    innerEdges; ...
                innerEdges,    right(1,:)];
        end
        
        function newOrientation = refineEdgeOrientation(obj, old)
            [T, F] = getTrueFalseArrays(obj, size(old,2));
            newOrientation = [old(3,:), old(2,:); ...
                old(1,:),        T; ...
                F, old(1,:)];
        end
    end
end
