% Bisec13 (subclass of AbstracBisection) Encapsulate bisection of an element
%   where the 1st and 2nd edge are bisected.

classdef Bisec13 < AbstractBisection
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        innerNodes = []
        nInnerEdges = 2
        nDescendants = 3
    end
    
    %% public methods
    methods (Access='public')
        function newElements = refineElement(~, oldNodes, edgeMidPts, ~)
            newElements = [edgeMidPts(1,:),   oldNodes(1,:),   oldNodes(2,:); ...
                oldNodes(3,:), edgeMidPts(1,:),   oldNodes(3,:); ...
                edgeMidPts(3,:), edgeMidPts(3,:), edgeMidPts(1,:)];
        end
        
        function newEdges = createNewEdges(~, oldNodes, edgeMidPts, ~)
            newEdges = [edgeMidPts(1,:), edgeMidPts(3,:); ...
                oldNodes(3,:), edgeMidPts(1,:)];
        end
        
        function newElement2Edges = refineEdgeConnectivity(~, oldEdges, left, right, innerEdges)
            newElement2Edges = [innerEdges(1,:),       left(1,:),   oldEdges(2,:); ...
                left(3,:), innerEdges(2,:), innerEdges(1,:); ...
                innerEdges(2,:),      right(3,:),      right(1,:)];
        end
        
        function newOrientation = refineEdgeOrientation(obj, old)
            [T, F] = getTrueFalseArrays(obj, size(old,2));
            newOrientation = [F, old(1,:), old(2,:); ...
                old(3,:),        T,        T; ...
                F, old(3,:), old(1,:)];
        end
    end
end
