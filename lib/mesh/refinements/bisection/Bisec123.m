% Bisec123 (subclass of AbstracBisection) Encapsulate bisection of an element
%   where all edges are bisected.

classdef Bisec123 < AbstractBisection
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        innerNodes = []
        nInnerEdges = 3
        nDescendants = 4
    end
    
    %% methods
    methods (Access='public')
        function newElements = refineElement(~, oldNodes, edgeMidPts, ~)
            newElements = [edgeMidPts(1,:),   oldNodes(1,:), edgeMidPts(1,:),   oldNodes(3,:); ...
                oldNodes(3,:), edgeMidPts(1,:),   oldNodes(2,:), edgeMidPts(1,:); ...
                edgeMidPts(3,:), edgeMidPts(3,:), edgeMidPts(2,:), edgeMidPts(2,:)];
        end
        
        function newEdges = createNewEdges(~, oldNodes, edgeMidPts, ~)
            newEdges = [edgeMidPts(1,:), edgeMidPts(3,:), edgeMidPts(2,:); ...
                oldNodes(3,:), edgeMidPts(1,:), edgeMidPts(1,:)];
        end
        
        function newElement2Edges = refineEdgeConnectivity(~, ~, left, right, innerEdges)
            newElement2Edges = [innerEdges(1,:),      left(1,:),      right(1,:), innerEdges(1,:); ...
                left(3,:),innerEdges(2,:),       left(2,:), innerEdges(3,:); ...
                innerEdges(2,:),     right(3,:), innerEdges(3,:),      right(2,:)];
        end
        
        function newOrientation = refineEdgeOrientation(obj, old)
            [T, F] = getTrueFalseArrays(obj, size(old,2));
            newOrientation = [F, old(1,:), old(1,:),        T; ...
                old(3,:),        T, old(2,:),        T; ...
                F, old(3,:),        F, old(2,:)];
        end
    end 
end
