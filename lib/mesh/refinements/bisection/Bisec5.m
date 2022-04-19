% Bisec5 (subclass of AbstracBisection) Encapsulate bisection of an element such
%   that an interior node is present.

classdef Bisec5 < AbstractBisection
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        innerNodes = [1;1;2]/4;
        nInnerEdges = 6
        nDescendants = 6
        interiorNodes (:,1) double
    end
    
    %% methods
    methods (Access=public)
        function newElements = refineElement(~, oldNodes, edgeMidPts, innerNodes)
            newElements = [edgeMidPts(3,:),   oldNodes(3,:),   oldNodes(1,:), edgeMidPts(1,:), edgeMidPts(2,:), edgeMidPts(1,:); ...
                edgeMidPts(1,:), edgeMidPts(3,:), edgeMidPts(1,:),   oldNodes(2,:),   oldNodes(3,:), edgeMidPts(2,:); ...
                innerNodes,      innerNodes, edgeMidPts(3,:), edgeMidPts(2,:),      innerNodes,     innerNodes];
        end
        
        function newEdges = createNewEdges(~, oldNodes, edgeMidPts, innerNodes)
            newEdges = [edgeMidPts(1,:), edgeMidPts(2,:), edgeMidPts(1,:), edgeMidPts(2,:), oldNodes(3,:), edgeMidPts(3,:); ...
                edgeMidPts(3,:), edgeMidPts(1,:),      innerNodes,      innerNodes,    innerNodes,      innerNodes];
        end
        
        function newElement2Edges = refineEdgeConnectivity(~, ~, left, right, innerEdges)
            newElement2Edges = [innerEdges(1,:),       left(3,:),       left(1,:),      right(1,:),      right(2,:), innerEdges(2,:); ...
                innerEdges(3,:), innerEdges(6,:), innerEdges(1,:),       left(2,:), innerEdges(5,:), innerEdges(4,:); ...
                innerEdges(6,:), innerEdges(5,:),      right(3,:), innerEdges(2,:), innerEdges(4,:), innerEdges(3,:)];
        end
        
        function newOrientation = refineEdgeOrientation(obj, old)
            [T, F] = getTrueFalseArrays(obj, size(old,2));
            newOrientation = [T, old(3,:), old(1,:), old(1,:), old(2,:), T; ...
                F,        F,        F, old(2,:),        F, F; ...
                T,        T, old(3,:),        F,        T, T];
        end
    end
end
