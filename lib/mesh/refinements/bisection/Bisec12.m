% Bisec12 (subclass of AbstracBisection) Encapsulate bisection of an element
%   where the 1st and 2nd edge are bisected.

classdef Bisec12 < AbstractBisection
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        innerNodes = []
        nInnerEdges = 2
        nDescendants = 3
    end
    
    %% methods
    methods (Access=public)
        function newElements = refineElement(~, oldNodes, edgeMidPts, ~)
            newElements = [oldNodes(3,:), edgeMidPts(1,:),   oldNodes(3,:); ...
                oldNodes(1,:),   oldNodes(2,:), edgeMidPts(1,:); ...
                edgeMidPts(1,:), edgeMidPts(2,:), edgeMidPts(2,:)];
        end
        
        function newEdges = createNewEdges(~, oldNodes, edgeMidPts, ~)
            newEdges = [edgeMidPts(1,:), edgeMidPts(2,:); ...
                oldNodes(3,:), edgeMidPts(1,:)];
            
        end
        
        function newElement2Edges = refineEdgeConnectivity(~, oldEdges, left, right, innerEdges)
            newElement2Edges = [oldEdges(3,:),      right(1,:), innerEdges(1,:); ...
                left(1,:),       left(2,:), innerEdges(2,:); ...
                innerEdges(1,:), innerEdges(2,:),      right(2,:)];
        end
        
        function newOrientation = refineEdgeOrientation(obj, old)
            [T, F] = getTrueFalseArrays(obj, size(old,2));
            newOrientation = [old(3,:), old(1,:),       T; ...
                old(1,:), old(2,:),       T; ...
                F,        F, old(2,:)];
        end
    end   
end
