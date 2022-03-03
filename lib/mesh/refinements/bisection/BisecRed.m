% BisecRed (subclass of AbstracBisection) Encapsulate red bisection of an
%   element.

classdef BisecRed < AbstractBisection
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        innerNodes = []
        nInnerEdges = 3
        nDescendants = 4
    end
    
    %% methods
    methods (Access='public')
        function newElements = refineElement(~, oldNodes, edgeMidPts, ~)
            newElements = [oldNodes(1,:), edgeMidPts(1,:), edgeMidPts(3,:), edgeMidPts(2,:); ...
                edgeMidPts(1,:),   oldNodes(2,:), edgeMidPts(2,:), edgeMidPts(3,:); ...
                edgeMidPts(3,:), edgeMidPts(2,:),   oldNodes(3,:), edgeMidPts(1,:)];
        end
        
        function newEdges = createNewEdges(~, ~, edgeMidPts, ~)
            newEdges = [edgeMidPts(1,:), edgeMidPts(2,:), edgeMidPts(3,:); ...
                edgeMidPts(3,:), edgeMidPts(1,:), edgeMidPts(2,:)];
        end
        
        function newElement2Edges = refineEdgeConnectivity(~, ~, left, right, innerEdges)
            newElement2Edges = [left(1,:),      right(1,:), innerEdges(3,:), innerEdges(3,:); ...
                innerEdges(1,:),       left(2,:),      right(2,:), innerEdges(1,:); ...
                right(3,:), innerEdges(2,:),       left(3,:), innerEdges(2,:)];
        end
        
        function newOrientation = refineEdgeOrientation(obj, old)
            [T, F] = getTrueFalseArrays(obj);
            newOrientation = [old(1,:), old(1,:),        F, T; ...
                F, old(2,:), old(2,:), T; ...
                old(3,:),        F, old(3,:), T];
        end
    end
end
