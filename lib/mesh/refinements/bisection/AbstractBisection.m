% AbstractBisection (value class) Abstract interface for bisection rules on 2D
%   triangles.

classdef AbstractBisection
    %% properties
    properties (Abstract, GetAccess='public', SetAccess='protected')
        innerNodes (3,:) double
        nInnerEdges (1,1) double
        nDescendants (1,1) double
    end
    
    properties (Dependent)
        nInnerNodes
        nMembers
    end
    
    properties (GetAccess='public', SetAccess='protected')
        elementIdx (:,1) double
    end
    
    %% methods
    methods (Access='public')
        function obj = AbstractBisection(elementIndices)
            obj.elementIdx = elementIndices;
        end
    end
    
    methods
        function n = get.nInnerNodes(obj)
            n = size(obj.innerNodes, Dim.Elements);
        end
        
        function n = get.nMembers(obj)
            n = length(obj.elementIdx);
        end
    end
    
    methods (Access='protected')
        function [T, F] = getTrueFalseArrays(obj)
            T = true(1, obj.nMembers);
            F = false(1, obj.nMembers);
        end
    end
    
    methods (Abstract, Access='public')
        newElements = refineElement(obj, elements, edge2newNodes, innerNodes)
        newEdges = createNewEdges(obj, elements, edge2newNodes, innerNodes)
        newElement2Edges = refineEdgeConnectivity(obj, element2edges, left, right, innerEdges)
        newOrientation = refineEdgeOrientation(obj, oldOrientation)
    end
end
