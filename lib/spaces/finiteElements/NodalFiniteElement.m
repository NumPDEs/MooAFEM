% NodalFiniteElement (abstract subclass of FiniteElement) Abstract class for
%   nodal Finite Elements, where the dofs can be located in the element (e.g.,
%   the dual basis is comprised only of point-evaluations). Locations must be
%   given as n-dimensional barycentric coordinates.

classdef NodalFiniteElement < FiniteElement
    %% properties
    properties (Abstract, GetAccess=public, SetAccess=protected)
        nodeDofLocations
        edgeDofLocations
        innerDofLocations
    end
    
    %% implementation of abstract superclass method
    methods
        function connectivity = getDofConnectivity(obj)
            % connectivity is given by number of dof locations
            connectivity = [length(obj.nodeDofLocations), ...
                size(obj.edgeDofLocations, Dim.Elements), ...
                size(obj.innerDofLocations, Dim.Elements)];
        end
        
        function bary = getDofLocations(obj)
            nDofsPerPart = getDofConnectivity(obj);
            nLocations = nDofsPerPart * [3;3;1];
            locations = zeros(3, nLocations);
            n = 0;
            
            nNodeDofs = nDofsPerPart(1);
            if nNodeDofs ~= 0
                locations(:,1:3) = eye(3);
                n = n + 3;
            end
            
            nEdgeDofs = nDofsPerPart(2);
            if nEdgeDofs ~= 0
                bary1d = Barycentric1D(obj.edgeDofLocations);
                for i = 1:3
                    bary2d = Barycentric2D.from1D(bary1d, i);
                    locations(:,(n+1):(n+nEdgeDofs)) = bary2d.coordinates;
                    n = n + nEdgeDofs;
                end
            end
            
            nInnerDofs = nDofsPerPart(3);
            if nInnerDofs ~= 0
                locations(:,(n+1):end) = obj.innerDofLocations;
            end
            
            bary = Barycentric2D(locations);
        end
        
        function val = evalEdgeShapeFunctions(obj, ~)
            error('No implementation of method ''evalEdgeShapeFunctions''.')
        end
    end
end
