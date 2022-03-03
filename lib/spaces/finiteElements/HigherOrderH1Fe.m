% HigherOrderH1Fe (subclass of NodalFiniteElement and BernsteinBezierPoly) Class
%   for globally continuous Lagrange Finite Elements of arbitrary order
%   implemented by Bernstein-Bezier polynomials.

classdef HigherOrderH1Fe < NodalFiniteElement & BernsteinBezierPoly
    %% implementation of abstract superclass properties
    properties (GetAccess='public', SetAccess='protected')
        order
        nodeDofLocations
        edgeDofLocations
        innerDofLocations
    end
    
    %% implementation of abstract superclass methods
    methods (Access = 'public')
        function obj = HigherOrderH1Fe(order)
            arguments
                order (1,1) double {mustBePositive}
            end
            
            obj.order = order;
            obj.nodeDofLocations = 1;
            
            if order >= 2
                obj.edgeDofLocations = (1:(order-1)) / order;
                obj.edgeDofLocations = [obj.edgeDofLocations; 1-obj.edgeDofLocations];
            end
            
            if order >= 3
                [x, y] = meshgrid((1:(order-1)) / order, (1:(order-1)) / order);
                z = x + y;
                ind = find(z < 1-10*eps);
                obj.innerDofLocations = [1-z(ind)'; x(ind)'; y(ind)'];
            end
            
            obj = obj.setExponents(round(getDofLocations(obj).coordinates*order));
        end
        
        function val = evalEdgeShapeFunctions(obj, bary)
            arguments
                obj
                bary Barycentric1D
            end
            
            bary2d = Barycentric2D.from1D(bary, 1);
            elementVal = evalShapeFunctions(obj, bary2d);
            nEdgeDofs = size(obj.edgeDofLocations, Dim.Elements);
            edgeDofs = [1,2,4:(3+nEdgeDofs)];
            val = elementVal(edgeDofs,:,:);
        end
    end
end
