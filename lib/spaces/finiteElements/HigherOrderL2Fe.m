% HigherOrderL2Fe (subclass of NodalFiniteElement and BernsteinBezierPoly) Class
%   for discontinuous Lagrange Finite Elements of arbitrary order implemented by
%   Bernstein-Bezier polynomials.

classdef HigherOrderL2Fe < NodalFiniteElement & BernsteinBezierPoly
    %% implementation of abstract superclass properties
    properties (GetAccess='public', SetAccess='protected')
        order
        nodeDofLocations
        edgeDofLocations
        innerDofLocations
    end
    
    %% implementation of abstract superclass methods
    methods (Access = 'public')
        %% constructor
        function obj = HigherOrderL2Fe(order)
            arguments
                order (1,1) double {mustBeNonnegative}
            end
            
            obj.order = order;
            obj.nodeDofLocations = [];
            obj.edgeDofLocations = [];
            
            p = order;
            [x, y] = meshgrid((0:p)/p, (0:p)/p);
            z = x + y;
            ind = find(z < 1+10*eps);
            obj.innerDofLocations = [1-z(ind)'; x(ind)'; y(ind)'];
            if order == 0
                obj.innerDofLocations = [1;1;1]/3;
            end
            
            obj = obj.setExponents(round(getDofLocations(obj).coordinates*order));
            
            % contract dof locations to make sure none is on element boundary
            delta = 0.95;
            obj.innerDofLocations = delta*obj.innerDofLocations + (1-delta)/3;
        end
    end
end
