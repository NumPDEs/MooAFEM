% FiniteElement (abstract class) Abstract class for Finite Elements. A FE must
%   know how to evaluate (derivatives of) shape functions on one element and how
%   the dofs are connected.

classdef FiniteElement
    methods (Abstract, Access=public)
        % evalShapeFunctions Evaluate local shapefunctions at barycentric
        %   coordinates
        evalShapeFunctions(obj, bary)
        
        % evalShapeFunctionDerivs Evaluate derivatives of local
        %   shapefunctions at barycentric coordinates
        evalShapeFunctionDerivs(obj, bary)
        
        % evalShapeFunctionHessian Evaluate Hessian matrix of local
        %   shapefunctions at barycentric coordinates
        evalShapeFunctionHessian(obj, bary)
        
        % evalEdgeShapeFunctions Evaluate local shapefunctions at barycentric
        %   coordinates on edges
        evalEdgeShapeFunctions(obj, bary)
        
        % getDofConnectivity Return how many degrees of freedom are at nodes /
        %   edges / elements
        getDofConnectivity(obj)
    end
end
