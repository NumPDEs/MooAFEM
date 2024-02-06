% LowestOrderL2Fe (value class, subclass of NodalFiniteElement) Class for
%   discontinuous Lagrange Finite Elements of lowest order. The implementation
%   is specialized to the lowest order case.

classdef LowestOrderL2Fe < NodalFiniteElement
    %% implementation of abstract superclass properties
    properties (GetAccess=public, SetAccess=protected)
        order
        nodeDofLocations
        edgeDofLocations
        innerDofLocations
    end
    
    %% implementation of abstract superclass methods
    methods (Access=public)
        function obj = LowestOrderL2Fe()
            obj.order = 0;
            obj.nodeDofLocations = [];
            obj.edgeDofLocations = [];
            obj.innerDofLocations = [1;1;1]/3;
        end
        
        function val = evalShapeFunctions(obj, bary)
            arguments
                obj %#ok<INUSA>
                bary Barycentric2D
            end
            
            val = ones(1, 1, bary.nNodes);
        end
        
        function val = evalShapeFunctionDerivs(obj, bary)
            arguments
                obj %#ok<INUSA>
                bary Barycentric2D
            end
            
            val = zeros(2, 1, bary.nNodes);
        end
        
        function val = evalShapeFunctionHessian(obj, bary)
            arguments
                obj %#ok<INUSA>
                bary Barycentric2D
            end
            
            val = zeros(4, 1, bary.nNodes);
        end
    end
end
