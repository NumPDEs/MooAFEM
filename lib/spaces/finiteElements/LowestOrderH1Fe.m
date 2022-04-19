% LowestOrderH1Fe (subclass of NodalFiniteElement) Class for globally continuous
%   Lagrange Finite Elements of lowest order. The implementation is specialized
%   to the lowest order case.

classdef LowestOrderH1Fe < NodalFiniteElement
    %% implementation of abstract superclass properties
    properties (GetAccess=public, SetAccess=protected)
        order
        nodeDofLocations
        edgeDofLocations
        innerDofLocations
    end
    
    %% implementation of abstract superclass methods
    methods (Access=public)
        function obj = LowestOrderH1Fe()
            obj.order = 1;
            obj.nodeDofLocations = 1;
            obj.edgeDofLocations = [];
            obj.innerDofLocations = [];
        end
        
        function val = evalShapeFunctions(obj, bary)
            arguments
                obj
                bary Barycentric2D
            end
            
            val = bary.arangeAlongDim(3);
        end
        
        function val = evalShapeFunctionDerivs(obj, bary)
            arguments
                obj
                bary Barycentric2D
            end
            
            val = repmat([-1;1;0;-1;0;1], 1, 1, bary.nNodes);
        end
        
        function val = evalShapeFunctionHessian(obj, bary)
            arguments
                obj
                bary Barycentric2D
            end
            
            val = zeros(12, 1, bary.nNodes);
        end
        
        function val = evalEdgeShapeFunctions(obj, bary)
            arguments
                obj
                bary Barycentric1D
            end
            
            val = bary.arangeAlongDim(3);
        end
    end 
end
