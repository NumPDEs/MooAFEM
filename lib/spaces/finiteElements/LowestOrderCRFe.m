% LowestOrderCRFe (subclass of NodalFiniteElement) Class for nonconforming 
%   Crouzeix-Raviart Finite Elements of lowest order. The implementation 
%   is specialized to the lowest order case.

classdef LowestOrderCRFe < NodalFiniteElement
    %% implementation of abstract superclass properties
    properties (GetAccess=public, SetAccess=protected)
        order
        nodeDofLocations
        edgeDofLocations
        innerDofLocations
    end
    
    %% implementation of abstract superclass methods
    methods (Access=public)
        function obj = LowestOrderCRFe()
            obj.order = 1;
            obj.nodeDofLocations = [];
            obj.edgeDofLocations = [0.5; 0.5];
            obj.innerDofLocations = [];
        end
        
        function val = evalShapeFunctions(obj, bary)
            arguments
                obj %#ok<INUSA>
                bary Barycentric2D
            end
            
            val = 1 - 2 * bary.arangeAlongDim(3);
            % reorder due to wrong local indexing of edges
            val = val([3, 1, 2],:,:);
        end
        
        function val = evalShapeFunctionDerivs(obj, bary)
            arguments
                obj %#ok<INUSA>
                bary Barycentric2D
            end

            % reorder due to wrong local indexing of edges
            gradientsS1 = [-1;1;0;-1;0;1];
            % gradientsNC = -2 * gradientsS1([5;1;3;6;2;4]);
            gradientsNC = -2 * gradientsS1;
            
            val = repmat(gradientsNC, 1, 1, bary.nNodes);
        end
        
        function val = evalShapeFunctionHessian(obj, bary)
            arguments
                obj %#ok<INUSA>
                bary Barycentric2D
            end
            
            val = zeros(12, 1, bary.nNodes);
        end
        
        function val = evalEdgeShapeFunctions(obj, bary)
            arguments
                obj %#ok<INUSA>
                bary Barycentric1D
            end
            
            val = 1 - 2 * bary.arangeAlongDim(3);
            val = val([3, 1, 2],:,:);
        end
    end 
end
