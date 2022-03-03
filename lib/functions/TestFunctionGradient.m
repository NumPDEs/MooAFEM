% TestFunctionGradient (subclass of Evaluable) (Column-wise) matrix of 
%   derivatives of test functions on a finite element space.
%
%   f = TestFunctionGradient(fes) initializes derivatives of a test function on
%       the finite element space fes. The evaluation of this function is a
%       (column-wise) matrix of the values of the derivatives of all shape
%       functions at the given barycentric coordinate.
%
%   See also: Evaluable

classdef TestFunctionGradient < Evaluable
    %% properties
    properties (SetAccess='protected', GetAccess='public')
        fes
        mesh
    end
    
    %% methods
    methods (Access='public')
        function obj = TestFunctionGradient(fes)
            arguments
                fes (1,1) FeSpace
            end
            
            obj.fes = fes;
            obj.mesh = fes.mesh;
        end
        
        function val = eval(obj, bary, idx)
            arguments
                obj
                bary Barycentric2D
                idx {mustBeIndexVector}= ':'
            end
            
            dofs = getDofs(obj.fes);
            trafo = getAffineTransformation(obj.fes.mesh);
            nLocDofs = size(dofs.element2Dofs, Dim.Vector);
            
            Dphi = evalShapeFunctionDerivs(obj.fes.finiteElement, bary);
            val = vectorProduct(Dphi, trafo.DFinv(:,idx), [nLocDofs,2], [2,2]);
        end
    end
end
