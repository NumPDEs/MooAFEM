% Hessian (subclass of Evaluable) Hessian of a finite element function.
%
%   f = Hessian(u) initializes the Hessian of the finite element function u.
%       Note that this only provides another view on the data of u and, in
%       particular, no data is copied.
%
%   See also: Evaluable

classdef Hessian < Evaluable
    %% properties
    properties (SetAccess=protected, GetAccess=public)
        u
        mesh
    end
    
    %% methods
    methods (Access=public)
        function obj = Hessian(u)
            arguments
                u (1,1) FeFunction
            end
            
            obj.u = u;
            obj.mesh = u.mesh;
        end
        
        function val = eval(obj, bary, idx)
            arguments
                obj
                bary Barycentric2D
                idx {mustBeIndexVector}= ':'
            end
            
            assert(obj.u.dataIsValid, 'Hessian:dataNotValid', ...
                'Data cannot be evaluated, since underlying Mesh / FeSpace has changed.')
            
            fes = obj.u.fes;
            dofs = getDofs(fes);
            trafo = getAffineTransformation(fes.mesh);
            nLocDofs = size(dofs.element2Dofs, Dim.Vector);
            
            D2phi = evalShapeFunctionHessian(fes.finiteElement, bary);
            I = dofs.element2Dofs(:,idx);
            uElements = reshape(obj.u.data(I), size(I));
            
            % D^2 phi = DF^{-T} D^2 \hat{phi} DF^{-1}
            val = vectorProduct(uElements, D2phi, [nLocDofs, 1]', [nLocDofs, 4]);
            val = vectorProduct(val, trafo.DFinv(:,idx), [2,2], [2,2]);
            val = vectorProduct(trafo.DFinv(:,idx), val, [2,2]', [2,2]);
        end
    end
end
