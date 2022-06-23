% Gradient (subclass of Evaluable) Gradient of a finite element function.
%
%   f = Gradient(u) initializes the gradient of the finite element function u.
%       Note that this only provides another view on the data of u and, in
%       particular, no data is copied.
%
%   See also: Evaluable

classdef Gradient < Evaluable
    %% properties
    properties (SetAccess=protected, GetAccess=public)
        u
        mesh
    end
    
    %% methods
    methods (Access=public)
        function obj = Gradient(u)
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
                idx {mustBeIndexVector} = ':'
            end
            
            assert(obj.u.dataIsValid, 'Gradient:dataNotValid', ...
                'Data cannot be evaluated, since underlying Mesh / FeSpace has changed.')
            
            fes = obj.u.fes;
            dofs = getDofs(fes);
            trafo = getAffineTransformation(fes.mesh);
            nLocDofs = size(dofs.element2Dofs, Dim.Vector);
            
            Dphi = evalShapeFunctionDerivs(fes.finiteElement, bary);
            I = dofs.element2Dofs(:,idx);
            uElements = reshape(obj.u.data(I), size(I));
            
            val = vectorProduct(uElements, Dphi, [nLocDofs, 1]', [nLocDofs, 2]);
            val = vectorProduct(val, trafo.DFinv(:,idx), [1,2], [2,2]);
        end
    end
end
