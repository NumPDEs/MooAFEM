% P1JacobiCascade (subclass of MultilevelSmoother) multilevel Jacobi smoother
%   for lowest order finite elements: smooth with diagonal of stiffness matrix
%   on changed patches on every level.
%
% See also: MultilevelSmoother, OptimalVcycleMultigridSolver


classdef P1JacobiCascade < MultilevelSmoother
    properties (Access=protected)
        P
        inverseDiagonal
        intergridMatrix
        freeVertices
        freeVerticesOld
    end
    
    %% event data
    properties (Access=private)
        listenerHandle
    end

    %% methods
    methods (Access=public)
        function obj = P1JacobiCascade(fes, blf, P)
            arguments
                fes
                blf
                P Prolongation
            end

            assert(fes.finiteElement.order == 1, ...
                'P1JacobiSmoother only works for lowest order finite elements.')
            
            obj = obj@MultilevelSmoother(fes, blf);
            obj.P = P;
        end

        function nonInvertedSmootherMatrix = setup(obj, A)
            setup@MultilevelSmoother(obj, A);
            
            L = obj.nLevels;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.fes);
            
            nonInvertedSmootherMatrix = A;
            obj.inverseDiagonal{L} = full(diag(A)).^(-1);
            
            if L >= 2
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(obj.freeVertices));
                obj.inverseDiagonal{L} = obj.inverseDiagonal{L}(obj.changedPatches{L});
            end
        end

        function Cx = smooth(obj, x, k)
            idx = obj.changedPatches{k};
            Cx = zeros(size(x));
            Cx(idx,:) = obj.inverseDiagonal{k} .* x(idx,:);
        end

        function Px = prolongate(obj, x, k)
            Px = obj.intergridMatrix{k} * x;
        end

        function Px = restrict(obj, x, k)
            Px = obj.intergridMatrix{k}' * x;
        end
    end
end