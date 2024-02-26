% CRJacobiCascade (subclass of MultilevelSmoother) multilevel Jacobi smoother
%   for lowest order nonconforming Crouzeix-Raviart finite elements: smooth 
%   with diagonal of stiffness matrix on changed patches on every level.
%
% See also: MultilevelSmoother, OptimalVcycleMultigridSolver


classdef CRJacobiCascade < MultilevelSmoother
    properties (Access=protected)
        P
        patchwiseA
        intergridMatrix
        freeDofs
        freeDofsOld
    end
    
    %% event data
    properties (Access=private)
        listenerHandle
    end

    %% methods
    methods (Access=public)
        function obj = CRJacobiCascade(fes, blf, P)
            arguments
                fes
                blf
                P Prolongation
            end

            % TODO assert CR element
            assert(class(fes.finiteElement) == "LowestOrderCRFe" && ...
                class(P.feType) == "CR", ...
                ['CRJacobiSmoother only works for lowest order nonconforming', ... 
                'Crouzeix-Raviart finite elements.'])
            assert(fes.finiteElement.order == 1, ...
                'CRJacobiSmoother only works for lowest order finite elements.')
            
            obj = obj@MultilevelSmoother(fes, blf);
            obj.P = P;
        end

        function nonInvertedSmootherMatrix = setup(obj, A)
            setup@MultilevelSmoother(obj, A);
            
            L = obj.nLevels;
            obj.freeDofsOld = obj.freeDofs;
            obj.freeDofs = getFreeDofs(obj.fes);
            nonInvertedSmootherMatrix = A;
            
            if L > 1
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeDofs, obj.freeDofsOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(obj.freeDofs));
                obj.patchwiseA{L} = assemblePatchwise(obj.blf, obj.fes, obj.changedPatches{L});
            end

        end

        function Cx = smooth(obj, x, k)
            Cx = obj.patchwiseA{k} \ x;
        end

        function Px = prolongate(obj, x, k)
            Px = obj.intergridMatrix{k} * x;
        end

        function Px = restrict(obj, x, k)
            Px = obj.intergridMatrix{k}' * x;
        end
    end
end
