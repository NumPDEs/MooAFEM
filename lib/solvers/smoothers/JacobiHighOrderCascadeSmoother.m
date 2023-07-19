% JacobiHighOrderCascadeSmoother (subclass of MultilevelSmoother) multilevel
%   Jacobi smoother for higher order finite elements: smooth locally with
%   patchwise stiffness matrix on changed patches on every level, except for
%   coarsest level, where global p1-solving is done. Because of the lifting from
%   p1 to higher order between the coarsest and the second level, all patches
%   are assumed to be changed on the second level
%   
%   Reference: [Innerberger, Miraçi, Praetorius, Streitberger; Optimal
%   computational costs of AFEM with optimal local hp-robust multigrid solver].
%
% See also: MultilevelSmoother, OptimalVcycleMultigridSolver


classdef JacobiHighOrderCascadeSmoother < MultilevelSmoother
    properties (Access=protected)
        P
        loFes
        patchwiseA
        inclusionMatrix
        intergridMatrix
        freeDofs
        freeDofsOld
    end

    %% methods
    methods (Access=public)
        function obj = JacobiHighOrderCascadeSmoother(fes, blf, P)
            obj = obj@MultilevelSmoother(fes, blf);
            
            assert(fes.finiteElement.order > 1, ...
                'JacobiHighOrderCascadeSmoother only works for higher order finite elements.')

            obj.P = P;
            obj.loFes = FeSpace(fes.mesh, LowestOrderH1Fe(), 'dirichlet', fes.bnd.dirichlet);
        end

        function nonInvertedSmootherMatrix = setup(obj, A)
            setup@MultilevelSmoother(obj, A);
            
            L = obj.nLevels;
            obj.freeDofsOld = obj.freeDofs;
            obj.freeDofs = getFreeDofs(obj.fes);
            freeVertices = getFreeDofs(obj.loFes);
            
            hoFes = obj.fes;
            if L == 1
                % set up FE inclusion and lowest order matrix on coarsest level
                obj.inclusionMatrix = SpaceProlongation(obj.loFes, hoFes) ...
                    .matrix(getFreeDofs(obj.loFes), getFreeDofs(hoFes));
                nonInvertedSmootherMatrix = assemble(obj.blf, obj.loFes);
                nonInvertedSmootherMatrix = nonInvertedSmootherMatrix(freeVertices, freeVertices);
            else
                % note: patchwiseA is not already indexed by free vertices -> use freeVertices here
                obj.changedPatches{L} = freeVertices(obj.changedPatches{L}(freeVertices));
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeDofs, obj.freeDofsOld);
                nonInvertedSmootherMatrix = A;
                if L == 2
                    obj.patchwiseA{L} = assemblePatchwise(obj.blf, hoFes, ':');
                else
                    obj.patchwiseA{L} = assemblePatchwise(obj.blf, hoFes, obj.changedPatches{L});
                end
            end
        end

        function Cx = smooth(obj, x, k)
            % higher order global (k=2) or local (else) patchwise smooting
            Cx = obj.patchwiseA{k} \ x;
        end

        % order of FESpace per level: 1 -> p -> ... -> p
        function Px = prolongate(obj, x, k)
            if k == 2
                x = obj.inclusionMatrix' * x;
            end
            Px = obj.intergridMatrix{k} * x;
        end

        function Px = restrict(obj, x, k)
            Px = obj.intergridMatrix{k}' * x;
            if k == 2
                Px = obj.inclusionMatrix * Px;
            end
        end
    end
end