% JacobiLowOrderCascadeSmoother (subclass of MultilevelSmoother) multilevel
%   Jacobi smoother for higher order finite elements: smooth locally with
%   diagonal of P1 stiffness matrix on changed patches on every level, except
%   for finest level, where global patchwise higher order smoothing is done.
%   
%   Reference: [Innerberger, Miraçi, Praetorius, Streitberger; Optimal
%   computational costs of AFEM with optimal local hp-robust multigrid solver].
%
% See also: MultilevelSmoother, OptimalVcycleMultigridSolver


classdef JacobiLowOrderCascadeSmoother < MultilevelSmoother
    properties (Access=protected)
        P
        loFes
        p1Smoother
        hoSmoother
        inclusionMatrix
        intergridMatrix
        freeVertices
        freeVerticesOld
    end

    %% methods
    methods (Access=public)
        function obj = JacobiLowOrderCascadeSmoother(fes, blf)
            obj = obj@MultilevelSmoother(fes, blf);
            
            assert(fes.finiteElement.order > 1, ...
                'JacobiLowOrderCascadeSmoother only works for higher order finite elements.')

            mesh = fes.mesh;
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', fes.bnd.dirichlet);
            obj.P = Prolongation.chooseFor(obj.loFes);
            
        end

        function nonInvertedSmootherMatrix = setup(obj, A)
            setup@MultilevelSmoother(obj, A);
            
            L = obj.nLevels;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.loFes);
            
            p1Matrix = assemble(obj.blf, obj.loFes);
            p1Matrix = p1Matrix(obj.freeVertices, obj.freeVertices);
            obj.p1Smoother{L} = full(diag(p1Matrix)).^(-1);
            nonInvertedSmootherMatrix = p1Matrix;
            
            if L >= 2
                % intermediate level smoothers and transfer operators
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(obj.freeVertices));
                obj.p1Smoother{L} = obj.p1Smoother{L}(obj.changedPatches{L});
                % finest level smoother and transfer operator
                obj.hoSmoother = assemblePatchwise(obj.blf, obj.fes, ':');
                obj.inclusionMatrix = SpaceProlongation(obj.loFes, obj.fes) ...
                    .matrix(getFreeDofs(obj.loFes), getFreeDofs(obj.fes));
            end
        end

        function Cx = smooth(obj, x, k)
            if k == obj.nLevels
                % higher order global patchwise smooting
                Cx = obj.hoSmoother \ x;
            else
                % local p1 patchwise smoothing
                idx = obj.changedPatches{k};
                Cx = zeros(size(x));
                Cx(idx,:) = obj.p1Smoother{k} .* x(idx,:);
            end
        end

        function Px = prolongate(obj, x, k)
            Px = obj.intergridMatrix{k} * x;
            if k == obj.nLevels
                Px = obj.inclusionMatrix' * Px;
            end
        end

        function Px = restrict(obj, x, k)
            if k == obj.nLevels
                x = obj.inclusionMatrix * x;
            end
            Px = obj.intergridMatrix{k}' * x;
        end
    end
end