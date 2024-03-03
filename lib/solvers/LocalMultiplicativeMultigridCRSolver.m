% LocalMultiplicativeMultigridCRSolver (subclass of IterativeSolver) Solves linear
%   equations iteratively, using a geometric multigrid solver with specified
%   MultiLevelSmoother.
%   One iteration = Vcycle(0,1) (no pre-/one post-smoothing step).
%
%   To obtain the solver described in [Innerberger, Miraçi, Praetorius,
%   Streitberger; Optimal computational costs of AFEM with optimal local
%   hp-robust multigrid solver], choose the smoother as DiagonalJacobi (p=1),
%   LowOrderBlockJacobi (p > 1, 1->1->p iteration), or HighOrderBlockJacobi
%   (p > 1, 1->p->p iteration).
%   This solver also comes with built-in estimator for the algebraic error,
%   which is also returned by this multigrid solver.
%
% See also IterativeSolver, MultigridSolver, MultiLevelSmoother, DiagonalJacobi,
%   LowOrderBlockJacobi, HighOrderBlockJacobi

classdef LocalMultiplicativeMultigridCRSolver < MultigridSolver
    %% methods
    methods
        function obj = LocalMultiplicativeMultigridCRSolver(smoother)
            obj = obj@MultigridSolver(smoother);
        end
    end
        
    methods (Access=protected)
        function [Cx, algEstimator] = cycle(obj, x)
            % wrapper for Vcycle algorithm to match interface of superclass
            [Cx, algEstimator] = obj.Vcycle(x);
        end

        function [Cx, algEstimator] = Vcycle(obj, x) 
            % one step of the Vcycle
            assert(isequal(size(x, 1), size(obj.A, 1)), ... 
                'Setup for multilevel iteration not complete!')
            
            % if there is only one coarse level: exact solve
            L = obj.nLevels;
            if L == 1
                Cx = obj.A \ x;
                algEstimator = sqrt(dot(Cx, obj.A * Cx, Dim.Vector));
                return
            end

            % descending cascade: no smoothing
            res{L} = x;
            for k = L:-1:2
                res{k-1} = obj.smoother.restrict(res{k}, k);
            end

            % exact solve on coarsest level to compute accumulative lifting of x
            Cx = obj.projectedMatrix{1} \ res{1};
            eta2 = dot(Cx, obj.projectedMatrix{1} * Cx, Dim.Vector);

            % ascending cascade: (local) smoothing and compute error estimator
            for k = 2:(L-1)
                Cx = obj.smoother.prolongate(Cx, k);
                updatedResidual = res{k} - obj.projectedMatrix{k} * Cx;
                rho = obj.smoother.smooth(updatedResidual, k);
                [etaUpdate2, xUpdate] = MultigridSolver.computeOptimalUpdate(obj.projectedMatrix{k}, updatedResidual, rho);
                Cx = Cx + xUpdate;
                eta2 = eta2 + etaUpdate2;
            end

            % final smoothing and residual update (with non-projected matrix)
            Cx = obj.smoother.prolongate(Cx, L);
            updatedResidual = res{L} - obj.A * Cx;
            rho = obj.smoother.smooth(updatedResidual, L);
            [etaUpdate2, xUpdate] = MultigridSolver.computeOptimalUpdate(obj.A, updatedResidual, rho);
            Cx = Cx + xUpdate;
            algEstimator = sqrt(eta2 + etaUpdate2);
        end
    end
end