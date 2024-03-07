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
            assert(isequal(size(x, 1), size(obj.A, 1)), ... 
                'Setup for multilevel iteration not complete!')
            
            % if there is only one coarse level: exact solve
            L = obj.nLevels;
            if L == 1
                Cx = obj.A \ x;
                algEstimator = sqrt(dot(Cx, obj.A * Cx, Dim.Vector));
                return
            end

            % sequential multiplicative V-cycle update
            algEstimator = 0;
            v = obj.x;
            for ell = 1:L
                xUpdate = obj.Vcycle(x, ell);
                % A priori step size
                if ell == L-1
                    rho = 1;
                else
                    rho = 0.8;
                end
                v = v + rho * xUpdate;
                % Compute optimal step size
                % updatedResidual = x - obj.A * Cx;
                % rho = obj.smoother.smooth(updatedResidual, ell+1);
                % [etaUpdate2, xUpdate] = MultigridSolver.computeOptimalUpdate(obj.A, updatedResidual, rho);
                % Cx = Cx + xUpdate;
                % Update estimator
                % algEstimator = algEstimator + etaUpdate2;
                % Update residual
                x = x - obj.A * rho * xUpdate;
            end
            Cx = v;
        end

        function Cx = Vcycle(obj, x, ell) 
            % one step of the Vcycle
            
            L = obj.nLevels;

            if ell == L
                Cx = obj.smoother.smooth(x, ell);
                return
            end

            % descending cascade: no smoothing
            res = x;
            for k = L:-1:ell+2
                res = obj.smoother.restrict(res, k);
            end
            res = obj.smoother.averageRestrict(res, ell+1);

            % smoothing step on level ell
            if ell == 1
                Cx = obj.projectedMatrix{1} \ res;
            else
                Cx = obj.smoother.smooth(res, ell);
            end
            Cx = obj.smoother.average(Cx, ell+1);

            % ascending cascade: no smoothing
            for k = ell+2:L
                Cx = obj.smoother.prolongate(Cx, k);
            end
        end
    end
end