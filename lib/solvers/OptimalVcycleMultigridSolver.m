% OptimalVcycleMultigridSolver (subclass of IterativeSolver) Solves linear
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
% See also IterativeSolver, MultiLevelSmoother, DiagonalJacobi,
%   LowOrderBlockJacobi, HighOrderBlockJacobi

classdef OptimalVcycleMultigridSolver < IterativeSolver
    %% properties 
    properties (SetAccess=protected, GetAccess=public)
        residualCNorm
        algEstimator
        residual
        Cresidual
    end

    properties (Dependent, Access=public)
        nLevels
    end
    
    properties (Access=protected)
        smoother
        projectedMatrix
        lhsNorm
    end
    
    methods
        function obj = OptimalVcycleMultigridSolver(smoother)
            obj = obj@IterativeSolver();
            obj.smoother = smoother;
        end

        function n = get.nLevels(obj)
            n = obj.smoother.nLevels;
        end

        function setupSystemMatrix(obj, A, varargin)
            % override superclass method, because ALL matrices are needed
            setupSystemMatrix@IterativeSolver(obj, A, varargin{:});

            projectedA = obj.smoother.setup(A, varargin{:});
            obj.projectedMatrix{obj.nLevels} = projectedA;
        end

        function setupRhs(obj, b, varargin)
            setupRhs@IterativeSolver(obj, b, varargin{:});
            % initialize residual & hierarchy
            obj.residual = b - obj.A * obj.x;
            [obj.Cresidual, obj.algEstimator] = obj.Vcycle(obj.residual);
            obj.residualCNorm = sqrt(dot(obj.residual, obj.Cresidual, Dim.Vector));
            obj.lhsNorm = vecnorm(b, 2, Dim.Vector);
        end
        
        % stopping criterion
        function tf = isConverged(obj)
            tf = ((obj.iterationCount >= obj.maxIter) ...
                | (obj.residualCNorm < obj.tol) ...
                | (obj.residualCNorm ./ obj.lhsNorm < obj.tol)); 
        end
    end
        
    methods (Access=protected)   
        % one step
        function computeUpdate(obj)
            % only iterate on active components
            idx = find(obj.activeComponents);
            
            % update solution & residual
            obj.x(:,idx) = obj.x(:,idx) + obj.Cresidual(:,idx);
            obj.residual(:,idx) = obj.residual(:,idx) - obj.A * obj.Cresidual(:,idx);

            % residual & update error correction
            [obj.Cresidual(:,idx), obj.algEstimator(idx)] = obj.Vcycle(obj.residual(:,idx));
            obj.residualCNorm(:,idx) = sqrt(dot(obj.residual(:,idx), obj.Cresidual(:,idx), Dim.Vector));
        end

        function [Cx, algEstimator] = Vcycle(obj, x)
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
                [etaUpdate2, xUpdate] = computeOptimalUpdate(obj.projectedMatrix{k}, updatedResidual, rho);
                Cx = Cx + xUpdate;
                eta2 = eta2 + etaUpdate2;
            end

            % final smoothing and residual update (with non-projected matrix)
            Cx = obj.smoother.prolongate(Cx, L);
            updatedResidual = res{L} - obj.A * Cx;
            rho = obj.smoother.smooth(updatedResidual, L);
            [etaUpdate2, xUpdate] = computeOptimalUpdate(obj.A, updatedResidual, rho);
            Cx = Cx + xUpdate;
            algEstimator = sqrt(eta2 + etaUpdate2);
        end
    end
end

%% auxiliary functions
% error correction with optimal stepsize 
function [etaUpdate2, resUpdate] = computeOptimalUpdate(A, res, rho)
    rhoArho = dot(rho, A*rho, Dim.Vector);
    if max(abs(rho)) < eps
        lambda = 1; 
    else
        lambda = dot(res, rho, Dim.Vector) ./ rhoArho;
    end
    resUpdate = lambda .* rho;
    etaUpdate2 = rhoArho .* lambda.^2;
    
    if any(lambda > 3)
       warning('MG step-sizes no longer bound by d+1. Optimality of step size cannot be guaranteed!')
    end
end
