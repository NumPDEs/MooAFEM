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
% See also IterativeSolver, MultiLevelSmoother, DiagonalJacobi,
%   LowOrderBlockJacobi, HighOrderBlockJacobi

classdef MultigridSolver < IterativeSolver
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
        function obj = MultigridSolver(smoother)
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

        function updateSystemMatrix(obj, A, varargin)
            % override superclass method, because ALL matrices are needed
            updateSystemMatrix@IterativeSolver(obj, A, varargin{:});

            projectedA = obj.smoother.update(A, varargin{:});
            obj.projectedMatrix{obj.nLevels} = projectedA;
        end

        function setupRhs(obj, b, varargin)
            setupRhs@IterativeSolver(obj, b, varargin{:});
            % initialize residual & hierarchy
            obj.residual = b - obj.A * obj.x;
            obj.lhsNorm = vecnorm(b, 2, Dim.Vector);
            obj.Cresidual = zeros(size(obj.residual));
            obj.algEstimator = Inf;
            obj.residualCNorm = Inf;
        end
        
        % stopping criterion
        function tf = isConverged(obj)
            tf = ((obj.iterationCount >= obj.maxIter) ...
                | (obj.residualCNorm < obj.tol) ...
                | (obj.residualCNorm ./ obj.lhsNorm < obj.tol)); 
        end
    end
        
    methods (Access=protected)
        function computeUpdate(obj)
            % one step

            % only iterate on active components
            idx = find(obj.activeComponents);

            % residual & update error correction
            [obj.Cresidual(:,idx), obj.algEstimator(idx)] = obj.cycle(obj.residual(:,idx));
            obj.residualCNorm(:,idx) = sqrt(dot(obj.residual(:,idx), obj.Cresidual(:,idx), Dim.Vector));   

            % update solution & residual
            obj.x(:,idx) = obj.x(:,idx) + obj.Cresidual(:,idx);
            obj.residual(:,idx) = obj.residual(:,idx) - obj.A * obj.Cresidual(:,idx);
        end

        % carry out one cycle
        [Cx, algEstimator] = cycle(obj, x)
    end

    methods (Static, Access=protected)
        function [etaUpdate2, resUpdate] = computeOptimalUpdate(A, res, rho)
            % error correction with optimal stepsize 
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
    end
end