% MGSolver (abstract subclass of IterativeSolver) Solves linear equations
%   iteratively, using a Vcycle geometric multigrid solver.
%
% Abstract methods to be implemented by subclasses:
%
%   solver.Vcycle(x) performs one Vcycle iteration.
%
% See also IterativeSolver

classdef MGSolver < IterativeSolver
    %% properties    
    properties (SetAccess=protected, GetAccess=public)
        residualCNorm
        algEstimator
    end
    
    properties (Access=protected)
        residual
        Cresidual
        normb
    end
    
    %% methods
    methods (Access=public)        
        % initialize iteration
        function setupSystemMatrix(obj, A, varargin)
            setupSystemMatrix@IterativeSolver(obj, A);
        end

        function setupRhs(obj, b, varargin)
            setupRhs@IterativeSolver(obj, b, varargin{:});
            % initialize residual & hierarchy
            obj.residual = b - obj.A*obj.x;
            [obj.Cresidual, obj.algEstimator] = obj.Vcycle(obj.residual);
            obj.residualCNorm = sum(obj.residual.*obj.Cresidual, 1);
            obj.normb = sqrt(sum(b.^2, 1));
        end
        
        % stopping criterion
        function tf = isConverged(obj)
            tf = ((obj.iterationCount >= obj.maxIter) ...
                        | (sqrt(obj.residualCNorm) < obj.tol) ...
                        | (sqrt(obj.residualCNorm)./obj.normb < obj.tol)); 
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
            obj.residualCNorm(:,idx) = sum(obj.residual(:,idx).*obj.Cresidual(:,idx), 1);
        end
    end
    
    methods (Abstract,Access=public)
        Vcycle(obj, x)
    end
end








