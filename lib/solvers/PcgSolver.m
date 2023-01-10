% PcgSolver (abstract subclass of IterativeSolver) Solves linear equations
%   iteratively by a preconditioned conjugate gradient (PCG) method.
%
% See also: IterativeSolver

classdef PcgSolver < IterativeSolver
    %% properties
    properties (SetAccess=protected, GetAccess=public)
        residualCNorm
    end
    
    properties (Access=protected)
        residual
        Cresidual
        searchDirection
        normb
    end
    
    %% extend superclass methods
    methods (Access=public)
        function setupSystemMatrix(obj, A)
            setupSystemMatrix@IterativeSolver(obj, A);
        end
        
        function setupRhs(obj, b, varargin)
            setupRhs@IterativeSolver(obj, b, varargin{:});
            
            % initialize residual & search direction
            obj.residual = b - obj.A*obj.x;
            obj.Cresidual = obj.preconditionAction(obj.residual);
            obj.searchDirection = obj.Cresidual;
            obj.residualCNorm = sum(obj.residual.*obj.Cresidual, 1);
            obj.normb = sqrt(sum(b.^2, 1));
        end
    end
    
    %% implement abstract superclass methods
    methods (Access=public)
        function tf = isConverged(obj)
            tf = ((obj.iterationCount >= obj.maxIter) ...
                        | (sqrt(obj.residualCNorm) < obj.tol) ...
                        | (sqrt(obj.residualCNorm)./obj.normb < obj.tol));
        end
    end
    
    methods (Access=protected)
        function computeUpdate(obj)
            % only iterate on active components
            idx = find(obj.activeComponents);
            
            % update solution
            AsearchDirection = obj.A * obj.searchDirection(:,idx);
            if sum(obj.searchDirection(:,idx).*AsearchDirection, 1) < eps
                alpha = 1;
            else
                alpha = obj.residualCNorm(:,idx) ./ sum(obj.searchDirection(:,idx).*AsearchDirection, 1);
            end
            obj.x(:,idx) = obj.x(:,idx) + alpha .* obj.searchDirection(:,idx);
            
            % update residual
            obj.residual(:,idx) = obj.residual(:,idx) - alpha .* AsearchDirection;
            obj.Cresidual(:,idx) = obj.preconditionAction(obj.residual(:,idx));
            residualCNormOld = obj.residualCNorm(:,idx);
            obj.residualCNorm(:,idx) = sum(obj.residual(:,idx).*obj.Cresidual(:,idx), 1);
            
            % update search direction
            beta = obj.residualCNorm(:,idx) ./ residualCNormOld;
            obj.searchDirection(:,idx) = obj.Cresidual(:,idx) + beta .* obj.searchDirection(:,idx);
        end
    end
    
    %% abstract preconditioning
    methods (Abstract, Access=public)
        preconditionAction(obj, x)
    end
end