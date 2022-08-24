% MGSolver (abstract subclass of Solver) Solves linear equations
%   iteratively, using a Vcycle.
%
% solver.setup(A, b, [x0]) sets up solver for the linear system A*x=b with
%   initial guess x0 (0 per default). The right-hand side b can have
%   multiple columns.
%
% solver.step() performs one MG step.
%
% isConverged(solver) returns state of convergence of the solver for each
%   column of the right-hand side.
%
% solver.solve(A, b, [x0]) performs an automatic loop of solver.step()
%   until convergence in each column of the right-hand side is reached.

classdef MGSolver < IterativeSolver
    %% properties    
    properties (SetAccess=protected, GetAccess=public)
        residualCNorm
        algEst2
    end
    
    properties (Access=protected)
        residual
        Cresidual
        normb
    end
    
    %% methods
    methods (Access=public)        
        % initialize iteration
        function setupLinearSystem(obj, A, b, x0)
            setupLinearSystem@IterativeSolver(obj, A, b, x0);
            
            % initialize residual & hierarchy
            obj.residual = b - obj.A*obj.x;
            [obj.Cresidual, obj.algEst2] = obj.Vcycle(obj.residual);
            obj.residualCNorm = sum(obj.residual.*obj.Cresidual, 1);
            obj.normb = sqrt(sum(b.^2, 1));
        end
        
        % stopping criterion
        function tf = isConverged(obj)
            tf = ((obj.iterationCount >= obj.maxIter) ...
                        | (sqrt(obj.residualCNorm)./obj.normb < obj.tol));  
        end
    end
        
    methods (Access=protected)   
        % one step
        function computeUpdate(obj)
            % only iterate on active components
            idx = find(obj.activeComponents);

            % residual & update error correction
            [obj.Cresidual(:,idx), obj.algEst2] = obj.Vcycle(obj.residual(:,idx));
            obj.residualCNorm(:,idx) = sum(obj.residual(:,idx).*obj.Cresidual(:,idx), 1);
            
            % update search direction & solution
            obj.x(:,idx) = obj.x(:,idx) + obj.Cresidual(:,idx);

            % update residual
            obj.residual(:,idx) = obj.residual(:,idx) - obj.A * obj.Cresidual(:,idx);
        end
    end
    
    methods (Abstract,Access=public)
        Vcycle(obj, x)
    end
end








