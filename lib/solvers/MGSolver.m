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

classdef MGSolver < Solver
    %% properties
    properties (Access=public)
        maxIter (1,:) double {mustBePositive}
        tol (1,:) double {mustBePositive}
    end
    
    properties (SetAccess=protected, GetAccess=public)
        A
        b
        x
        residualCNorm
        iterationCount
        algEst2
    end
    
    properties (Access=protected)
        residual
        Cresidual
        searchDirection
        normb
        activeComponents
    end
    
    %% methods
    methods (Access=public)
        function obj = MGSolver(varargin)
            obj.maxIter = 100;
            obj.tol = 1e-8;
        end
        
        % complete iteration
        function x = solve(obj, A, b, varargin)
            % adaptive loop
            obj.setup(A, b, varargin{:});
            while ~all(isConverged(obj))
                obj.step()
            end
            
            if any(obj.iterationCount >= obj.maxIter)
                warning('MG iterations exhausted! maxIter = %d', obj.maxIter)
            end
            
            x = obj.x;
        end
        
        % initialize iteration
        function setup(obj, A, b, x0, c)
            arguments
                obj
                A (:,:) double {Solver.mustBeQuadratic}
                b (:,:) double {Solver.mustBeCompatible(A, b)}
                x0 (:,:) double {Solver.mustBeCompatible(A, x0)} = zeros(size(b));
            end
            
            arguments (Repeating)
                c % catching possible additional arguments
            end
            
            assert(isequal(size(b,2), size(x0,2)), 'Right-hand size and initial guess must be of compatible size.')
            
            obj.A = A;
            obj.b = b;
            obj.x = x0;
            obj.residual = b - obj.A*obj.x;
            [obj.Cresidual obj.algEst2] = obj.Vcycle(obj.residual);
            obj.searchDirection = obj.Cresidual;
            obj.residualCNorm = sum(obj.residual.*obj.Cresidual, 1);
            obj.iterationCount = zeros(1, size(b,2));
            obj.normb = sqrt(sum(b.^2, 1));
            obj.activeComponents = true(1, size(b,2));
        end
        
        
        % one step
        function step(obj)
            % only iterate on active components
            idx = ~isConverged(obj);
            stillActive = idx(obj.activeComponents);
            obj.activeComponents = obj.activeComponents & idx;

            % residual & update error correction
            obj.residual = obj.residual(:,stillActive); 
            [obj.Cresidual obj.algEst2] = obj.Vcycle(obj.residual);
            obj.residualCNorm = sum(obj.residual.*obj.Cresidual, 1);
            
            % update search direction & solution
            obj.searchDirection = obj.Cresidual;
            obj.searchDirection = obj.searchDirection(:,stillActive);
            obj.x(:,idx) = obj.x(:,idx) + obj.searchDirection;

            % update residual
            AsearchDirection = obj.A * obj.searchDirection;
            obj.residual = obj.residual(:,stillActive) -  AsearchDirection;
            
            
            obj.iterationCount(:,idx) = obj.iterationCount(:,idx)+1;



%             % only iterate on active components
%             idx = ~isConverged(obj);
%             stillActive = idx(obj.activeComponents);
%             obj.activeComponents = obj.activeComponents & idx;
%             
%             % update solution
%             obj.searchDirection = obj.searchDirection(:,stillActive);
%             AsearchDirection = obj.A * obj.searchDirection;
%             obj.x(:,idx) = obj.x(:,idx) + obj.searchDirection;
%             
%             % update residual
%             obj.residual = obj.residual(:,stillActive) -  AsearchDirection;
%             [obj.Cresidual obj.algEst2] = obj.Vcycle(obj.residual);
%             obj.residualCNorm = sum(obj.residual.*obj.Cresidual, 1);
%             
%             % update search direction
%             obj.searchDirection = obj.Cresidual;
%             obj.iterationCount(:,idx) = obj.iterationCount(:,idx)+1;
        end
        
        % stopping criterion
        function tf = isConverged(obj)
                        tf = ( (sqrt(obj.residualCNorm)./obj.normb < obj.tol));        
        end
    end
    
    methods (Abstract,Access=public)
        Vcycle(obj, x)
    end
end








