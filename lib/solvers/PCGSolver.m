% PCGSolver (abstract subclass of Solver) Solves linear equations
%   iteratively, using a preconditioner.
%
% 

classdef PCGSolver < Solver
    %% properties
    properties (SetAccess=protected, GetAccess=public)
        A
        b
        x
        maxIter
        tol
        residualCNorm
        iterationCount
    end
    
    properties (Access=protected)
        residual
        Cresidual
        searchDirection
        normb
    end
    
    %% methods
    methods
        function obj = PCGSolver(varargin)
            obj = obj@Solver(varargin);
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
                warning('(P)CG iterations exhausted! maxIter = %d', obj.maxIter)
            end
            
            x = obj.x;
        end
        
        % initialize iteration
        function setup(obj, A, b, x0)
            arguments
                obj
                A (:,:) double {Solver.mustBeQuadratic}
                b (:,:) double {Solver.mustBeCompatible(A, b)}
                x0 (:,:) double {Solver.mustBeCompatible(A, x0)} = zeros(size(b));
            end
            
            assert(isequal(size(b,2), size(x0,2)), 'Right-hand size and initial guess must be of compatible size.')
            
            obj.A = A;
            obj.b = b;
            obj.x = x0;
            obj.residual = b - obj.A*obj.x;
            obj.Cresidual = obj.preconditionAction(obj.residual);
            obj.searchDirection = obj.Cresidual;
            obj.residualCNorm = sum(obj.residual.*obj.Cresidual, 1);
            obj.iterationCount = zeros(1, size(b,2));
            obj.normb = sqrt(sum(b.^2, 1));
        end
        
        
        % one step
        function step(obj)
            % TODO: handle 0-solution in a correct way!
            AsearchDirection = obj.A * obj.searchDirection;
            alpha = obj.residualCNorm ./ sum(obj.searchDirection.*AsearchDirection, 1);
            obj.x = obj.x + alpha .* obj.searchDirection;
            obj.residual = obj.residual - alpha .* AsearchDirection;
            obj.Cresidual = obj.preconditionAction(obj.residual);
            residualCNormOld = obj.residualCNorm;
            obj.residualCNorm = sum(obj.residual.*obj.Cresidual, 1);
            beta = obj.residualCNorm ./ residualCNormOld;
            obj.searchDirection = obj.Cresidual + beta .* obj.searchDirection;
            obj.iterationCount = obj.iterationCount+1;
        end
        
        % stopping criterion
        function tf = isConverged(obj)
            tf = (obj.iterationCount > 0);
            tf = tf & ((obj.iterationCount >= obj.maxIter) ...
                        | (sqrt(obj.residualCNorm)./obj.normb < obj.tol));
        end
    end
    
    methods(Abstract)
        preconditionAction(obj, x)
    end
end