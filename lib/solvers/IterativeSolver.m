% IterativeSolver (abstract handle class) Interface for iterative solvers.
%
% solver.setup(A, b, [x0]) sets up solver for the linear system A*x=b with
%   initial guess x0 (0 per default). The right-hand side b can have
%   multiple columns.
%
% solver.solve() performs an automatic loop of solver.step() until
%   convergence in each column of the right-hand side is reached, checking
%   the stopping criterion after every step.
%
% solver.step() performs one step of the iterative solver on active
%   components of the right-hand side and automatoically keeps track of
%   iteration numbers.
%
% tf = solver.applyStoppingCriterion() checks stopping criterion and set
%   converged components to inactive.
%
% Abstract methods to implement:
%
% solver.computeUpdate() computes the update of active components of the
%   solution.
%
% isConverged(solver) returns state of convergence of the solver for each
%   column of the right-hand side.


classdef IterativeSolver < handle
    %% properties
    properties (Access=public)
        maxIter (1,:) double {mustBePositive}
        tol (1,:) double {mustBePositive}
    end
    
    properties (SetAccess=protected, GetAccess=public)
        A
        b
        x
        iterationCount
        activeComponents
    end
    
    methods (Abstract, Access=public)
        isConverged(obj)
    end
    
    methods (Abstract, Access=protected)
        computeUpdate(obj)
    end
    
    %% generic routines for iterative solvers
    % overwrite to include additional steps
    methods (Access=public)
        function obj = IterativeSolver()
            obj.maxIter = 100;
            obj.tol = 1e-8;
        end
        
        function setupLinearSystem(obj, A, b, x0)
            arguments
                obj
                A (:,:) double {IterativeSolver.mustBeQuadratic}
                b (:,:) double {IterativeSolver.mustBeCompatible(A, b)}
                x0 (:,:) double {IterativeSolver.mustBeCompatible(A, x0)} = zeros(size(b));
            end
            
            assert(isequal(size(b,2), size(x0,2)), ...
                'Right-hand size and initial guess must be of compatible size.')
            
            obj.A = A;
            obj.b = b;
            obj.x = x0;
            obj.iterationCount = zeros(1, size(b,2));
            obj.activeComponents = true(1, size(b,2));
        end
        
        function x = solve(obj)
            while ~all(applyStoppingCriterion(obj))
                obj.step()
            end
            
            if any(obj.iterationCount >= obj.maxIter)
                warning('Maximal number of iterations exhausted! maxIter = %d', obj.maxIter)
            end
            
            x = obj.x;
        end
        
        function step(obj)
            obj.computeUpdate();
            obj.iterationCount(obj.activeComponents) = ...
                obj.iterationCount(obj.activeComponents) + 1;
        end
        
        function tf = applyStoppingCriterion(obj)
            tf = isConverged(obj);
            obj.activeComponents = obj.activeComponents & ~tf;
        end
    end
    
    %% validation functions to use within this class
    methods (Static, Access=protected)
        function mustBeQuadratic(A)
            if ~(size(A,1) == size(A,2))
                eidType = 'Solver:matrixNotQuadratic';
                msgType = 'Matrix must be quadratic.';
                throwAsCaller(MException(eidType,msgType))
            end
        end
        
        function mustBeCompatible(A, b)
            if ~(size(A,2) == size(b,1))
                eidType = 'Solver:dataNotCompatible';
                msgType = 'Vector must have size compatible with matrix.';
                throwAsCaller(MException(eidType,msgType))
            end
        end
    end
end