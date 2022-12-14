% DirectSolver (subclass of IterativeSolver) Solves linear equations
%   directly. This is mainly intended for testing purposes!
%
% See also: IterativeSolver

classdef DirectSolver < IterativeSolver
    %% implement abstract superclass methods
    methods (Access=public)
        function tf = isConverged(obj)
            % converges in the first step
            tf = (obj.iterationCount > 0);
        end
    end
    
    methods (Access=protected)
        function computeUpdate(obj)
            % direct solve
            obj.x = obj.A \ obj.b;
        end
    end
end