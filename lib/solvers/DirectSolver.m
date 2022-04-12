% DirectSolver (subclass of Solver) Solves via backslash operator.

classdef DirectSolver < Solver
    methods
        function u = solve(~, A, b)
            arguments
                ~
                A (:,:) double {Solver.mustBeQuadratic} 
                b (:,1) double {Solver.mustBeCompatible(A, b)}
            end
            u = A \ b;
        end
    end
end