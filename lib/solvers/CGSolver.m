% GSolver (subclass of PCGSolver) Solves linear equations iteratively using
%   the CG method without preconditioner.

classdef CGSolver < PCGSolver    
    %% methods
    methods
        % preconditioner = identity
        function Cx = preconditionAction(obj, x)
            Cx = x;
        end
    end
end