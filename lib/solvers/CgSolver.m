% CgSolver (subclass of PCGSolver) Solves linear equations iteratively using
%   the CG method without preconditioner.

classdef CgSolver < PcgSolver
    %% methods
    methods (Access=public)
        % preconditioner = identity
        function Cx = preconditionAction(obj, x)
            Cx = x;
        end
    end
end