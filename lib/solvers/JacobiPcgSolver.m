% JacobiPcgSolver (subclass of PcgSolver) inear equations iteratively using
%   the CG method with diagonal preconditioner.

classdef JacobiPcgSolver < PcgSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        C (:,1) double
    end
    
    %% methods
    methods (Access=public)
        % preconditioner action: inverse of diagonal
        function setupLinearSystem(obj, A, b, x0)
            obj.C = 1./full(diag(A));
            setupLinearSystem@PcgSolver(obj, A, b, x0)
        end
        
        function Cx = preconditionAction(obj, x)
            Cx = obj.C .* x;
        end
    end
end