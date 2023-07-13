% DiagonalJacobi (subclass of Preconditioner) diagonal preconditioner.
%
% See also: Preconditioner, PcgSolver

classdef DiagonalJacobi < Preconditioner
    %% properties
    properties (Access=private)
        C
    end
    
    %% methods
    methods (Access=public)
        % preconditioner action: inverse of diagonal
        function setup(obj, A)
            obj.C = full(diag(A)).^(-1);
        end
           
        function Cx = apply(obj, x)
            Cx = obj.C .* x;
        end
    end
end