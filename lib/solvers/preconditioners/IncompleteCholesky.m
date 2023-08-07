% IncompleteCholesky (subclass of Preconditioner) incomplete Cholesky
%   decomposition preconditioner.
%
% See also: Preconditioner, PcgSolver

classdef IncompleteCholesky < Preconditioner
    %% properties
    properties (Access=private)
        C
    end
    
    %% methods
    methods (Access=public)
        function setup(obj, A)
            obj.C = ichol(A);
        end

        function Cx = apply(obj, x)
            Cx = obj.C' \ (obj.C \ x);
        end
    end
end