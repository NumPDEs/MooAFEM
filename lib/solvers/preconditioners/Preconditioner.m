% Preconditioner (abstract handle class) Interface for preconditioners.
%
% pc.setup(A, ...) hook to set up the preconditioner for the matrix A.
%   Default implementation does nothing.
%
% Abstract methods to implement:
%
% Cx = pc.apply(x) applies the preconditioner to the vector x.
%
% See also: PcgSolver


classdef Preconditioner < handle
    methods (Access=public)
        function setup(~, ~)
            % Default implementation does nothing.
        end
    end
    
    methods (Abstract, Access=public)
        apply(obj, x)
    end
end