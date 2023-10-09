% NoPreconditioner (subclass of Preconditioner) identity preconditioner (=no-op).
%
% See also: Preconditioner, PcgSolver


classdef NoPreconditioner < Preconditioner
    methods (Access=public)
        function Cx = apply(~, x)
            Cx = x;
        end
    end
end