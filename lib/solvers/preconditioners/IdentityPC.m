% IdentityPC (subclass of Preconditioner) identity preconditioner (=no-op).
%
% See also: Preconditioner


classdef IdentityPC < Preconditioner
    methods (Access=public)
        function Cx = apply(~, x)
            Cx = x;
        end
    end
end