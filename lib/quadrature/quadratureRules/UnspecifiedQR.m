% UnspecifiedQR Null object as placeholder for not-set quadrature rule.

classdef UnspecifiedQR < QuadratureRule
    %% methods
    methods (Access=public)
        function obj = UnspecifiedQR()
            obj = obj@QuadratureRule(Barycentric1D([0; 1]), 1);
            obj.bary = [];
            obj.weights = [];
            obj.nNodes = 0;
            obj.dim = NaN;
        end
    end
end