% BernsteinBezierPoly Provide higher order polynomials
%
%   The Bernstein-Bezier polynomials of order p are sums of the basis polynomials
%       C * lambda_1^i lambda_2^j lambda_3^k,
%   where 0 < i+j+k <= p. This class provides efficient evaluation of these on
%   and their derivatives the reference triangle.

classdef BernsteinBezierPoly
    %% properties
    properties (Access='protected')
        constants
        exponents
        nFun
    end
    
    %% methods
    methods (Access='public')
        function obj = setExponents(obj, exponents)
            arguments
                obj
                exponents (3,:) double {mustBeNonnegative}
            end
            
            obj.exponents = exponents;
            n = sum(exponents, 1);
            obj.constants = factorial(n)./prod(factorial(exponents), 1);
        end
        
        function val = evalShapeFunctions(obj, bary)
            arguments
                obj
                bary Barycentric2D
            end
            val = pagetranspose(obj.constants .* ...
                                prod(bary.arangeAlongDim(3).^obj.exponents, 1));
        end
        
        function val = evalShapeFunctionDerivs(obj, bary)
            arguments
                obj
                bary Barycentric2D
            end
            d1 = obj.exponents(1,:).*obj.constants .* ...
                prod(bary.arangeAlongDim(3).^(max(obj.exponents-[1;0;0], 0)), 1);
            d2 = obj.exponents(2,:).*obj.constants .* ...
                prod(bary.arangeAlongDim(3).^(max(obj.exponents-[0;1;0], 0)), 1);
            d3 = obj.exponents(3,:).*obj.constants .* ...
                prod(bary.arangeAlongDim(3).^(max(obj.exponents-[0;0;1], 0)), 1);
            
            val = pagetranspose(cat(2, d2-d1, d3-d1));
        end
        
        function val = evalShapeFunctionHessian(obj, bary)
            arguments
                obj
                bary Barycentric2D
            end
            d11 = obj.exponents(1,:).*(obj.exponents(1,:)-1).*obj.constants .* ...
                prod(bary.arangeAlongDim(3).^(max(obj.exponents-[2;0;0], 0)), 1);
            d22 = obj.exponents(2,:).*(obj.exponents(2,:)-1).*obj.constants .* ...
                prod(bary.arangeAlongDim(3).^(max(obj.exponents-[0;2;0], 0)), 1);
            d33 = obj.exponents(3,:).*(obj.exponents(3,:)-1).*obj.constants .* ...
                prod(bary.arangeAlongDim(3).^(max(obj.exponents-[0;0;2], 0)), 1);
            d12 = obj.exponents(1,:).*obj.exponents(2,:).*obj.constants .* ...
                prod(bary.arangeAlongDim(3).^(max(obj.exponents-[1;1;0], 0)), 1);
            d13 = obj.exponents(1,:).*obj.exponents(3,:).*obj.constants .* ...
                prod(bary.arangeAlongDim(3).^(max(obj.exponents-[1;0;1], 0)), 1);
            d23 = obj.exponents(2,:).*obj.exponents(3,:).*obj.constants .* ...
                prod(bary.arangeAlongDim(3).^(max(obj.exponents-[0;1;1], 0)), 1);
            
            dxy = d11+d23-d12-d13;
            val = pagetranspose(cat(2, d11+d22-2*d12, dxy, dxy, d11+d33-2*d13));
        end
    end
end
