classdef TestQuadrature < matlab.unittest.TestCase
    
properties (TestParameter)
    p = struct('low', 1, 'mediumLow', 5, 'mediumHigh', 13, 'high', 22);
end

methods (Test)
    function barycentricReshapingHasCorrectDimension(testCase)
        nodes = [ones(3,1)/3, eye(3)];
        bary = Barycentric2D(nodes);
        testCase.verifyEqual(size(bary.arangeAlongDim(3)), [3, 1, 4])
    end

    function integrityCheckThrowsErrorIfSumNotOne(testCase)
        testCase.verifyError(@() Barycentric2D([1;1;1]), 'Barycentric:notBarycentric')
    end

    function quadrature1DExactForOrderP(testCase, p)
        qr = QuadratureRule.ofOrder(p, '1D');
        x = qr.bary.coordinates(1,:);
        integral = (x.^p) * qr.weights';
        testCase.verifyEqual(integral, 1/(p+1), 'RelTol', eps*1e2)
    end

    function quadrature2DExactForMonomials(testCase, p)
        qr = QuadratureRule.ofOrder(p);
        zxy = qr.bary.coordinates;
        integral = (zxy.^p) * qr.weights' / 2;
        testCase.verifyEqual(integral, 1/((p+1)*(p+2))*[1;1;1], 'RelTol', eps*1e2)
    end

    function quadrature2DExactForMultinomials(testCase, p)
        qr = QuadratureRule.ofOrder(p);
        x = qr.bary.coordinates(2,:);
        y = qr.bary.coordinates(3,:);
        a = ceil(p/2);
        b = floor(p/2);
        integral = (x.^a .* y.^b) * qr.weights' / 2;
        expected = factorial(a) * factorial(b) / factorial(p + 2);
        testCase.verifyEqual(integral, expected, 'RelTol', eps*1e2)
    end
end

end
