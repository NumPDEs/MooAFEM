classdef TestIntegration < matlab.unittest.TestCase
    
properties (TestParameter)
    order = struct('low', 1, 'mediumLow', 5, 'mediumHigh', 12, 'high', 23);
end

methods (Test)
    function integrationOfConstantScalesProperlyWithMesh(testCase, order)
        mesh = Mesh.loadFromGeometry('unitsquare');
        qr = QuadratureRule.ofOrder(order);
        a = Constant(mesh, pi);
        actual = ones(1,4);
        for i = 1:4
            mesh.refineUniform(1);
            actual(i) = sum(integrateElement(a, qr), Dim.Elements);
        end
        testCase.verifyEqual(actual, [1,1,1,1]*pi, 'RelTol', eps*1e2)
    end
    
    function nonscalarFunctionIsTreatedProperly(testCase, order)
        mesh = Mesh.loadFromGeometry('unitsquare');
        qr = QuadratureRule.ofOrder(order);
        f = MeshFunction(mesh, @(x) x);
        actual = ones(2,4);
        for i = 1:4
            mesh.refineUniform(1);
            actual(:,i) = sum(integrateElement(f, qr), Dim.Elements);
        end
        testCase.verifyEqual(actual, ones(2,4)/2, 'RelTol', eps*1e2)
    end
    
    function integrationOverEdgesScalesWithEdgeLengths(testCase, order)
        mesh = Mesh.loadFromGeometry('unitsquare');
        qr = QuadratureRule.ofOrder(order, '1D');
        A = Constant(mesh, [1;2;3;4]);
        integral = sum(integrateEdge(A, qr), Dim.Elements);
        trafo = getAffineTransformation(mesh);
        testCase.verifyEqual(integral, [1;2;3;4]*sum(trafo.ds), 'RelTol', sqrt(eps))
    end
    
    function edgeJumpIntegralIsBoundaryIntegralForConstants(testCase, order)
        mesh = Mesh.loadFromGeometry('unitsquare');
        qr = QuadratureRule.ofOrder(order, '1D');
        a = Constant(mesh, pi);
        actual = ones(1,4);
        for i = 1:4
            mesh.refineUniform(1);
            actual(i) = sum(abs(integrateJump(a, qr)), Dim.Elements);
        end
        testCase.verifyEqual(actual, [4,4,4,4]*pi, 'RelTol', eps*1e2)
    end
    
    function normalJumpPostprocessingGivesCorrectResults(testCase, order)
        mesh = Mesh.loadFromGeometry('unitsquare');
        qr = QuadratureRule.ofOrder(order, '1D');
        f = MeshFunction(mesh, @(x) x);
        integral = sum(integrateNormalJump(f, qr, @(j) abs(j), {}, ':'), Dim.Elements);
        trafo = getAffineTransformation(mesh);
        testCase.verifyEqual(integral, 2, 'RelTol', sqrt(eps))
    end
end
    
end
