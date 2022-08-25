classdef TestAssembly < matlab.unittest.TestCase
    
properties
    mesh
	fes
    blf
    lf
end

methods (TestMethodSetup)
    function initializeForms(testCase)
        % expected values below are computed by hand for unit triangle!
        testCase.mesh = Mesh.loadFromGeometry('unittriangle');
        testCase.fes = FeSpace(testCase.mesh, LowestOrderH1Fe());
        testCase.blf = BilinearForm();
        testCase.lf = LinearForm();
    end
end

methods (Test)
    function assembleWithoutCoefficientThrowsError(testCase)
        testCase.blf.a = [];
        testCase.verifyError(@() assemble(testCase.blf, testCase.fes), 'BilinearForm:NoCoefficients')
        testCase.verifyError(@() assemble(testCase.lf, testCase.fes), 'LinearForm:NoCoefficients')
    end
    
    function validateStiffnessScalar(testCase)
        testCase.blf.a = Constant(testCase.mesh, 1);
        expected = [2,-1,-1;-1,1,0;-1,0,1]/2;
        testCase.verifyEqual(full(assemble(testCase.blf, testCase.fes)), expected, 'RelTol', eps);
    end
    
    function validateStiffnessMatrix(testCase)
        testCase.blf.a = Constant(testCase.mesh, [1;0;1;2]);
        expected = [4,-1,-3;-2,1,1;-2,0,2]/2;
        testCase.verifyEqual(full(assemble(testCase.blf, testCase.fes)), expected, 'RelTol', eps);
    end
    
    function validateConvection(testCase)
        testCase.blf.a = [];
        testCase.blf.b = Constant(testCase.mesh, [1;1]);
        expected = [-2,1,1;-2,1,1;-2,1,1]/6;
        testCase.verifyEqual(full(assemble(testCase.blf, testCase.fes)), expected, 'RelTol', eps);
    end
    
    function validateReactionWithHigherOrderQuadrature(testCase)
        testCase.blf.c = MeshFunction(testCase.mesh, @(x) 120*sum(x, Dim.Vector));
        expected = [4,3,3;3,8,4;3,4,8];
        testCase.verifyNotEqual(full(assemble(testCase.blf, testCase.fes)), expected);
        testCase.blf.qrc = QuadratureRule.ofOrder(2);
        testCase.verifyNotEqual(full(assemble(testCase.blf, testCase.fes)), expected);
        testCase.blf.qrc = QuadratureRule.ofOrder(3);
        testCase.verifyEqual(full(assemble(testCase.blf, testCase.fes)), expected, 'RelTol', eps*10);
    end
    
    function validateL2rhsScalar(testCase)
        testCase.lf.f = Constant(testCase.mesh, 1);
        expected = [1;1;1]/6;
        testCase.verifyEqual(assemble(testCase.lf, testCase.fes), expected, 'RelTol', eps);
    end
    
    function validateH1rhsMeshFunctionWithHigherOrderQuadrature(testCase)
        testCase.lf.fvec = MeshFunction(testCase.mesh, ...
            @(x) [120;0].*(prod(x, Dim.Vector).*(1-sum(x, Dim.Vector))));
        expected = [-1;1;0];
        testCase.verifyNotEqual(assemble(testCase.lf, testCase.fes), expected);
        testCase.lf.qrfvec = QuadratureRule.ofOrder(2);
        testCase.verifyNotEqual(assemble(testCase.lf, testCase.fes), expected);
        testCase.lf.qrfvec = QuadratureRule.ofOrder(3);
        testCase.verifyEqual(assemble(testCase.lf, testCase.fes), expected, 'RelTol', eps*10);
    end
    
    function validateNeumann(testCase)
        neumannFes = FeSpace(testCase.mesh, LowestOrderH1Fe(), 'neumann', ':');
        nlf = LinearForm();
        nlf.neumann = Constant(testCase.mesh, pi);
        expected = [2;1+sqrt(2);1+sqrt(2)]*(pi/2);
        testCase.verifyEqual(assemble(nlf, neumannFes), expected, 'RelTol', eps*10);
    end
    
    function validateRobin(testCase)
        robinFes = FeSpace(testCase.mesh, LowestOrderH1Fe(), 'dirichlet', 2, 'robin', 1);
        rlf = LinearForm();
        rblf = BilinearForm();
        rblf.robin = Constant(testCase.mesh, pi);
        rblf.qrRobin = QuadratureRule.ofOrder(2, '1D');
        expected = [4,1,1;1,2,0;1,0,2]*(pi/6);
        testCase.verifyEqual(full(assemble(rblf, robinFes)), expected, 'RelTol', eps*10);
    end
end
    
end
