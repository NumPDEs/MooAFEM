classdef TestFiniteElementFunctions < matlab.unittest.TestCase
    
properties
    u
    phi
    Dphi
    vec1
    vec2
end

properties (ClassSetupParameter)
    fe = struct('H1p1', LowestOrderH1Fe(), 'L2p1', LowestOrderL2Fe(), ...
        'H1p3', HigherOrderH1Fe(4), 'L2p7', HigherOrderL2Fe(7));
end

methods (TestClassSetup)
    function createConstantFeFunction(testCase, fe)
        mesh = Mesh.loadFromGeometry('unitsquare');
        fes = FeSpace(mesh, fe);
        testCase.u = FeFunction(fes);
        testCase.u.setData(pi);
        testCase.phi = TestFunction(fes);
        testCase.Dphi = TestFunctionGradient(fes);
        testCase.vec1 = Constant(mesh, [1;2]/3);
        testCase.vec2 = Constant(mesh, [1;-1]/2);
    end
end

methods (Test)
    function integralOfConstantFunctionIsConstant(testCase)
        qr = QuadratureRule.ofOrder(1);
        integral = sum(integrateElement(testCase.u, qr), Dim.Elements);
        testCase.verifyEqual(integral, pi, 'RelTol', eps*100)
    end
    
    function gradientOfConstantFunctionIsZero(testCase)
        qr = QuadratureRule.ofOrder(1);
        f = CompositeFunction(@(Du,vec) abs(vectorProduct(Du, vec)), ...
            Gradient(testCase.u), testCase.vec1);
        integral = sum(integrateElement(f, qr), Dim.Elements);
        testCase.verifyEqual(integral, 0, 'AbsTol', sqrt(eps))
    end
    
    function hessianOfConstantFunctionIsZero(testCase)
        qr = QuadratureRule.ofOrder(1);
        f = CompositeFunction( ...
            @(D2u,vec1,vec2) abs(vectorProduct(vec2, vectorProduct(D2u, vec1, [2,2], [2,1]))), ...
            Hessian(testCase.u), testCase.vec1, testCase.vec2);
        integral = sum(integrateElement(f, qr), Dim.Elements);
        testCase.verifyEqual(integral, 0, 'AbsTol', sqrt(eps))
    end
    
    function testFunctionAndGradientHaveCorrectDimension(testCase)
        p = testCase.phi.fes.finiteElement.order;
        nElements = testCase.phi.fes.mesh.nElements;
        bary = Barycentric2D(ones(3,4)/3);
        correctElement = 4 * (p+1)*(p+2)/2;
        testCase.verifyNumElements(eval(testCase.phi, bary), correctElement);
        testCase.verifyNumElements(eval(testCase.Dphi, bary), 2*nElements*correctElement);
    end
    
    function setFreeDataWorksWithScalars(testCase)
        scalar = 3;
        testCase.u.setFreeData(scalar);
        nFreeDofs = length(getFreeDofs(testCase.u.fes));
        testCase.verifyEqual(nnz(testCase.u.data == scalar), nFreeDofs)
    end
    
    function prolongationPreservesConstants(testCase)
        mesh = Mesh.loadFromGeometry('unitsquare');
        fes = FeSpace(mesh, testCase.u.fes.finiteElement);
        v = FeFunction(fes);
        testCase.refineAndInterpolate(mesh, fes, v, pi);
        qr = QuadratureRule.ofOrder(1);
        integral = sum(integrateElement(v, qr), Dim.Elements);
        testCase.verifyEqual(integral, pi, 'AbsTol', sqrt(eps))
    end
    
    function generalProlongationPreservesBasisFunctions(testCase)
        mesh = Mesh.loadFromGeometry('unitsquare');
        fes = FeSpace(mesh, testCase.u.fes.finiteElement);
        qr = QuadratureRule.ofOrder(max(fes.finiteElement.order,1));
        v = FeFunction(fes);
        P = GeneralFeProlongation(fes);
        for k = 1:min(5, size(getDofs(fes).element2Dofs, 1))
            testCase.setToElementwiseBasisFunction(v, k);
            before = sum(integrateElement(v, qr));
            mesh.refineLocally(1, 'NVB');
            v.setData(prolongate(P, v));
            after = sum(integrateElement(v, qr));
            testCase.verifyEqual(before, after, 'AbsTol', sqrt(eps))
        end
    end
end

methods (Access='private')
    function refineAndInterpolate(~, mesh, fes, v, val)
        v.setData(val);
        if isa(fes.finiteElement, 'LowestOrderH1Fe')
            P = LoH1Prolongation(fes);
        elseif isa(fes.finiteElement, 'LowestOrderL2Fe')
            P = LoL2Prolongation(fes);
        else
            P = GeneralFeProlongation(fes);
        end
        mesh.refineUniform();
        v.setData(prolongate(P, v));
    end
    
    function setToElementwiseBasisFunction(~, v, k)
        v.setData(0);
        data = v.data;
        data(getDofs(v.fes).element2Dofs(k,:)) = 1;
        v.setData(data);
    end
end
    
end
