classdef TestPatchwiseMatrix < matlab.unittest.TestCase
    
properties (TestParameter)
    p = {1, 2, 3, 4, 5};
end

methods (Test)
    function lowestOrderIsDiagonalInverse(testCase)
        mesh = Mesh.loadFromGeometry('Lshape');
        fes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', [], 'neumann', ':');
        blf = BilinearForm();
        blf.a = Constant(mesh, 1);
        blf.c = Constant(mesh, 1);
        for ell = 1:5
            mesh.refineUniform();
            A = assemble(blf, fes);
            B = assemblePatchwise(blf, fes);
            actual = B \ full(diag(A));
            testCase.verifyEqual(actual, ones(getDofs(fes).nDofs, 1), 'RelTol', 10*eps)
        end
    end
    
    function dofsAreInCorrectlyManyPatches(testCase, p)
        mesh = Mesh.loadFromGeometry('Lshape');
        for ell = 1:5, mesh.refineLocally(1); end
        fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', 1, 'neumann', 2);
        blf = BilinearForm();
        blf.a = Constant(mesh, 1);
        blf.qra = QuadratureRule.ofOrder(max(2*(p-1),1));
        A = assemblePatchwise(blf, fes);
        countOccurrence = accumarray(getPatchDofs(A, ':'), 1, [getDofs(fes).nDofs, 1]);
        
        isFree = false(getDofs(fes).nDofs, 1);
        isFree(getFreeDofs(fes)) = true;
        % vertex dofs should occur exactly once (if free)
        vertexDofs = createVertexDofs(fes, ':');
        testCase.verifyEqual(countOccurrence(vertexDofs), 1*isFree(vertexDofs));
        % edge dofs should occur exactly twice (if free)
        edgeDofs = createInnerEdgeDofs(fes, ':');
        testCase.verifyEqual(countOccurrence(edgeDofs), 2*isFree(edgeDofs));
        % inner dofs should occur exactly three times
        elementDofs = createInnerElementDofs(fes, ':');
        testCase.verifyEqual(countOccurrence(elementDofs), 3*isFree(elementDofs));
    end
    
    function patchwiseRelatesToFullMatrix(testCase, p)
        mesh = Mesh.loadFromGeometry('unitsquare');
        for ell = 1:5, mesh.refineLocally(1); end
        fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', [], 'neumann', ':');
        blf = BilinearForm();
        blf.a = Constant(mesh, 1);
        blf.qra = QuadratureRule.ofOrder(max(2*(p-1),1));
        blf.c = MeshFunction(mesh, @(x) 1+x(1,:,:));
        blf.qrc = QuadratureRule.ofOrder(2*p+1);
        A = assemble(blf, fes);
        B = assemblePatchwise(blf, fes);
        for k = 1:10:51
            vertexDof = createVertexDofs(fes, k);
            idx = getPatchDofs(B, k);
            expected = zeros(getDofs(fes).nDofs, 1);
            expected(vertexDof) = 1;
            actual = (B \ expected);
            actual(idx) = A(idx,idx) * actual(idx);
            testCase.verifyEqual(actual, expected, 'AbsTol', sqrt(eps))
        end
    end
end
    
end
