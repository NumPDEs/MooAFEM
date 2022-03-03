classdef TestFeSpace < matlab.unittest.TestCase
    
properties
    dummyMesh
end

properties (TestParameter)
    dirichlet = {1, 2, [1,2], [1,2,1,2,1,1], 'all'};
    expectedDirichlet = {1, 2, [1,2], [1,2], [1,2]};
    fe = struct('H1p1', LowestOrderH1Fe(), 'L2p1', LowestOrderL2Fe(), ...
        'H1p3', HigherOrderH1Fe(4), 'L2p7', HigherOrderL2Fe(7));
    nDof = {3, 1, 15, 36};
    nFreeDof = {0, 1, 3, 36};
end

methods (TestClassSetup)
    function createDummyMesh(testCase)
        testCase.dummyMesh = Mesh([],[], {[]; []});
    end
end

methods (Test, ParameterCombination='sequential')
    function implicitDirichletTakesAllBoundaries(testCase)
        dummyFe = LowestOrderH1Fe();
        fes = FeSpace(testCase.dummyMesh, dummyFe);
        testCase.verifyEqual(fes.dirichlet, [1,2])
    end
    
    function allOptionsForDirichletGiveCorrectBoundaries(testCase, dirichlet, expectedDirichlet)
        dummyFe = LowestOrderH1Fe();
        fes = FeSpace(testCase.dummyMesh, dummyFe, 'dirichlet', dirichlet);
        testCase.verifyEqual(fes.dirichlet, expectedDirichlet);
    end
    
    function wrongDirichletDataThrowsError(testCase)
        dummyFe = LowestOrderH1Fe();
        testCase.verifyError( ...
            @() FeSpace(testCase.dummyMesh, dummyFe, 'dirichlet', [1,3,1,2,1,1]), ...
            'FeSpace:invalidDirichletBoundary')
        testCase.verifyError( ...
            @() FeSpace(testCase.dummyMesh, dummyFe, 'dirichlet', 'whole'), ...
            'FeSpace:invalidDirichletBoundary')
    end
    
    function dofDataIsComputedCorrectly(testCase, fe, nDof, nFreeDof)
        mesh = Mesh.loadFromGeometry('unittriangle');
        fes = FeSpace(mesh, fe);
        testCase.verifyEqual(getDofs(fes).nDofs, nDof)
        testCase.verifyEqual(nnz(getFreeDofs(fes)), nFreeDof)
    end
    
    function dataGetsRecomputedAtMeshRefinement(testCase)
        mesh = Mesh.loadFromGeometry('unitsquare');
        fes = FeSpace(mesh, LowestOrderH1Fe());
        oldDof = getDofs(fes).nDofs;
        oldFree = nnz(getFreeDofs(fes));
        mesh.refineUniform(1);
        newDof = getDofs(fes).nDofs;
        newFree = nnz(getFreeDofs(fes));
        testCase.verifyGreaterThan(newDof, oldDof);
        testCase.verifyGreaterThan(newFree, oldFree);
    end
end
    
end