classdef TestIterativeSolver < matlab.unittest.TestCase
    
properties (TestParameter)
    p = {1, 2, 3, 4, 5};
end

methods (Test)
    function normDecreases(testCase, p)
        mesh = Mesh.loadFromGeometry('unitsquare');
        fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', ':');
        u = FeFunction(fes);
        blf = BilinearForm();
        lf = LinearForm();
        blf.a = Constant(mesh, 1);
        blf.c = Constant(mesh, 1);
        lf.f = Constant(mesh, 1);
        solver = chooseIterativeSolver(fes, blf, 'multigrid', 'lowOrderVcycle');
        solver.tol = 1e-8;
        solver.maxIter = 100;
        norm = [];
        for ell = 1:2
            mesh.refineLocally(1);
            A = assemble(blf, fes);
            rhs = assemble(lf, fes);
            freeDofs = getFreeDofs(fes);
            A = A(freeDofs,freeDofs);
            u.setData(0);
            solver.setupSystemMatrix(A);
            solver.setupRhs(rhs(freeDofs), u.data(freeDofs)');
            exact = A \ rhs(freeDofs);
            solver.solve();
            norm(ell) = (exact - solver.x)' * A * (exact - solver.x);
        end
        normDifference = norm(1) - norm(2);
        testCase.verifyGreaterThan(normDifference, 0, 'RelTol', 10*eps);
    end
end
    
end
