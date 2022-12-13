classdef TestIterativeSolver < matlab.unittest.TestCase
    
properties (TestParameter)
    p = struct('low', 1, 'medium', 2, 'high', 5);
    variant = struct('CG', ["cg", ""], ...
        'jacobiPCG', ["pcg", "jacobi"], ...
        'additiveSchwarzPCG', ["pcg", "additiveSchwarz"], ...
        'lowMG', ["multigrid", "lowOrderVcycle"], ...
        'highMG', ["multigrid", "highOrderVcycle"]);
end

methods (Test)
    function oneStepDecreasesNorm(testCase, p, variant)
        [mesh, fes, blf, lf] = setupProblem(testCase, p);
        s = variant(1); v = variant(2);
        solver = chooseIterativeSolver(fes, blf, s, v);
        
        for k = 1:3
            [xstar, A, F] = assembleData(testCase, blf, lf, fes);
            solver.setupSystemMatrix(A);
            solver.setupRhs(F);
            if k < 3, mesh.refineLocally(1); end
        end
        
        normBefore = sqrt((xstar - solver.x)' * A * (xstar - solver.x));
        solver.step();
        normAfterwards = sqrt((xstar - solver.x)' * A * (xstar - solver.x));
        
        testCase.verifyGreaterThan(normBefore, normAfterwards);
    end
end

methods (Access=private)
    function [mesh, fes, blf, lf] = setupProblem(~, p)
        mesh = Mesh.loadFromGeometry('unitsquare');
        fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', ':');
        blf = BilinearForm();
        lf = LinearForm();
        blf.a = Constant(mesh, 1);
        blf.qra = QuadratureRule.ofOrder(max(2*(p-1), 1));
        blf.c = Constant(mesh, 1);
        blf.qrc = QuadratureRule.ofOrder(2*p);
        lf.f = Constant(mesh, 1);
    end
    
    function [xstar, A, F] = assembleData(~, blf, lf, fes)
        A = assemble(blf, fes);
        F = assemble(lf, fes);
        freeDofs = getFreeDofs(fes);
        A = A(freeDofs,freeDofs);
        F = F(freeDofs);
        xstar = A \ F;
    end
end
    
    
end
