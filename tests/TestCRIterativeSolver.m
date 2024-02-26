classdef TestCRIterativeSolver < matlab.unittest.TestCase
    
properties (TestParameter)
    variant = struct('CG', ["cg", ""],...
        'direct', ["direct", ""], ...
        'MG', ["multigrid", "lowOrderVcycle"]);
end

methods (Test)
    function firstStepDecreasesNorm(testCase, variant)
        [~, fes, blf, lf] = setupProblem(testCase);
        s = variant(1); v = variant(2);
        solver = IterativeSolver.choose(fes, blf, s, v);
        
        [xstar, A, F] = assembleData(testCase, blf, lf, fes);
        solver.setupSystemMatrix(A);
        solver.setupRhs(F);

        normBefore = sqrt((xstar - solver.x)' * A * (xstar - solver.x));
        solver.step();
        normAfterwards = sqrt((xstar - solver.x)' * A * (xstar - solver.x));
        
        testCase.verifyGreaterThan(normBefore, normAfterwards);
    end

    function laterStepDecreasesNorm(testCase, variant)
        [mesh, fes, blf, lf] = setupProblem(testCase);
        s = variant(1); v = variant(2);
        solver = IterativeSolver.choose(fes, blf, s, v);
        
        for k = 1:3
            [xstar, A, F] = assembleData(testCase, blf, lf, fes);
            solver.setupSystemMatrix(A);
            if k < 3, mesh.refineUniform(); end
        end
        
        solver.setupRhs(F);
        normBefore = sqrt((xstar - solver.x)' * A * (xstar - solver.x));
        solver.step();
        normAfterwards = sqrt((xstar - solver.x)' * A * (xstar - solver.x));
        
        testCase.verifyGreaterThan(normBefore, normAfterwards);
    end

    function solverIsLinear(testCase, variant)
        [mesh, fes, blf, lf] = setupProblem(testCase);
        s = variant(1); v = variant(2);
        solver = IterativeSolver.choose(fes, blf, s, v);
        
        for k = 1:3
            [xstar, A, F] = assembleData(testCase, blf, lf, fes);
            solver.setupSystemMatrix(A);
            if k < 3, mesh.refineLocally(k); end
        end
        
        solver.setupRhs([F, -pi*F], 0*[F, F]);
        solver.solve();

        testCase.verifyEqual(-pi*solver.x(:,1), solver.x(:,2), 'RelTol', 2e-5);
        testCase.verifyEqual(solver.x(:,1), xstar, 'RelTol', 2e-5);
    end
end

methods (Access=private)
    function [mesh, fes, blf, lf] = setupProblem(~)
        mesh = Mesh.loadFromGeometry('unitsquare');
        fes = FeSpace(mesh, LowestOrderCRFe(), 'dirichlet', ':');
        blf = BilinearForm();
        lf = LinearForm();
        blf.a = Constant(mesh, 1);
        blf.c = Constant(mesh, 1);
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
