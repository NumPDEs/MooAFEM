classdef TestIterativeSolver < matlab.unittest.TestCase
    
properties (TestParameter)
    p = struct('low', 1, 'medium', 2, 'high', 5);
    variant = struct('CG', ["cg", ""],...
        'direct', ["direct", ""], ...
        'jacobiPCG', ["pcg", "jacobi"], ...
        'iChol', ["pcg", "iChol"], ...
        'additiveSchwarzPCG', ["pcg", "additiveSchwarzLowOrder"], ...
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
            if k < 3, mesh.refineUniform(); end
        end
        
        solver.setupRhs(F);
        normBefore = sqrt((xstar - solver.x)' * A * (xstar - solver.x));
        solver.step();
        normAfterwards = sqrt((xstar - solver.x)' * A * (xstar - solver.x));
        
        testCase.verifyGreaterThan(normBefore, normAfterwards);
    end
    
    function solverIsLinear(testCase, p, variant)
        [mesh, fes, blf, lf] = setupProblem(testCase, p);
        s = variant(1); v = variant(2);
        solver = chooseIterativeSolver(fes, blf, s, v);
        
        for k = 1:3
            [xstar, A, F] = assembleData(testCase, blf, lf, fes);
            solver.setupSystemMatrix(A);
            if k < 3, mesh.refineLocally(k); end
        end
        
        solver.setupRhs([F, -pi*F], 0*[F, F]);
        solver.solve();
%         for i = 1:60
%             solver.step();
%         end
        %xmat = cgs(A, F, 1e-16, 100);
        
         testCase.verifyEqual(-pi*solver.x(:,1), solver.x(:,2), 'RelTol', 2e-5);
         testCase.verifyEqual(solver.x(:,1), xstar, 'RelTol', 2e-5);
         %testCase.verifyEqual(solver.x(:, 1), xmat(:, 1), 'RelTol', 2e-5)
        %testCase.verifyEqual(xmat(:,1), xstar, 'RelTol', 2e-5);
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
%         blf.c = Constant(mesh, 1);
%         blf.qrc = QuadratureRule.ofOrder(2*p);
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
