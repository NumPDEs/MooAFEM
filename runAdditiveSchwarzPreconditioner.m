function leveldata = runAdditiveSchwarzPreconditioner(pMax)

    % Proceed input
    if nargin < 1
        pMax = 1e1;
    end

    % Number of refinement levels
    nLevels = 2;
 
    % Initialize LevelData
    leveldata = LevelData('results');
    leveldata.problem = 'Poisson';
    leveldata.domain = 'UnitSquare';
    leveldata.method = 'Sp_PCG_AdditiveSchwarz';

    % Run experiment for various p
    for p = 1:pMax
        % Load mesh
        mesh = Mesh.loadFromGeometry('unitsquare');
        %mesh.refineUniform();
    
        % Create FE space
        fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', ':');
        % u = FeFunction(fes);
    
        % Setup variational problem
        blf = BilinearForm();
        lf = LinearForm();
        blf.a = Constant(mesh, 1);
        blf.qra = QuadratureRule.ofOrder(max(2*(p-1), 1));
        lf.f = Constant(mesh, 1);
    
        % Choose iterative solver
        solver = chooseIterativeSolver(fes, blf, "pcg", "additiveSchwarzLowOrder");
        % solver = chooseIterativeSolver(fes, blf, "pcg", "additiveSchwarzHighOrder");
        % solver = chooseIterativeSolver(fes, blf, "pcg", "iChol");
        % solver = chooseIterativeSolver(fes, blf, "multigrid", "lowOrderVcycle");
    
        for k = 1:nLevels-1
            % Assemble matrices on intermediate meshes
            A = assemble(blf, fes);
            F = assemble(lf, fes);
            freeDofs = getFreeDofs(fes);
            solver.setupSystemMatrix(A(freeDofs,freeDofs));
            solver.setupRhs(F(freeDofs));  % initializes x0 as well
            % Mesh refinement
            mesh.refineUniform();
        end
    
        % Assemble matrices on fine mesh
        A = assemble(blf, fes);
        if ~issymmetric(A)
            %warning(['symmetrizing, norm(asym(A)) = ', num2str(norm(A - A', 1))])
            A = (A + A') / 2;
        end
        F = assemble(lf, fes);
        freeDofs = getFreeDofs(fes);
        solver.setupSystemMatrix(A(freeDofs,freeDofs));
        solver.setupRhs(F(freeDofs));  % initializes x0 as well

        A = A(freeDofs,freeDofs);

        figure(2);
        spy(A);
        title('Matrix A');

        % Compute C
        C = zeros(length(freeDofs));
        for j = 1:length(freeDofs)
            ej = zeros(length(freeDofs), 1);
            ej(j) = 1;
            C(:,j) = solver.preconditionAction(ej);  
        end
        if ~issymmetric(C)
            % warning(['symmetrizing, norm(asym(C)) = ', num2str(norm(C - C', 1))])
            C = (C + C') / 2;
        end

        % figure(3);
        % spy(C);
        % title('Matrix C');
        % 
        % asymC = (C - C');

        % figure(4);
        % largeAsymC = asymC;
        % largeAsymC(abs(largeAsymC) < max(abs(largeAsymC(:))) / 10) = 0;
        % spy(largeAsymC);
        % title('Large entries of asym(C)');

        % figure(5);
        % asymC(abs(asymC) < 1e-10) = 0;
        % spy(asymC);
        % title('Nonzero entries of asym(C)');

        % Compute conditional number
        CA = C*A;
        lambda = eig(full(CA));
        condition = max(abs(lambda)) / min(abs(lambda));
        % lambdaMin = eigs(@(x)applyCA(x, A, solver), length(freeDofs), 1, 'smallestreal');
        % lambdaMax = eigs(@(x)applyCA(x, A, solver), length(freeDofs), 1, 'largestreal');
        % lambdaMin = eigs(CA, 1, 'smallestabs');
        % lambdaMax = eigs(CA, 1, 'largestabs');
        % condition = lambdaMax / lambdaMin;

        % Print result to commandline
        leveldata.append('p', uint32(p), 'condition', condition);
        leveldata.printLevel();
    end

    % Plot condition
    figure();
    leveldata.plot('p');
end


function CAx = applyCA(x, A, solver)
    % CAx = solver.preconditionAction(A * x);
    % CAx = (A * x);
    CAx = solver.preconditionAction(x);
end
