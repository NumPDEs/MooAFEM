function leveldata = runIterativeSolver(maxNiter)

    % Proceed input
    if nargin < 1
        maxNiter = 1e2;
    end

    % Polynomial degree
    p = 4;

    % Number of refinement levels
    nLevels = 6;

    % Load mesh
    mesh = Mesh.loadFromGeometry('unitsquare');
    mesh.refineUniform();

    % Create FE space
    fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', ':');
    u = FeFunction(fes);
    ustar = FeFunction(fes);

    % Setup variational problem
    blf = BilinearForm();
    lf = LinearForm();
    blf.a = Constant(mesh, 1);
    blf.qra = QuadratureRule.ofOrder(max(2*(p-1), 1));
    lf.f = Constant(mesh, 1);

    % Choose iterative solver
    solver = IterativeSolver.choose(fes, blf, "pcg", "additiveSchwarzLowOrder");
    % solver = IterativeSolver.choose(fes, blf, "pcg", "additiveSchwarzHighOrder");
    % solver = IterativeSolver.choose(fes, blf, "pcg", "iChol");
    % solver = IterativeSolver.choose(fes, blf, "multigrid", "lowOrderVcycle");

    % Initialize LevelData
    leveldata = LevelData('results');
    leveldata.problem = 'Poisson';
    leveldata.domain = 'UnitSquare';
    leveldata.method = 'Sp_PCG_AdditiveSchwarz';

    for k = 1:nLevels-1
        % Output
        fprintf('Assemble data for level = %d\n', k);
        % Assemble matrices on intermediate meshes
        A = assemble(blf, fes);
        F = assemble(lf, fes);
        freeDofs = getFreeDofs(fes);
        solver.setupSystemMatrix(A(freeDofs,freeDofs));
        solver.setupRhs(F(freeDofs));  % initializes x0 as well
        % Mesh refinement
        mesh.refineUniform();
    end

    fprintf('Assemble data for level = %d\n', nLevels);

    % Assemble matrices on fine mesh
    A = assemble(blf, fes);
    F = assemble(lf, fes);
    freeDofs = getFreeDofs(fes);
    solver.setupSystemMatrix(A(freeDofs,freeDofs));
    solver.setupRhs(F(freeDofs));  % initializes x0 as well

    % Compute exact solution
    fprintf('Do exact solve\n');
    xstar = A(freeDofs,freeDofs) \ F(freeDofs);

    % Refinement loop
    nIter = 1;
    while true
        ndof = getDofs(fes).nDofs;

        % Compute energy norm
        errorU = sqrt((xstar - solver.x)' * A(freeDofs,freeDofs) * (xstar - solver.x));

        % Compute contraction factor for energy norm
        if nIter == 1
            contraction = Inf;
        else
            contraction = errorU / errorUPrevious;
        end

        % Store variables
        leveldata.append('ndof', uint32(ndof), ...
                         'errorU', errorU);
        leveldata.setAbsolute(nIter, 'nIter', uint32(nIter), ...
                                     'contraction', contraction);

        % Print level information
        leveldata.printLevel();

        % Plot solution
        % if p == 2
        %     x = zeros(mesh.nCoordinates + mesh.nEdges, 1);
        %     x(freeDofs) = solver.x;
        %     figure(1);
        %     plotS2(mesh.clone(), x);
        % 
        %     figure(2);
        %     x(freeDofs) = xstar;
        %     plotS2(mesh.clone(), x);
        %     title("Exact algebraic solution")
        % else
        %     figure(1);
        %     u.setFreeData(solver.x);
        %     plot(u);
        % 
        %     figure(2);
        %     ustar.setFreeData(xstar);
        %     plot(ustar);
        % end
        % pause(2)

        % Break condition
        % if ndof > maxNdof
        %     break;
        % end

        % Mesh refinement
        % mesh.refineUniform();

        if nIter > maxNiter
            break;
        end

        % One iterative solver step
        solver.step();

        % Update iteration counter
        nIter = nIter + 1;
        errorUPrevious = errorU;
    end

    % Plot solution
    u.setFreeData(solver.x);
    plot(u);

    % Plot error
    figure();
    leveldata.plot('nIter');

    % Plot contraction
    figure();
    leveldata.plotAbsolute('nIter');
end
