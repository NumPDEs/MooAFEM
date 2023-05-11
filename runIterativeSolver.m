function leveldata = runIterativeSolver(maxNiter)

    % Proceed input
    if nargin < 1
        maxNiter = 1e2;
    end

    % Polynomial degree
    p = 2;

    % Load mesh
    mesh = Mesh.loadFromGeometry('unitsquare');

    % Create FE space
    fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', ':');
    u = FeFunction(fes);

    % Setup variational problem
    blf = BilinearForm();
    lf = LinearForm();
    blf.a = Constant(mesh, 1);
    blf.qra = QuadratureRule.ofOrder(max(2*(p-1), 1));
    lf.f = Constant(mesh, 1);

    % Choose iterative solver
    solver = chooseIterativeSolver(fes, blf, "pcg", "additiveSchwarz");

    % Initialize LevelData
    leveldata = LevelData('results');
    leveldata.problem = 'Poisson';
    leveldata.domain = 'UnitSquare';
    leveldata.method = 'Sp_PCG_AdditiveSchwarz';

    for k = 1:4
        mesh.refineUniform();
        % Preparation
        A = assemble(blf, fes);
        F = assemble(lf, fes);
        freeDofs = getFreeDofs(fes);
        solver.setupSystemMatrix(A(freeDofs,freeDofs));
        solver.setupRhs(F(freeDofs));  % initializes x0 as well
    end

    % Refinement loop
    nIter = 1;
    while true
        ndof = getDofs(fes).nDofs;

        % Compute exact solution
        xstar = A(freeDofs,freeDofs) \ F(freeDofs);

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
        leveldata.setAbsolute(nIter, 'nIter', nIter, 'contraction', contraction);

        % Print level information
        leveldata.printLevel();

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