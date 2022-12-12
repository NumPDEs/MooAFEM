% ******************************************************************************
% Adaptive FEM algorithm with known solution and convergence rates for higher
% order finite elements.
% ******************************************************************************

%% parameters
nDofsMax = 1e4;
theta = 0.5;
pmax = 5;
mu = 0.1;

%% initialization for every polynomial degree
refError2 = cell(pmax,1);
[nElem, nDofs, errEst, time, directSolveTime] = deal(zeros(pmax, 1000));
for p = 1:pmax
    %% setup geometry & spaces
    printLogMessage('*** p = %d (of %d) ***', p, pmax)
    mesh = Mesh.loadFromGeometry('Lshape');
    fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', ':');
    u = FeFunction(fes);
    u.setData(0);
    
    %% set problem data for -\Delta u = 1 on L-shape
    blf = BilinearForm();
    lf = LinearForm();
    
    blf.a = Constant(mesh, 1);
    blf.qra = QuadratureRule.ofOrder(max(2*p-2, 1));
    
    lf.f = MeshFunction(mesh, @(x) 1);
    lf.qrf = QuadratureRule.ofOrder(2*p);
    
    %% set up solver and operator for nested iteration
    P = FeProlongation(fes);
    solver = chooseIterativeSolver(fes, blf, 'multigrid', 'lowOrderVcycle');
    solver.tol = 1e-8;
    solver.maxIter = 100;
    
    %% adaptive loop
    ell = 0;
    meshSufficientlyFine = false;
    while ~meshSufficientlyFine
        %% setup
        tic;
        ell = ell + 1;
        A = assemble(blf, fes);
        rhs = assemble(lf, fes);
        
        freeDofs = getFreeDofs(fes);
        A = A(freeDofs,freeDofs);
        solver.setupSystemMatrix(A);
        solver.setupRhs(rhs(freeDofs), u.data(freeDofs)');
        
        % exact solution as reference
        tic;
        x_ex = A \ rhs(freeDofs);
        directSolveTime(p, ell) = toc;
        refError2{p} = [refError2{p}, (x_ex - solver.x)'*A*(x_ex - solver.x)];

        %% iterative solve
        solverIsConverged = false;
        while ~solverIsConverged
            solver.step();
            refError2{p} = [refError2{p}, (x_ex - solver.x)'*A*(x_ex - solver.x)];
            algEta2 = solver.algEstimator;
            u.setFreeData(solver.x);
            eta2 = estimate(blf, lf, u);
            
            % adaptive stopping criterion with squares
            solverIsConverged = all(algEta2 < mu^2*sum(eta2));
        end
        
        %% estimate error and store data
        nDofs(p,ell) = getDofs(fes).nDofs;
        errEst(p,ell) = sqrt(sum(eta2));
        nElem(p,ell) = mesh.nElements;
        printLogMessage('number of dofs: %d, estimator: %.2e after %d iterations', ...
            nDofs(p,ell), errEst(p,ell), solver.iterationCount);
        
        %% stoping criterion
        meshSufficientlyFine = (nDofs(p,ell) > nDofsMax);
        
        %% refine mesh
        if ~meshSufficientlyFine
            marked = markDoerflerSorting(eta2, theta);
            mesh.refineLocally(marked, 'NVB');
            u.setData(prolongate(P, u));
        end
        time(p,ell) = toc;
    end
end

%% plot convergence rates
figure()
for p = 1:pmax
    idx = find(nDofs(p,:) > 0);
    loglog(nDofs(p,idx), errEst(p,idx), '-o', 'LineWidth', 2, 'DisplayName', ['p=',num2str(p)]);
    hold on
end
for p = unique([1,pmax])
    x = nDofs(p,idx) / nDofs(p,1);
    loglog(nDofs(p,1)*x, errEst(p,1)*x.^(-p/2), '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', ['\alpha = ', num2str(p), '/2'])
end
legend
xlabel('#T_k')
ylabel('\eta')
title('error estimator over number of dofs')

figure()
for p = 1:pmax
    idx = find(nDofs(p,:) > 0);
    loglog(nDofs(p,idx), directSolveTime(p,idx)./nDofs(p,idx), '-^', 'LineWidth', 2, 'DisplayName', ['[direct solve] p=',num2str(p)]);
    hold on
    remainder = time(p,idx)-directSolveTime(p,idx);
    loglog(nDofs(p,idx), remainder./nDofs(p,idx), '-o', 'LineWidth', 2, 'DisplayName', ['[everything else] p=',num2str(p)]);
end
x = nDofs(p,idx) / nDofs(p,1);
loglog(nDofs(p,1)*x, directSolveTime(p,1)*x./x, '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', '\alpha = 1')
legend
xlabel('#T_k')
ylabel('t [s]')
title('runtime over number of dofs')

figure()
for p = 1:pmax
    idx = 1:numel(refError2{p});
    semilogy(idx, refError2{p}, '-o', 'LineWidth', 2, 'DisplayName', ['p=',num2str(p)]);
    hold on
end
legend
xlabel('# solver step')
ylabel('difference to real solution')

%% local function for residual a posteriori error estimation
% \eta(T)^2 = h_T^2 * || \Delta u ||_{L^2(T)}^2
%               + h_T * || [[(Du - fvec) * n]] ||_{L^2(E) \cap \Omega}^2
function indicators = estimate(blf, lf, u)
p = u.fes.finiteElement.order;
mesh =  u.fes.mesh;
trafo = getAffineTransformation(mesh);

% compute volume residual element-wise
% For p=1, the diffusion term vanishes in the residual.
if p == 1
    f = lf.f;
else
    f = CompositeFunction(@(v,f) (v(1,:,:) + v(4,:,:) + f).^2, Hessian(u), lf.f);
end
qr = QuadratureRule.ofOrder(2*p);
volumeResidual = integrateElement(f, qr);

% compute edge residual edge-wise
qr = QuadratureRule.ofOrder(max(p-1, 1), '1D');
edgeResidual = integrateNormalJump(Gradient(u), qr, ...
    @(j) zeros(size(j)), {}, vertcat(mesh.boundaries{:}), ... % no jump on dirichlet edges
    @(j) j.^2, {}, ':');                                      % normal jump on inner edges

% combine the residuals suitably
hT = sqrt(trafo.area);
indicators = hT.^2 .* volumeResidual + ...
    hT .* sum(edgeResidual(mesh.element2edges), Dim.Vector);
end
