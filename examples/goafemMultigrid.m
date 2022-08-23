% ******************************************************************************
% Example of a GOAFEM algorithm with iterative solver.
% ******************************************************************************

%% paramters
nDofsMax = 1e4;
theta = 0.5;
p = 4;

%% setup geometry & spaces
[nElem, nDofs, nIterPrimal, nIterDual, goalErrEst] = deal(zeros(1, 1000));

printLogMessage('*** GOAFEM with p = %d and iterative solver ***', p)
mesh = Mesh.loadFromGeometry('Lshape');
mesh.refineUniform(1, 'RGB');
fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', ':');
u = FeFunction(fes);
z = FeFunction(fes);
u.setData(0);
z.setData(0);

%% set problem data
% primal: -\Delta u = "lorentzian peak 1"
% dual:   -\Delta z = "lorentzian peak 2"
blf = BilinearForm(fes);
blf.a = Constant(mesh, 1);
blf.qra = QuadratureRule.ofOrder(max(2*p-2, 1));

lfF = LinearForm(fes);
lfF.f = MeshFunction(mesh, @(x) lorentzian(x, [0.7;0.7], 1e-1));
lfF.qrf = QuadratureRule.ofOrder(2*p);

lfG = LinearForm(fes);
lfG.f = MeshFunction(mesh, @(x) lorentzian(x, [0.2;0.3], 1e-2));
lfG.qrf = QuadratureRule.ofOrder(2*p);

%% set up solver and lifting operator for nested iteration
P = FeProlongation(fes);
solver = pLoc_MG(P);
solver.tol = 1e-4;
solver.maxIter = 1000;

%% adaptive loop
ell = 0;
meshSufficientlyFine = false;
while ~meshSufficientlyFine
    ell = ell + 1;

    %% assemble & solve FEM system iteratively
    freeDofs = getFreeDofs(fes);
    p1freedofs = freeDofs(freeDofs<=mesh.nCoordinates);
    mesh.freeVert{end+1} = p1freedofs;
    A = assemble(blf);
    rhs = [assemble(lfF), assemble(lfG)];
    uz0 = [u.data', z.data'];
    if ell == 1
        uz0(freeDofs,:) = A(freeDofs,freeDofs)\rhs(freeDofs,:); %direct solver on coarse mesh
        u.setFreeData(uz0(freeDofs,1));
        z.setFreeData(uz0(freeDofs,2));
    else
        solver.setupLinearSystem(A, rhs(freeDofs,:), uz0(freeDofs,:));
        while ~all(isConverged(solver))
            solver.step();
        end

%         exsol = A(freeDofs,freeDofs)\rhs(freeDofs,:);
%         u.setFreeData(exsol(:,1));
%         z.setFreeData(exsol(:,2));
         u.setFreeData(solver.x(:,1));
         z.setFreeData(solver.x(:,2));

    end

    %% estimate error and store data
    eta2 = estimate(blf, lfF, u);
    zeta2 = estimate(blf, lfG, z);
    nDofs(ell) = getDofs(fes).nDofs;
    nElem(ell) = mesh.nElements;
    if ell == 1
        nIterPrimal(ell) = 1;
        nIterDual(ell) = 1;
    else
        nIterPrimal(ell) = solver.iterationCount(1);
        nIterDual(ell) = solver.iterationCount(2);
    end
    goalErrEst(ell) = sqrt(sum(eta2)*sum(zeta2));
    printLogMessage('number of dofs: %d, iterations (primal/dual): %d / %d, estimator: %.2e', ...
        nDofs(ell), nIterPrimal(ell), nIterDual(ell), goalErrEst(ell));

    %% stoping criterion
    meshSufficientlyFine = (nDofs(ell) > nDofsMax);

    %% refine mesh and transfer solutions to finer mesh for nested iteration
    if ~meshSufficientlyFine
        marked = markGoafemMS(eta2, zeta2, theta);
        clear meshtemp;
        save('meshtemp','mesh');
        meshtemp = load('meshtemp.mat');
        meshtemp = meshtemp.mesh;
        meshtemp.intergridMatrix(marked, 'NVB');
        interMat = meshtemp.intergrid;
        vertices = meshtemp.locVert;
        mesh.refineLocally(marked, 'NVB');
        mesh.intergrid = interMat;
        mesh.locVert = vertices;
        u.setData(prolongate(P, u));
        z.setData(prolongate(P, z));



    end
end

delete meshtemp.mat

%% plot convergence rates
figure()
idx = 1:ell;
x = cumsum(nDofs(idx).*max(nIterPrimal(idx), nIterDual(idx)));
yyaxis left
loglog(x, goalErrEst(idx), '-x', 'LineWidth', 2, 'MarkerSize', 8);
hold on
loglog(x, goalErrEst(1)*(x/x(1)).^(-p), '--', 'LineWidth', 2, 'Color', 'k')
ylabel('goal error estimator of final iterates')
yyaxis right
stem(x, nIterPrimal(idx), 's', 'filled')
stem(x, nIterDual(idx), 'o', 'filled')
ylabel('number of iterations')
legend('\eta_l \zeta_l', ['\alpha=', num2str(p)], 'primal iterations', 'dual iteration')
xlabel('total work')
title(['GOAFEM p=', num2str(p)])

%% local function for residual a posteriori error estimation
% \eta(T)^2 = h_T^2 * || \Delta u ||_{L^2(T)}^2
%               + h_T * || [[(Du - fvec) * n]] ||_{L^2(E) \cap \Omega}^2
function indicators = estimate(blf, lf, u)
p = blf.fes.finiteElement.order;
mesh =  blf.fes.mesh;
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

% Lorentz-peak 1/(eps+||x-x0||^2)
function y = lorentzian(x, x0, eps)
dx = x - x0;
y = 1 ./ (eps + vectorProduct(dx,dx));
end