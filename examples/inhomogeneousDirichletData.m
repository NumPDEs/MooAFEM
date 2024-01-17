% ******************************************************************************
% Kellogg interface benchmark problem with inhomogeneous Dirichlet boundary
% conditions
% ******************************************************************************

%% parameters
nDofsMax = 1e5;
p = 1;
theta = 0.5;
lambda = 1e-1;

%% setup geometry & spaces
printLogMessage('*** p = %d ***', p)
mesh = Mesh.loadFromGeometry('bigsquare');
mesh.refineUniform(4);
fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', ':');
u = FeFunction(fes);
u.setData(0);
uExact = FeFunction(fes);
uLifting = FeFunction(fes);
uLifting.setFreeData(0);

duD_dt = MeshFunction(mesh, @(x) duD_dt_fct(x));

%% set problem data
blf = BilinearForm();
lf = LinearForm();

ncFes = FeSpace(mesh, LowestOrderL2Fe());
a = FeFunction(ncFes);
blf.a = a;
blf.qra = QuadratureRule.ofOrder(max(2*p-2, 1));

lf.f = Constant(mesh, 0);
lf.qrf = QuadratureRule.ofOrder(p);

%% set up solver and operator for nested iteration
[solver, P] = IterativeSolver.choose(fes, blf, "multigrid", "lowOrderVcycle");
solver.tol = 1e-8;
solver.maxIter = 100;

%% adaptive loop
ell = 0;
cost = 0;
time = 0;
meshSufficientlyFine = false;
while ~meshSufficientlyFine
    tic;
    %% setup
    ell = ell + 1;

    % evaluate piecewise constant coefficent function
    a.setData(nodalInterpolation(MeshFunction(mesh, @(x) diffusion(x)), ncFes));

    % assemble linear system
    A = assemble(blf, fes);
    F = assemble(lf, fes);
    freeDofs = getFreeDofs(fes);
    solver.setupSystemMatrix(A(freeDofs, freeDofs));
    nDofs(ell) = getDofs(fes).nDofs; %#ok<SAGROW>

    % include Dirichlet boundary conditions
    [rhs, uDirData] = setDirichletData(@uExact_fct, fes, A, F);
    u.setFixedData(uDirData);

    % exact algebraic solution
    uEllExact = A(freeDofs, freeDofs) \ rhs(freeDofs);

    % iterative algebraic solver step
    solver.setupRhs(rhs(freeDofs), u.data(freeDofs)');

    solverIsConverged = false;
    while ~solverIsConverged
        cost = cost + nDofs(ell);
        solver.step();

        jIncrementNorm = energyNorm(A(freeDofs, freeDofs), u.data(freeDofs)' - solver.x);

        u.setFreeData(solver.x);
        eta2 = estimate(blf, lf, u, duD_dt);
        estimator = sqrt(sum(eta2));

        % algebraic stopping criterion from AISFEM
        solverIsConverged = all(jIncrementNorm < lambda * estimator);
    end

    %% estimate error and store data
    errEst(ell) = estimator; %#ok<SAGROW>
    complexity(ell) = cost; %#ok<SAGROW>

    printLogMessage('number of dofs: %d, estimator: %.2e after %d iterations', ...
        nDofs(ell), errEst(ell), solver.iterationCount);

    %% stoping criterion
    meshSufficientlyFine = (nDofs(ell) > nDofsMax);

    %% refine mesh
    if ~meshSufficientlyFine
        marked = markDoerflerBinning(eta2, theta);
        mesh.refineLocally(marked, 'NVB');
        u.setData(prolongate(P, u));
    end
    time = time + toc;
    times(ell) = time; %#ok<SAGROW>
end
idx = find(nDofs > 0);
estimator = errEst(idx);
ndofs = nDofs(idx);

%% plot convergence rates
figure()
for k = p
    idx = find(complexity > 0);
    loglog(complexity(idx), errEst(idx), '-o', 'LineWidth', 2, 'DisplayName', ['k=',num2str(k)]);
    hold on
end
k=1;
x = complexity(k,idx) / complexity(k,1);
loglog(complexity(k,1)*x, errEst(k,1)*x.^(-p/2), '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', ['\alpha = ', num2str(p), '/2'])

legend
xlabel('\(\sum_{|\ell^{\prime}, k^{\prime}, j^{\prime}| \le |\ell, \underline{k}, \underline{j}|} \mathrm{dim} \mathcal{X}_{\ell^\prime}\)', 'interpreter', 'latex')
ylabel('\(\eta_\ell(u_{\ell}^{\underline{k}, \underline{j}})\)', 'interpreter', 'latex')
title('error estimator over computational costs')

figure;
plot(mesh);

figure;
plot(u);

function indicators = estimate(blf, lf, u, duD_dt)
p = u.fes.finiteElement.order;
mesh =  u.fes.mesh;
trafo = getAffineTransformation(mesh);

% compute volume residual element-wise
f = CompositeFunction(@(a, v, f) (a .* (v(1,:,:) + v(4,:,:)) + f).^2, blf.a, Hessian(u), lf.f);
qr = QuadratureRule.ofOrder(2*p);
volumeResidual = integrateElement(f, qr);

% compute edge residual edge-wise
f = CompositeFunction(@(a, Du) a .* Du, blf.a, Gradient(u));
qr = QuadratureRule.ofOrder(max(p-1, 1), '1D');
edgeResidual = integrateNormalJump(f, qr, ...
    @(j) zeros(size(j)), {}, vertcat(mesh.boundaries{:}), ... % no jump on dirichlet edges
    @(j) j.^2, {}, ':');                                      % normal jump on inner edges

% compute boundary-data oscillations
ncfes = FeSpace(mesh, LowestOrderL2Fe());
duD_dt_mean = FeFunction(ncfes);
duD_dt_mean.setData(0);
idx = getCombinedBndEdges(u.fes.mesh, u.fes.bnd.dirichlet);
dirichletEdgeData = integrateEdge(duD_dt, qr, idx);
% WARNING
% The following only works if triangles have at most one boundary edge!
edgeData = duD_dt_mean.data;
edgeData(mesh.boundaryElems) = dirichletEdgeData;
duD_dt_mean.setData(edgeData);
f = CompositeFunction(@(dudt, mean) (dudt - mean), duD_dt, duD_dt_mean);
boundaryDataError = integrateJump(f, qr, @(j) j.^2, {}, idx);

% maybe prefer something like:
% boundaryDataError = zeros(size(edgeResidual));
% boundaryDataError(idx) = integrateEdge(f, qr, idx);

% combine the terms suitably
hT = sqrt(trafo.area);
indicators = hT.^2 .* volumeResidual + ...
    hT .* sum(edgeResidual(mesh.element2edges), Dim.Vector) + ...
    hT .* sum(boundaryDataError(mesh.element2edges), Dim.Vector);
end

function e = energyNorm(A, v)
e = sqrt(v' * A * v);
end

%% PROBLEM INPUT DATA
function val = diffusion(x)
    % Constants
    R = 161.4476387975881;
    % Evaluation
    Ind1 = 0 < prod(x, Dim.Vector);
    Ind2 = prod(x, Dim.Vector) <= 0;
    val = zeros(size(Ind1));
    val(Ind1) = R;
    val(Ind2) = 1;
end

%% EXACT SOLUTION
% Auxiliary functions
function val = mu(theta)
    % Constants
    gamma = 0.1;
    rho = pi / 4;
    sigma = -14.92256510455152;

    val = zeros(size(theta));
    Ind1 = (0 <= theta & theta <= pi/2);
    Ind2 = (pi/2 <= theta & theta <= pi) & ~Ind1;
    Ind3 = (pi <= theta & theta <= 3/2*pi) & ~Ind1 & ~Ind2;
    Ind4 = (3/2*pi <= theta & theta <= 2*pi) & ~Ind1 & ~Ind2 & ~Ind3;
    val(Ind1) = cos((pi/2 - sigma) * gamma) .* cos((theta(Ind1) - pi/2 + rho) * gamma);
    val(Ind2) = cos(rho * gamma) * cos((theta(Ind2) - pi + sigma) * gamma);
    val(Ind3) = cos(sigma * gamma) * cos((theta(Ind3) - pi - rho) * gamma);
    val(Ind4) = cos((pi/2 - rho) * gamma) * cos((theta(Ind4) - 3/2*pi - sigma) * gamma);
end
% Solution
function val = uExact_fct(x)
    % Constants
    gamma = 0.1;
    % Transform to polar coordinates
    [phi, r] = cart2pol(x(1,:,:), x(2,:,:));
    Ind = phi < 0;
    phi(Ind) = phi(Ind) + 2 * pi;
    % Evaluation
    val = r.^gamma .* mu(phi);
end
% Derivatives
function val = dmu_dtheta(theta)
    % Constants
    gamma = 0.1;
    rho = pi / 4;
    sigma = -14.92256510455152;
    % Evaluation
    val = zeros(size(theta));
    Ind1 = (0 <= theta & theta <= pi/2);
    Ind2 = (pi/2 <= theta & theta <= pi) & ~Ind1;
    Ind3 = (pi <= theta & theta <= 3/2*pi) & ~Ind1 & ~Ind2;
    Ind4 = (3/2*pi <= theta & theta <= 2*pi) & ~Ind1 & ~Ind2 & ~Ind3;
    val(Ind1) = - gamma * cos((pi/2 - sigma) * gamma) .* sin((theta(Ind1) - pi/2 + rho) * gamma);
    val(Ind2) = - gamma * cos(rho * gamma) * sin((theta(Ind2) - pi + sigma) * gamma);
    val(Ind3) = - gamma * cos(sigma * gamma) * sin((theta(Ind3) - pi - rho) * gamma);
    val(Ind4) = - gamma * cos((pi/2 - rho) * gamma) * sin((theta(Ind4) - 3/2*pi - sigma) * gamma);
end
function val = gradUExact_fct(x) %#ok<DEFNU>
    % Constants
    gamma = 0.1;
    % Transform to polar coordinates
    [phi, r] = cart2pol(x(1,:,:), x(2,:,:));
    Ind = (phi < 0);
    phi(Ind) = phi(Ind) + 2 * pi;
    % Evaluation
    du_dr = gamma * r.^(gamma-1) .* mu(phi);
    du_dphi = r.^gamma .* dmu_dtheta(phi);
    du_dx = du_dr .* cos(phi) - du_dphi .* sin(phi) ./ r;
    du_dy = du_dr .* sin(phi) + du_dphi .* cos(phi) ./ r;
    du_dx(r == 0) = Inf;
    du_dy(r == 0) = Inf;
    val = [du_dx, du_dy];
end

function val = duD_dt_fct(x)
    % Determine boundary parts
    right = (x(1, :, :) == 1);
    left = (x(1, :, :) == -1);
    top = (x(2, :, :) == 1);
    bottom = (x(2, :, :) == -1);
    % Constants
    gamma = 0.1;
    % Transform to polar coordinates
    [phi, r] = cart2pol(x(1,:,:), x(2,:,:));
    Ind = (phi < 0);
    phi(Ind) = phi(Ind) + 2 * pi;
    % Evaluation
    du_dr = gamma * r.^(gamma-1) .* mu(phi);
    du_dphi = r.^gamma .* dmu_dtheta(phi);
    du_dx = du_dr .* cos(phi) - du_dphi .* sin(phi) ./ r;
    du_dy = du_dr .* sin(phi) + du_dphi .* cos(phi) ./ r;
    du_dx(r == 0) = Inf;
    du_dy(r == 0) = Inf;
    % Assemble
    val = zeros(size(x(1, :, :)));
    val(right) = du_dy(right);
    val(left) = - du_dy(left);
    val(top) = - du_dx(top);
    val(bottom) = - du_dx(bottom);
end


% function y = dudn(x)
%     % determine boundary parts
%     right  = x(1, :, :) == 1;
%     left   = x(1, :, :) == -1;
%     top    = x(2, :, :) == 1;
%     bottom = x(2, :, :) == -1;
% 
% 
%     [th, r] = cart2pol(x(1,:,:), x(2,:,:));
%     th = th + 2*pi*(th < 0);
% 
%     gamma = 0.0009;
%     %jump = 2001405.429972;
%     rho = 0.785398;
%     sigma = -1744.543854;
%     mu = zeros(size(th));
%     dmu = zeros(size(th));
%     for i = 1:length(th)
%         if th(i) <= pi/2
%             mu(i) = cos(gamma*pi/2 - gamma*sigma)*cos((th(i)-pi/2+rho)*gamma);
%             dmu(i) = -gamma*cos(gamma*pi/2 - gamma*sigma)*sin((th(i)-pi/2+rho)*gamma);
%         elseif (th(i) > pi/2) && (th(i) <= pi)
%             mu(i) = cos(rho*gamma)*cos((th(i)-pi+sigma)*gamma);
%             dmu(i) = -gamma*cos(rho*gamma)*sin((th(i)-pi+sigma)*gamma);
%         elseif (th(i) > pi) && (th(i) <= 3/2*pi)
%             mu(i) = cos(sigma*gamma)*cos((th(i)-pi+rho)*gamma);
%             dmu(i) = -gamma*cos(sigma*gamma)*sin((th(i)-pi+rho)*gamma);
%         else
%             mu(i) = cos(gamma*pi/2 - gamma*rho)*cos((th(i)-3*pi/2-sigma)*gamma);
%             dmu(i) = -gamma*cos(gamma*pi/2 - gamma*rho)*sin((th(i)-3*pi/2-sigma)*gamma);
%         end
%     end
%     dudr = (r.^(gamma-1)).* gamma .* mu;
%     dudth = (r.^(gamma)).* dmu ;
%     dudx = (cos(th).*dudr - (sin(th)./r) .* dudth);
%     dudy = (sin(th).*dudr + (cos(th)./r) .* dudth);
% 
%     y = zeros(size(x(1, :, :)));
%     y(right)  = dudx(right);
%     y(left)   =-dudx(left);
%     y(top)    = dudy(top);
%     y(bottom) =-dudy(bottom);
% end


