% ******************************************************************************
% Adaptive FEM algorithm with known solution and convergence rates for higher
% order finite elements.
% ******************************************************************************

%% paramters
nDofsMax = 1e4;
theta = 0.5;
pmax = 8;

%% initialization for every polynomial degree
[nElem, nDofs, errEst, h1Err] = deal(zeros(pmax, 1000));
for p = 1:pmax
    %% setup geometry & spaces
    printLogMessage('*** p = %d (of %d) ***', p, pmax)
    mesh = Mesh.loadFromGeometry('Lshape');
    fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', 1, 'neumann', 2);
    u = FeFunction(fes);
    uex = FeFunction(fes);
    
    %% set problem data for -\Delta u = 1 on L-shape
    blf = BilinearForm();
    lf = LinearForm();
    
    blf.a = Constant(mesh, 1);
    blf.qra = QuadratureRule.ofOrder(max(2*p-2, 1));
    
    lf.neumann = MeshFunction(mesh, @exactSolutionNeumannData);
    lf.qrNeumann = QuadratureRule.ofOrder(2*p, '1D');

    %% adaptive loop
    ell = 1;
    meshSufficientlyFine = false;
    while ~meshSufficientlyFine
        %% assemble & solve FEM system
        ell = ell + 1;
        A = assemble(blf, fes);
        F = assemble(lf, fes);
        free = getFreeDofs(fes);
        u.setFreeData(A(free,free) \ F(free));
        uex.setData(nodalInterpolation(MeshFunction(mesh, @exactSolution), fes));

        %% estimate error and store data
        eta2 = estimate(blf, lf, u);
        nDofs(p,ell) = getDofs(fes).nDofs;
        errEst(p,ell) = sqrt(sum(eta2));
        deltaU = u.data - uex.data;
        h1Err(p,ell) = sqrt(deltaU * A * deltaU');
        nElem(p,ell) = mesh.nElements;
        printLogMessage('number of dofs: %d, estimator: %.2e', nDofs(p,ell), errEst(p,ell));

        %% stoping criterion
        meshSufficientlyFine = (nDofs(p,ell) > nDofsMax);

        %% refine mesh
        if ~meshSufficientlyFine
            marked = markDoerflerSorting(eta2, theta);
            mesh.refineLocally(marked, 'NVB');
        end
    end
end

%% plot convergence rates
plotData(nDofs, 'number of dofs', errEst, 'error estimator')
plotData(nDofs, 'number of dofs', h1Err, 'H^1 error')

%% helper function for plotting
function plotData(xdata, xlab, ydata, ylab)
    pmax = size(xdata, 1);
    figure()
    for p = 1:pmax
        idx = find(xdata(p,:) > 0);
        loglog(xdata(p,idx), ydata(p,idx), '-o', 'LineWidth', 2, 'DisplayName', ['p=',num2str(p)]);
        hold on
    end
    for p = unique([1,pmax])
        x = xdata(p,idx) / xdata(p,1);
        loglog(xdata(p,1)*x, ydata(p,1)*x.^(-p/2), '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', ['\alpha = ', num2str(p), '/2'])
    end
    legend
    xlabel(xlab)
    ylabel(ylab)
    title([ylab, ' over ', xlab])
end

%% local function for residual a posteriori error estimation
% \eta(T)^2 = h_T^2 * || \Delta u + f ||_{L^2(T)}^2
%               + h_T * || [[Du * n]] ||_{L^2(E) \cap \Omega}^2
%               + h_T * || Du * n - phi ||_{L^2(E) \cap \Gamma_N}^2
function indicators = estimate(~, lf, u)
    fes = u.fes;
    p = fes.finiteElement.order;
    mesh =  fes.mesh;
    
    % compute volume residual element-wise
    % For p=1, the diffusion term vanishes in the residual.
    if p == 1
        volumeRes = 0;
    else
        f = CompositeFunction(@(D2u) (D2u(1,:,:) + D2u(4,:,:)).^2, Hessian(u));
        qr = QuadratureRule.ofOrder(max(2*(p-2), 1));
        volumeRes = integrateElement(f, qr);
    end
    
    % compute edge residual edge-wise
    qr = QuadratureRule.ofOrder(p, '1D');
    dirichletBnd = getCombinedBndEdges(mesh, fes.bnd.dirichlet);
    neumannBnd = getCombinedBndEdges(mesh, fes.bnd.neumann);
    edgeRes = integrateNormalJump(Gradient(u), qr, ...
        @(j) zeros(size(j)), {}, dirichletBnd, ...          % no jump on dirichlet edges
        @(j,phi) j-phi, {lf.neumann}, neumannBnd,...        % neumann jump on neumann edges
        @(j) j.^2, {}, ':');                                % normal jump on inner edges
    
    % combine the resdiuals suitably
    hT = sqrt(getAffineTransformation(mesh).area);
    indicators = hT.^2 .* volumeRes + ...
        hT .* sum(edgeRes(mesh.element2edges), Dim.Vector);
end

%% local functions for exact solution
function y = exactSolution(x)
    [phi, r] = cart2pol(x(1,:,:), x(2,:,:));
    phi = phi + 2*pi*(phi < 0);
    y = r.^(2/3) .* sin(2/3 * phi);
end

function y = exactSolutionNeumannData(x)
    x1 = x(1,:,:);
    x2 = x(2,:,:);
    
    % determine boundary partsexactSolutionNeumannData
    right  = (x1 > 0) & (abs(x1) > abs(x2));
    left   = (x1 < 0) & (abs(x1) > abs(x2));
    top    = (x2 > 0) & (abs(x1) < abs(x2));
    bottom = (x2 < 0) & (abs(x1) < abs(x2));

    % compute d/dn uex
    [phi, r] = cart2pol(x(1,:,:), x(2,:,:));
    Cr = 2/3 * r.^(-4/3);
    Cphi = 2/3 * (phi + 2*pi*(phi < 0));
    dudx = Cr .* (x1.*sin(Cphi) - x2.*cos(Cphi));
    dudy = Cr .* (x2.*sin(Cphi) + x1.*cos(Cphi));

    y = zeros(size(x1));
    y(right)  = dudx(right);
    y(left)   =-dudx(left);
    y(top)    = dudy(top);
    y(bottom) =-dudy(bottom);
end
