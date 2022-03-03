% ******************************************************************************
% Adaptive FEM algorithm with known solution and convergence rates.
% This was used to generate data for Figure 7.
% ******************************************************************************

%% paramters
nDofsMax = 1e6;
theta = 0.5;
p = 1;

%% initialization
[nElem, nDofs, errEst] = deal(zeros(1, 1000));
time = zeros(6,1000);

%% setup geometry & spaces
printLogMessage('*** p = %d ***', p)
mesh = Mesh.loadFromGeometry('Lshape');
fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', 1);
u = FeFunction(fes);
uex = FeFunction(fes);

%% set problem data for -\Delta u = 1 on L-shape
blf = BilinearForm(fes);
lf = LinearForm(fes);

blf.a = Constant(mesh, 1);
blf.qra = QuadratureRule.ofOrder(max(2*p-2, 1));

lf.neumann = MeshFunction(mesh, @exactSolutionNeumannData);
lf.qrNeumann = QuadratureRule.ofOrder(2*p, '1D');
lf.bndNeumann = 2;

%% adaptive loop
i = 1;
tic
while 1
    %% assemble & solve FEM system
    A = assemble(blf);
    time(1,i) = toc;
    F = assemble(lf);
    time(2,i) = toc;
    free = getFreeDofs(fes);
    u.setFreeData(A(free,free) \ F(free));
    time(3,i) = toc;

    %% estimate error and store data
    eta2 = estimate(blf, lf, u);
    errEst(i) = sqrt(sum(eta2));
    time(4,i) = toc;
    nDofs(i) = getDofs(fes).nDofs;
    nElem(i) = mesh.nElements;
    printLogMessage('number of dofs: %d, estimator: %.2e', nDofs(i), errEst(i));

    %% refine mesh
    marked = markDoerflerSorting(eta2, theta);
    time(5,i) = toc;
    mesh.refineLocally(marked, 'NVB');
    time(6,i) = toc;
    
    %% stoping criterion
    if nDofs(i) > nDofsMax
        break
    end
    i = i+1;
end

%% store fine grained computation times
idx = find(nDofs > 0);
time = time(:,idx);
time = reshape([time(1,1); diff(time(:))], [], length(idx));
fileID = fopen(['timingP', num2str(p), '.dat'], 'w');
fprintf(fileID, 'nElements,nDofs,tAssembleA,tAssembleF,tSolve,tEstimate,tMark,tRefine,tCombined\n');
fprintf(fileID, '%d,%d,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n', ...
    [nElem(idx); nDofs(idx); time; sum(time, 1)]);
fclose(fileID);

%% local function for residual a posteriori error estimation
function indicators = estimate(blf, lf, u)
    p = blf.fes.finiteElement.order;
    mesh =  blf.fes.mesh;
    
    % compute volume residual element-wise
    if p == 1
        volumeRes = 0;
    else
        f = CompositeFunction(@(D2u) (D2u(1,:,:) + D2u(4,:,:)).^2, Hessian(u));
        qr = QuadratureRule.ofOrder(max(2*(p-2), 1));
        volumeRes = integrateElement(f, qr);
    end
    
    % compute edge residual edge-wise
    qr = QuadratureRule.ofOrder(p, '1D');
    edgeRes = integrateNormalJump(Gradient(u), qr, ...
        @(j) zeros(size(j)), {}, mesh.boundaries{1}, ...    % no jump on dirichlet edges
        @(j,phi) j-phi, {lf.neumann}, mesh.boundaries{2},...% neumann jump on neumann edges
        @(j) j.^2, {}, ':');                                % normal jump on inner edges
    
    % combine the resdiuals suitably
    hT = sqrt(getAffineTransformation(mesh).area);
    indicators = hT.^2 .* volumeRes + ...
        hT .* sum(edgeRes(mesh.element2edges), Dim.Vector);
end

%% exact solution for Lshape
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
