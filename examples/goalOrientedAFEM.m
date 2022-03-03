% ******************************************************************************
% Example of a GOAFEM algorithm taken from [Mommer, Stevenson; 2009].
% ******************************************************************************

%% paramters
nDofsMax = 1e4;
theta = 0.5;
pmax = 3;

%% initialization for every polynomial degree
[nElem, nDofs, goalErrEst] = deal(zeros(pmax, 1000));
for p = 1:pmax
    %% setup geometry & spaces
    printLogMessage('*** p = %d (of %d) ***', p, pmax)
    mesh = Mesh.loadFromGeometry('unitsquare');
    mesh.refineUniform(1, 'RGB');
    ncFes = FeSpace(mesh, LowestOrderL2Fe);
    v = FeFunction(ncFes);
    w = FeFunction(ncFes);
    fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', 1);
    u = FeFunction(fes);
    z = FeFunction(fes);
    
    %% set problem data given in [Mommer, Stevenson; 2009]
    blf = BilinearForm(fes);
    blf.a = Constant(mesh, 1);
    blf.qra = QuadratureRule.ofOrder(max(2*p-2, 1));
    
    chiT1 = MeshFunction(mesh, @(x) sum(x, Dim.Vector) < 1/2);
    v.setData(nodalInterpolation(chiT1, ncFes));
    lfF = LinearForm(fes);
    lfF.fvec = CompositeFunction(@(v) [v;zeros(size(v))], v);
    lfF.qrfvec = QuadratureRule.ofOrder(max(p-1, 1));
    
    chiT2 = MeshFunction(mesh, @(x) sum(x, Dim.Vector) > 3/2);
    w.setData(nodalInterpolation(chiT2, ncFes));
    lfG = LinearForm(fes);
    lfG.fvec = CompositeFunction(@(w) [-w;zeros(size(w))], w);
    lfG.qrfvec = QuadratureRule.ofOrder(max(p-1, 1));
    
    %% set up lifting operators for rhs FEM-data
    Pv = LoL2Prolongation(v);
    Pw = LoL2Prolongation(w);

    %% adaptive loop
    i = 1;
    while 1
        %% assemble & solve FEM system
        A = assemble(blf);
        rhs = [assemble(lfF), assemble(lfG)];
        freeDofs = getFreeDofs(fes);
        uz = A(freeDofs,freeDofs) \ rhs(freeDofs,:);
        u.setFreeData(uz(:,1));
        z.setFreeData(uz(:,2));

        %% estimate error and store data
        eta2 = estimate(blf, lfF, u);
        zeta2 = estimate(blf, lfG, z);
        nDofs(p,i) = getDofs(fes).nDofs;
        nElem(p,i) = mesh.nElements;
        goalErrEst(p,i) = sqrt(sum(eta2)*sum(zeta2));
        printLogMessage('number of dofs: %d, estimator: %.2e', nDofs(p,i), goalErrEst(p,i));

        %% stoping criterion
        if nDofs(p,i) > nDofsMax
            break
        end

        %% refine mesh
        marked = markGoafemMS(eta2, zeta2, theta);
        mesh.refineLocally(marked, 'NVB');
        i = i+1;
        v.setData(Pv.interpolatedData);
        w.setData(Pw.interpolatedData);
    end
end

%% plot convergence rates
pmax = size(nDofs, 1);
figure()
for p = 1:pmax
    idx = find(nDofs(p,:) > 0);
    loglog(nDofs(p,idx), goalErrEst(p,idx), '-o', 'LineWidth', 2, 'DisplayName', ['p=',num2str(p)]);
    hold on
end
for p = unique([1,pmax])
    x = nDofs(p,idx) / nDofs(p,1);
    loglog(nDofs(p,1)*x, goalErrEst(p,1)*x.^(-p), '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', ['\alpha = ', num2str(p)])
end
legend
xlabel('number of dofs')
ylabel('goal error estimator')
title(['goal error estimator over number of dofs'])

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
        volumeResidual = 0;
    else
        f = CompositeFunction(@(v) (v(1,:,:) + v(4,:,:)).^2, Hessian(u));
        qr = QuadratureRule.ofOrder(max(2*p-4, 1));
        volumeResidual = integrateElement(f, qr);
    end
    
    % compute edge residual edge-wise
    qr = QuadratureRule.ofOrder(max(p-1, 1), '1D');
    f = CompositeFunction(@(p,fvec) p-fvec, Gradient(u), lf.fvec);
    edgeResidual = integrateNormalJump(f, qr, ...
        @(j) zeros(size(j)), {}, mesh.boundaries{1}, ...    % no jump on dirichlet edges
        @(j) j.^2, {}, ':');                                % normal jump on inner edges
    
    % combine the resdiuals suitably
    hT = sqrt(trafo.area);
    indicators = hT.^2 .* volumeResidual + ...
        hT .* sum(edgeResidual(mesh.element2edges), Dim.Vector);
end
