%% paramters
nElemMax = 1e5;
theta = 0.5;
lambda = 0.1;
deltaZ = 0.1;
linearization = 'kacanov';
geometry = 'Lshape';
p = 3;

%% initialisation
[nElem, eta, E] = deal(zeros(1,1000));

%% setup geometry & spaces
mesh = Mesh.loadFromGeometry(geometry);
mesh.refineUniform(2);
fes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', ':');
u = FeFunction(fes);

%% set problem data given in [Heid, Praetorius, Wihler; 2021]
Du = Gradient(u);
g = Constant(mesh, 1);
[blf, lf] = setProblemData(mesh, fes, Du, g, linearization, p);
eDensity = CompositeFunction(@(t) muIntegral(vectorProduct(t, t),p), Du);
u.setFreeData(rand(1,length(getFreeDofs(u.fes))));

%% nested iteration
P = LoMeshProlongation(fes);

%% adaptive algorithm
i = 1;
Eold = 0;
qr = QuadratureRule.ofOrder(1); % tweak?
while 1
    % the volume residual for lowest order is independent of the approximate
    % solution, so it can be computed once
    vol = estimateVolume(g);
        
    %% inner iteration to solve the nonlinear problem
    nIterations = 1;
    while true
        %% solve linearized problem
        A = assemble(blf, fes);
        F = assemble(lf, fes);
        updateDataU(u, deltaZ, A, F, linearization);
        
        %% compute remainder of error estimator and energy
        edge = estimateEdge(u, mesh, p);
        eta2 = combineEstimators(vol, edge, mesh);
        eta(i) = sqrt(sum(eta2));
        nElem(i) = mesh.nElements;
            
        if i > 1, Eold = E(i-1); end
        E(i) = sum(integrateElement(eDensity, qr)) + u.data * F;
            
        %% stopping criterion of the iterative solver
        if sqrt(E(i) - Eold) <= lambda * eta(i) % norm(e) <= fehler
            break
        end
        nIterations = nIterations + 1;
        i = i+1;
    end

    %% stoping criterion for the AFEM algorithm
    printLogMessage('number of elements: %d, iterations: %d, estimator: %.2e', ...
        nElem(i), nIterations, eta(i));
    if nElem(i) > nElemMax
        break
    end

    %% refine mesh
    marked = markDoerflerSorting(eta2, theta);
    mesh.refineLocally(marked, 'NVB');
    u.setData(prolongate(P, u));
    i = i+1;
end

%% plot convergence rates
plotData(nElem, 'number of elements', eta, 'error estimator', linearization,p)

%% helper function for plotting
function plotData(xdata, xlab, ydata, ylab, name,p)
    figure()
    for k = 1:size(xdata,1)
        idx = find(xdata(k,:) > 0);
        loglog(xdata(k,idx), ydata(k,idx), '-o', 'LineWidth', 2, 'DisplayName', strcat(name,', p=', num2str(p)));
        hold on
    end
    x = xdata(1,idx) / xdata(1,1);
    loglog(xdata(1,1)*x, ydata(1,1)*x.^(-1/2), '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', '\alpha = 1/2')
    legend
    xlabel(xlab)
    ylabel(ylab)
    title([ylab, ' over ', xlab])
end

%% local function for problem data
function [blf, lf] = setProblemData(mesh, fes, Du, g, linearization, p)
    blf = BilinearForm();
    lf = LinearForm();

    switch linearization
        case "zarantonello"
            blf.a = Constant(mesh, 1);
            lf.f = g;
            lf.fvec = CompositeFunction(@(t) -mu(vectorProduct(t, t),p) .* t, Du);
        case "kacanov"
            blf.a = CompositeFunction(@(t) mu(vectorProduct(t, t),p), Du);
            lf.f = g;
            lf.fvec = [];
        case "newton"
            blf.a = CompositeFunction(@(t) mu(vectorProduct(t, t),p) .* [1;0;0;1] ...
                + 2*muDerivative(vectorProduct(t, t),p).*vectorProduct(t, t, [2,1], [2,1]'), Du);
            lf.f = g;
            lf.fvec = CompositeFunction(@(t) -mu(vectorProduct(t, t),p) .* t, Du);
    end

    [blf.qra, lf.qrf, lf.qrfvec] = deal(QuadratureRule.ofOrder(1));
end

function updateDataU(u, deltaZ, A, F, linearization)
    freeDofs = getFreeDofs(u.fes);
    v = FeFunction(u.fes);
    switch linearization
        case "zarantonello"
            v.setFreeData(A(freeDofs,freeDofs) \ F(freeDofs));
            u.setData(u.data + deltaZ*v.data);
        case "kacanov"
            %u.setData(A \ F);
            u.setFreeData(A(freeDofs,freeDofs) \ F(freeDofs));
        case "newton"
            v.setFreeData(A(freeDofs,freeDofs) \ F(freeDofs));
            u.setData(u.data + v.data);
    end
end

%% local function for residual a posteriori error estimation
% \eta(T)^2 = h_T^2 * || g ||_{L^2(T)}^2
%               + h_T * || [[ mu(|Du|^2)*Du * n]] ||_{L^2(E) \cap \Omega}^2
function vol = estimateVolume(g)
    qr = QuadratureRule.ofOrder(1);
    vol = integrateElement(g, qr);
end
   
function edge = estimateEdge(u, mesh, p)
    qr = QuadratureRule.ofOrder(1, '1D');
    dirichlet = vertcat(mesh.boundaries{:});
    
    f = CompositeFunction(@(t) mu(vectorProduct(t, t),p).*t, Gradient(u));
    edge = integrateNormalJump(f, qr, ...
        @(j) zeros(size(j)), {}, dirichlet, ...     % no jump on dirichlet edges
        @(j) j.^2, {}, ':');                        % normal jump on inner edges
    
end

function indicators = combineEstimators(vol, edge, mesh)
    trafo = getAffineTransformation(mesh);
    hT = sqrt(trafo.area);
    indicators = hT.^2.*vol + hT.*sum(edge(mesh.element2edges), Dim.Vector);
end

%% nonlinearities
function y = mu(t,p)
    y = t.^(p/2 - 1);
end

function y = muIntegral(t,p)
    y = p/4 * t.^(p/2);
end

function y = muDerivative(t,p)
    y = (p/2 - 1) * t.^(p/2 - 2);
end
