% ******************************************************************************
% Example from Section 5.3 of [Heid, Praetorius, Wihler; 2021].
% ******************************************************************************

%% paramters
nElemMax = 1e4;
theta = 0.5;
lambda = 0.1;
deltaZ = 0.1;
linearizations = ["zarantonello", "kacanov", "newton"];

%% initialization for every linearization method
[nElem, eta, E] = deal(zeros(length(linearizations), 1000));
for k = 1:length(linearizations)
    %% setup geometry & spaces
    printLogMessage('*** starting %s iteration', linearizations(k))
    mesh = Mesh.loadFromGeometry('Lshape');
    mesh.refineUniform(2);
    fes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', ':');
    u = FeFunction(fes);

    %% set problem data given in [Heid, Praetorius, Wihler; 2021]
    Du = Gradient(u);
    g = Constant(mesh, 1);
    [blf, lf] = setProblemData(mesh, Du, g, linearizations(k));
    eDensity = CompositeFunction(@(p) 1/2*muIntegral(vectorProduct(p, p)), Du);
    u.setData(0);    

    %% nested iteration
    P = LoFeProlongation(fes);

    %% adaptive algorithm
    ell = 0;
    Eold = 0;
    qr = QuadratureRule.ofOrder(1);
    meshSufficientlyFine = false;
    while ~meshSufficientlyFine
        % the volume residual for lowest order is independent of the approximate
        % solution, so it can be computed once
        vol = estimateVolume(g);
        
        %% inner iteration to solve the nonlinear problem
        nIterations = 0;
        solverIsConverged = false;
        while ~solverIsConverged
            %% solve linearized problem
            ell = ell + 1;
            nIterations = nIterations + 1;
            
            A = assemble(blf, fes);
            F = assemble(lf, fes);
            updateDataU(u, deltaZ, A, F, linearizations(k));
            
            %% compute remainder of error estimator and energy
            edge = estimateEdge(u, mesh);
            eta2 = combineEstimators(vol, edge, mesh);
            eta(k,ell) = sqrt(sum(eta2));
            nElem(k,ell) = mesh.nElements;
            E(k,ell) = sum(integrateElement(eDensity, qr)) + u.data * F;
            
            %% stopping criterion of the iterative solver
            solverIsConverged = (sqrt(E(k,ell) - Eold) <= lambda * eta(k,ell));
            Eold = E(k,ell);
        end

        %% stoping criterion for the AFEM algorithm
        printLogMessage('number of elements: %d, iterations: %d, estimator: %.2e', ...
            nElem(k,ell), nIterations, eta(k,ell));
        meshSufficientlyFine = (nElem(k,ell) > nElemMax);

        %% refine mesh
        if ~meshSufficientlyFine
            marked = markDoerflerSorting(eta2, theta);
            mesh.refineLocally(marked, 'NVB');
            u.setData(prolongate(P, u));
        end
    end
end

%% plot convergence rates
plotData(nElem, 'number of elements', eta, 'error estimator', linearizations)

%% helper function for plotting
function plotData(xdata, xlab, ydata, ylab, names)
    figure()
    for k = 1:size(xdata,1)
        idx = find(xdata(k,:) > 0);
        loglog(xdata(k,idx), ydata(k,idx), '-o', 'LineWidth', 2, 'DisplayName', names(k));
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
function [blf, lf] = setProblemData(mesh, Du, g, linearization)
    blf = BilinearForm();
    lf = LinearForm();

    switch linearization
        case "zarantonello"
            blf.a = Constant(mesh, 1);
            lf.f = g;
            lf.fvec = CompositeFunction(@(p) -mu(vectorProduct(p, p)) .* p, Du);
        case "kacanov"
            blf.a = CompositeFunction(@(p) mu(vectorProduct(p, p)), Du);
            lf.f = g;
            lf.fvec = [];
        case "newton"
            blf.a = CompositeFunction(@(p) mu(vectorProduct(p, p)) .* [1;0;0;1] ...
                + 2*muDerivative(vectorProduct(p, p)).*vectorProduct(p, p, [2,1], [2,1]'), Du);
            lf.f = g;
            lf.fvec = CompositeFunction(@(p) -mu(vectorProduct(p, p)) .* p, Du);
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
   
function edge = estimateEdge(u, mesh)
    qr = QuadratureRule.ofOrder(1, '1D');
    dirichlet = getCombinedBndEdges(mesh, u.fes.bnd.dirichlet);
    
    f = CompositeFunction(@(p) mu(vectorProduct(p,p)).*p, Gradient(u));
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
function y = mu(t)
    y = 1 + exp(-t);
end

function y = muIntegral(t)
    y = 1/2 * (1 + t - exp(-t));
end

function y = muDerivative(t)
    y = -exp(-t);
end
