% ******************************************************************************
% Adaptive FEM algorithm with known solution and convergence rates for higher
% order finite elements.
% ******************************************************************************

%% parameters
nDofsMax = 1e4;
theta = 0.5;
pmax = 5;
mu = 0.1;
directSol = 0;

%for plotting
allAlgErr = [];
allAlgEst = [];
allContr = [];
time = [];


%% initialization for every polynomial degree
[nElem, nDofs, errEst, h1Err] = deal(zeros(pmax, 1000));
for p = 1:pmax
    %% setup geometry & spaces
    printLogMessage('*** p = %d (of %d) ***', p, pmax)
    mesh = Mesh.loadFromGeometry('Lshape');
    fes = FeSpace(mesh, HigherOrderH1Fe(p), 'dirichlet', ':');
    u = FeFunction(fes);
    u.setData(0);
    uex = FeFunction(fes);
    
    %% set problem data for -\Delta u = 1 on L-shape
    blf = BilinearForm();
    lf = LinearForm();
    
    blf.a = Constant(mesh, 1);
    blf.qra = QuadratureRule.ofOrder(max(2*p-2, 1));
    
    
    lf.f = MeshFunction(mesh, @(x) 1);
    lf.qrf = QuadratureRule.ofOrder(2*p);
    
    %% set up solver and operator for nested iteration
    P = FeProlongation(fes);
    solver = pLoc_MG(fes, blf);
    solver.tol = 1e-8;
    solver.maxIter = 100;
    
    %% adaptive loop
    i = 0;
    ell = 0;
    meshSufficientlyFine = false;
    
    
    while 1
        
        ell = ell + 1; i = i+1;
        
        freeDofs = getFreeDofs(fes);
        A = assemble(blf, fes);
        rhs = assemble(lf, fes);
        
        
        
        
        
        if ~directSol
            
            %freeDofs = getFreeDofs(fes);
            p1freedofs = freeDofs(freeDofs<=mesh.nCoordinates);
            mesh.freeVert{end+1} = p1freedofs;
            %A = assemble(blf);
            %rhs = assemble(lfF);
            u0 = u.data';
            iter = 0; algEta2 = 10;
            %setup : interpolation matrices; patches for high-order
            solver.setupLinearSystem(A, rhs(freeDofs), u0(freeDofs));
            
            if ell == 1
                nIterPrimal(ell) = 1;
                solver.step();
                algEta2 = solver.algEst2;
                
            else
                
                %exact solve
                x_ex = A(freeDofs,freeDofs)\rhs(freeDofs);
                alg_error2 = (x_ex - solver.x)'*A(freeDofs,freeDofs)*(x_ex - solver.x);
                
                tic;
                while (algEta2 >= mu^2*sum(eta2)) %adaptive stopping criterion with squares
                    
                    iter = iter+1;
                    solver.step();
                    oldalg_error2 = alg_error2;
                    alg_error2 = (x_ex - solver.x)'*A(freeDofs,freeDofs)*(x_ex - solver.x);
                    algEta2 = solver.algEst2;
                    
                    allAlgErr = [allAlgErr;sqrt(oldalg_error2)];
                    allAlgEst = [allAlgEst;sqrt(algEta2)];
                    allContr = [allContr;sqrt(alg_error2)/sqrt(oldalg_error2)];
                end
                
                t = toc;
                time = [time,t];
                nIterPrimal(ell) = solver.iterationCount(1);
            end
            
            u.setFreeData(solver.x);
            
        else
            
            %direct solve
            x_ex = A(freeDofs,freeDofs)\rhs(freeDofs);
            u.setFreeData(x_ex);
            nIterPrimal(ell) = 1;
        end
        
        
        
        
        
        
        
        
        
        %% estimate error and store data
        eta2 = estimate(blf, lf, u);
        nDofs(p,i) = getDofs(fes).nDofs;
        errEst(p,i) = sqrt(sum(eta2));
        nElem(p,i) = mesh.nElements;
        printLogMessage('number of dofs: %d, estimator: %.2e', nDofs(p,i), errEst(p,i));
        
        %% stoping criterion
        if nDofs(p,i) > nDofsMax
            break
        end
        
        %% refine mesh
        %marked = markDoerflerSorting(eta2, theta);
        marked = markDoerflerBinning(eta2, theta);
        
        meshtemp = clone(mesh);
        meshtemp.intergridMatrix(marked, 'NVB');
        mesh.refineLocally(marked, 'NVB');
        mesh.locVert{end+1} = meshtemp.locVert{end};
        
        u.setData(prolongate(P, u));
    end
end

%% plot convergence rates
plotData(nDofs, 'number of dofs', errEst, 'error estimator')
%plotData(nDofs, 'number of dofs', h1Err, 'H^1 error')

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
