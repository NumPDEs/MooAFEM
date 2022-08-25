% ******************************************************************************
% Example of an adaptive FEM algorithm.
% ******************************************************************************

%% paramters
nEmax = 1e4;
theta = 0.5;

%% initialization
mesh = Mesh.loadFromGeometry('Lshape');
fes = FeSpace(mesh, LowestOrderH1Fe());
ncFes = FeSpace(mesh, LowestOrderL2Fe());
midpoint = [1;1;1]/3;
u = FeFunction(fes);

%% problem data
blf = BilinearForm();
lf = LinearForm();
blf.a = Constant(mesh, [1;0;0;2]);
blf.b = Constant(mesh, [1;1]);
blf.c = MeshFunction(mesh, @(x) sum(x.^2, Dim.Vector));
lf.f = Constant(mesh, 1);
cutoff = MeshFunction(mesh, @(x) (diff(x, [], Dim.Vector) > 1));
w = FeFunction(ncFes);
lf.fvec = CompositeFunction(@(w) [w;-w], w);

%% adaptive loop
meshSufficientlyFine = false;
while ~meshSufficientlyFine
    %% compute data for rhs coefficient and assemble forms
    w.setData(eval(cutoff, midpoint));
    A = assemble(blf, fes);
    F = assemble(lf, fes);
    
    %% solve FEM system
    freeDofs = getFreeDofs(fes);
    u.setFreeData(A(freeDofs,freeDofs) \ F(freeDofs));
    
    %% compute refinement indicators
    eta2 = estimate(blf, lf, u);
    printLogMessage('number of elements: %d, estimator: %.2e', mesh.nElements, sqrt(sum(eta2)));
    
    %% stoping criterion
    meshSufficientlyFine = (mesh.nElements > nEmax);
    
    %% refine mesh
    if ~meshSufficientlyFine
        marked = markDoerflerBinning(eta2, theta);
        mesh.refineLocally(marked, 'NVB');
    end
end

%% plot solution and estimator
if mesh.nElements < 1e5
    plot(u);
    w.setData(log(eta2));
    plot(w);
end

%% local function for residual a posteriori error estimation
% This function shows how to build a residual error estimator with the provided
% tools. Not hiding this in its own class has the advantage that estimation is
% very flexible with respect to the involved terms.
% \eta(T)^2 = h_T^2 * || -div(a*Du) + b*Du + c*u - f + div(fvec) ||_{L^2(T)}^2
%               + h_T * || [ a*Du - fvec ] ||_{L^2(E)}^2
function indicators = estimate(blf, lf, u)
    fes = u.fes;
    mesh =  fes.mesh;
    
    %% compute volume residual element-wise
    % Here, one can optimize the estimator for the given situation. In the case
    % p=1 with constant diffusion, the diffusion term vanishes in the residual.
    % Also div(fvec) vanishes if fvec is element-wise constant.
    f = CompositeFunction(...
        @(b, c, f, u, Du) (vectorProduct(b, Du) + c .* u - f).^2, ...
        blf.b, blf.c, lf.f, u, Gradient(u));
    qrTri = QuadratureRule.ofOrder(4);
    volumeRes = integrateElement(f, qrTri);
    
    %% compute edge residual edge-wise
    % Here, none of the terms vanish, so they must be computed. This can be done
    % with a lower order than the element-residuals.
    f = CompositeFunction(...
        @(a, fvec, Du) vectorProduct(a, Du, [2,2], [2,1]) - fvec, ...
        blf.a, lf.fvec, Gradient(u));
    qrEdge = QuadratureRule.ofOrder(1, '1D');
    edgeRes = integrateNormalJump(f, qrEdge, @(j) j.^2, {}, ':');
    dirichlet = getCombinedBndEdges(mesh, fes.bnd.dirichlet);
    edgeRes(dirichlet) = 0;
    
    %% combine the resdiuals suitably
    hT = sqrt(getAffineTransformation(mesh).area);
    indicators = hT.^2 .* volumeRes + ...
        hT .* sum(edgeRes(mesh.element2edges), Dim.Vector);
end
