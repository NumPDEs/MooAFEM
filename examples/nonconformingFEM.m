% *************************************************************************

% *************************************************************************

fprintf("#\n# Adaptive NCFEM - Mesh and solution\n#\n")

%% parameters
nDofsMax = 10000;
theta = 1;
lambda = 1e-1;

%% main run

% setup geometry & spaces
mesh = Mesh.loadFromGeometry('unitsquare');
% mesh.refineUniform();
fes = FeSpace(mesh, LowestOrderCRFe(), 'dirichlet', ':');
u = FeFunction(fes);
u.setData(0);

% set problem data
blf = BilinearForm();
lf = LinearForm();

blf.a = Constant(mesh, 1);

lf.f = Constant(mesh, 1);

% set up leveldata for output
leveldata = LevelData('results');
leveldata.metaData("problem") = "Poisson";
leveldata.metaData("domain") = "square";
leveldata.metaData("method") = "CR";
leveldata.metaData("identifier") = "S_nDofsMax" + num2str(nDofsMax) + ...
                                   "_theta" + num2str(fix(100*theta)) + ...
                                   "_lambda0";

%% adaptive loop
meshSufficientlyFine = false;  % initialize break condition
ell = 0;  % no of level
cost = 0;
runtime = 0;  % cumulative runtime
while ~meshSufficientlyFine
    %% setup
    ell = ell + 1;

    % assemble linear system
    A = assemble(blf, fes);
    F = assemble(lf, fes);
    freeDofs = getFreeDofs(fes);
    nDofs = getDofs(fes).nDofs;
    cost = cost + nDofs;

    % include Dirichlet boundary conditions
    % [rhs, uDirData] = setDirichletData(uExact_fct, fes, A, F);
    % u.setFixedData(uDirData);
    
    %% exact algebraic solution
    u.setFreeData(A(freeDofs,freeDofs) \ F(freeDofs));
    % eta2 = estimate(blf, lf, u);
    % estimator = sqrt(sum(eta2));
    estimator = nan;
    
    % store level-oriented data
    leveldata.append('ell', int32(ell), ...
                    'nDofs', int32(nDofs), ...
                    'nElem', int32(mesh.nElements), ...
                    'cost', int32(cost), ...
                    'estimator', estimator, ...
                    'case', 'R');
    
    %% stoping criterion
    meshSufficientlyFine = (nDofs > nDofsMax);
    
    %% refine mesh
    if ~meshSufficientlyFine
        if theta == 1
            mesh.refineUniform();
        else
            error('Not implemented yet.')
            % marked = markDoerflerBinning(eta2, theta);
            % mesh.refineLocally(marked, 'NVB');
        end
    end
    
    %% commandline output in case of mesh refinement
    leveldata.printLevel();
end


%% save data to files
leveldata.saveToTable();

%% plot mesh
h = figure(1);
cla;
plot(u.fes.mesh);
axis equal tight;
grid off;
ax = gca;
ax.Children.FaceAlpha = 0;
% saveFigure(h, [], "Fig1_Kellogg_problem_mesh");

%% production
% Remove axis for plot in tikz axis
% axis off;
% saveFigure(h, leveldata.foldername, leveldata.filename + '_mesh', ...
%            [0 0 5 5], false);

%% plot solution
h = figure(2);
cla;
plot(u);
grid off;
view(110, 30);
ax = gca;
ax.Children.EdgeColor = 'none';
% saveFigure(h, [], "Fig1_Kellogg_problem_solution");

%% production
% Remove axis for plot in tikz axis
% axis off;
% saveFigure(h, leveldata.foldername, leveldata.filename + '_solution', [], false);



function indicators = estimate(blf, lf, u)
    
    % abbreviate data structures
    fes = u.fes;
    mesh =  fes.mesh;
    trafo = getAffineTransformation(mesh);
    
    % compute volume residual element-wise
    res_handle = @(a, v,f) (a .* (v(1,:,:) + v(4,:,:)) + f).^2;
    res = CompositeFunction(res_handle, blf.a, Hessian(u), lf.f);
    qr = QuadratureRule.ofOrder(2);
    volumeResidual = integrateElement(res, qr);

    % compute edge residual edge-wise
    qr = QuadratureRule.ofOrder(1, '1D');
    f = CompositeFunction(@(p, a) a .* p, Gradient(u), blf.a);
    dirichletBnd = getCombinedBndEdges(mesh, fes.bnd.dirichlet);
    edgeResidual = integrateTangentialJump(f, qr, ...
        @(j) zeros(size(j)), {}, dirichletBnd, ...    % no jump on dirichlet edges
        @(j) j.^2, {}, ':');                          % tangential jump on inner edges
    
    % combine the resdiuals suitably
    hT = sqrt(trafo.area);
    indicators = ... hT.^2 .* volumeResidual + ...
        hT .* sum(edgeResidual(mesh.element2edges), Dim.Vector);
end