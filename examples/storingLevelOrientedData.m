function leveldata = storingLevelOrientedData(doPlots)
%%MAINLEVELDATA Example of an adaptive FEM algorithm with storing and
%displaying level-oriented dataa

    %% proceed optional input
    if nargin < 1
        doPlots = true;
    end

    %% paramters
    nEmax = 1e4;
    theta = 0.5;
    
    %% initialization
    mesh = Mesh.loadFromGeometry('Lshape');
    fes = FeSpace(mesh, LowestOrderH1Fe());
    u = FeFunction(fes);

    %% initialize output data structure
    pathToStorage = 'results';
    leveldata = LevelData(pathToStorage);
    leveldata.problem = 'Poisson';
    leveldata.domain = 'Lshape';
    leveldata.method = 'S1';
    leveldata.identifier = 'example';

    %% problem data
    blf = BilinearForm();
    lf = LinearForm();
    blf.a = Constant(mesh, [1;0;0;1]);
    blf.b = Constant(mesh, [0;0]);
    blf.c = Constant(mesh, 0);
    lf.f = Constant(mesh, 1);
    lf.fvec = Constant(mesh, [0;0]);
    
    % reference value for energy norm of exact solution
    normExactSquared = 0.2140750232;

    %% adaptive loop
    meshSufficientlyFine = false;
    while ~meshSufficientlyFine
        %% compute data for rhs coefficient and assemble forms
        A = assemble(blf, fes);
        F = assemble(lf, fes);
        
        %% solve FEM system
        freeDofs = getFreeDofs(fes);
        tic;
        u.setFreeData(A(freeDofs,freeDofs) \ F(freeDofs));
        runtimeSolve = toc;

        %% approximate error
        err = sqrt(normExactSquared - u.data * A * u.data');

        %% compute refinement indicators
        tic;
        eta2 = estimate(blf, lf, u);
        eta = sqrt(sum(eta2));
        runtimeEstimate = toc;

        %% compute efficiency index
        eff = err / eta;

        %% storing error quantities and general data
        leveldata.append('ndof', uint32(length(freeDofs)), ...
                         'eta', eta, ...
                         'err', err, ...
                         'coordinates', mesh.coordinates, ...
                         'elements', mesh.elements);

        %% storing results of timing
        leveldata.setTime(leveldata.nLevel, ...
                          'runtimeSolve', runtimeSolve, ...
                          'runtimeEstimate', runtimeEstimate);

        %% storing absolute values
        leveldata.setAbsolute(leveldata.nLevel, 'eff', eff);

        %% print information on current level to command line
        leveldata.printLevel();

        %% save intermediate results to file 
        % (prevents data loss in case of crashes)
        leveldata.saveToFile();
      
        %% stoping criterion
        meshSufficientlyFine = (mesh.nElements > nEmax);
        
        %% refine mesh
        if ~meshSufficientlyFine
            marked = markDoerflerBinning(eta2, theta);
            mesh.refineLocally(marked, 'NVB');
        end
    end

    %% post-processing variables
    leveldata.setTime(':', 'runtimeTotal', ...
                      sum(leveldata.get(':', 'runtimeSolve', 'runtimeEstimate'), 2));

    %% plot results
    if doPlots
        % plot error quantities (in general, converging to zeros)
        figure();
        leveldata.plot('ndof');
    
        % plot runtimes (in general, diverging to infinity)
        figure();
        leveldata.plotTime('ndof');
    
        % plot absolute values in semi-logarithmic plot
        figure();
        leveldata.plotAbsolute('ndof');
    
        % plot of all triangulations
        figure();
        leveldata.plotTriangulation();
    
        %% export error plot to file
        leveldata.plotToFile('ndof');

        %% export refinement video
        % (generates video of about 110 MB)
        %leveldata.plotTriangulationToFile();
    end

    %% save all values to comma-separated table
    leveldata.saveToTable();
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
