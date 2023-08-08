% Example of an adaptive FEM algorithm with storing and displaying
% level-oriented data

function leveldata = storingLevelOrientedData(doPlots, doStore)
    arguments
        doPlots (1,1) logical = true
        doStore (1,1) logical = false
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
    leveldata.metaData('problem') = 'Poisson';
    leveldata.metaData('domain') = 'Lshape';
    leveldata.metaData('method') = 'S1';
    leveldata.metaData('identifier') = 'example';
    
    %% problem data
    blf = BilinearForm();
    lf = LinearForm();
    blf.a = Constant(mesh, 1);
    lf.f = Constant(mesh, 1);
    
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
        if doStore
            leveldata.saveToFile();
        end
    
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
    end

    %% additional visualization/storage capabilities
    if doStore
        % export error plot to file
        leveldata.plotToFile('ndof');

        % save all values to comma-separated table
        leveldata.saveToTable();

        % plot of all triangulations one after another
        figure();
        leveldata.plotTriangulation();
    
        % export refinement video (generates video of about 110 MB)
        leveldata.plotTriangulationToFile();
    end
end



%% local function for residual a posteriori error estimation
% \eta(T)^2 = h_T^2 * || \Delta u + f ||_{L^2(T)}^2 + h_T * || [Du] ||_{L^2(E)}^2
function indicators = estimate(~, ~, u)
    fes = u.fes;
    mesh =  fes.mesh;
    
    %% compute volume residual element-wise
    f = CompositeFunction(@(Du) vectorProduct(Du, Du), Gradient(u));
    qrTri = QuadratureRule.ofOrder(1);
    volumeRes = integrateElement(f, qrTri);
    
    %% compute edge residual edge-wise
    qrEdge = QuadratureRule.ofOrder(1, '1D');
    edgeRes = integrateNormalJump(Gradient(u), qrEdge, @(j) j.^2, {}, ':');
    dirichlet = getCombinedBndEdges(mesh, fes.bnd.dirichlet);
    edgeRes(dirichlet) = 0;
    
    %% combine the resdiuals suitably
    hT = sqrt(getAffineTransformation(mesh).area);
    indicators = hT.^2 .* volumeRes + ...
        hT .* sum(edgeRes(mesh.element2edges), Dim.Vector);
end
