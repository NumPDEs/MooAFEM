% OptimalLocalMGSolver (subclass of MGSolver) Solves linear equations
%   iteratively using high-order p-optimal multilevel local smoothing (additive
%   Schwarz/block Jacobi).
%   One iteration = Vcycle(0,1) (no pre-/one post-smoothing step).
%
% With the solver comes an 

classdef LocalMgLowOrderVcycle < MGSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end

    properties (Access=protected)
        blf
        hoFes
        loFes
        p1Matrix
        p1Smoother
        P
        intergridMatrix
        inclusionMatrix
        freeVertices
        freeVerticesOld
        changedPatches
        patchwiseA
    end
    
    %% event data
    properties (Access=private)
        listenerHandle
    end

    %% methods
    methods (Access=public)
        function obj = LocalMgLowOrderVcycle(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end
            obj = obj@MGSolver();
            
            assert(fes.finiteElement.order > 1, ...
                'LocalMgLowOrderVcycle only works for higher order finite elements.')
            assert(isempty(blf.b), ...
                'Multigrid solvers only tested for symmetric problems.')
            
            obj.nLevels = 0;
            
            mesh = fes.mesh;
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', ':');
            obj.hoFes = fes;
            obj.P = Prolongation.chooseFor(obj.loFes);
            obj.blf = blf;
            
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end

        function setupSystemMatrix(obj, A)
            obj.nLevels = obj.nLevels + 1;
            
            L = obj.nLevels;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.loFes);
            freeVerticesHigher = getFreeDofs(obj.hoFes);
            
            obj.p1Matrix{L} = assemble(obj.blf, obj.loFes);
            obj.p1Matrix{L} = obj.p1Matrix{L}(obj.freeVertices, obj.freeVertices);
            obj.p1Smoother{L} = full(diag(obj.p1Matrix{L})).^(-1);
            
            if L >= 2
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(obj.freeVertices));
                obj.patchwiseA = assemblePatchwise(obj.blf, obj.hoFes);
                obj.inclusionMatrix = spaceProlongation(obj.loFes, obj.hoFes);
                obj.inclusionMatrix = obj.inclusionMatrix(obj.freeVertices, freeVerticesHigher);
            end
            setupSystemMatrix@MGSolver(obj, A);
        end

        function setupRhs(obj, b, varargin)
            setupRhs@MGSolver(obj, b, varargin{:});
        end

        % Geometric MultiGrid
        function [Cx, algError2] = Vcycle(obj, res)
            assert(isequal(size(res, 1), size(obj.A, 1)), ...
                'Setup for multilevel iteration not complete!')
            
            % if there is only one coarse level: exact solve
            if obj.nLevels == 1
                Cx = obj.A \ res;
                algError2 = scalarProduct(Cx, obj.A*Cx);
                return
            end

            L = obj.nLevels;
            algError2 = zeros(1, size(res, 2)); % built-in estimator of the algebraic error
            
            % descending cascade in P1: no smoothing
            residual{L} = obj.inclusionMatrix * res; %interpolateFreeData(res, obj.hoFes, obj.loFes);
            for k = L:-1:2
                residual{k-1}  = obj.intergridMatrix{k}'*residual{k};
            end

            % exact solve on coarsest level to compute accumulative lifting
            % of the residual (sigma)
            sigma = obj.p1Matrix{1} \ residual{1};
            algError2 = algError2 + scalarProduct(sigma, obj.p1Matrix{1}*sigma);

            % ascending cascade in P1 AND local smoothing
            for k = 2:(L-1)
                sigma = obj.intergridMatrix{k}*sigma;
                uptres = residual{k} - obj.p1Matrix{k}*sigma; %updated residual 
                rho = p1LocalSmoothing(obj, k, uptres);
                [aeUpd, sUpd] = computeOptimalUpdate(obj.p1Matrix{k}, uptres, rho);
                sigma = sigma + sUpd;
                algError2 = algError2 + aeUpd;
            end
            
            % smoothing on finest level dependent on p
            sigma = obj.intergridMatrix{L}*sigma;
            sigma = obj.inclusionMatrix' * sigma; %interpolateFreeData(sigma, obj.loFes, obj.hoFes);
            uptres = res - obj.A*sigma;
            rho = hoGlobalSmoothing(obj, uptres);
            [aeUpd, sUpd] = computeOptimalUpdate(obj.A, uptres, rho);
            sigma = sigma + sUpd;
            algError2 = algError2 + aeUpd;

            Cx = sigma;
        end
    end
    
    methods (Access=protected)
        % callback function to be executed before mesh refinement:
        % -> patches change for all vertices adjacent to refined edges and
        %    all new vertices
        function getChangedPatches(obj, mesh, bisecData)
            nCOld = mesh.nCoordinates;
            nCNew = mesh.nCoordinates + nnz(bisecData.bisectedEdges) + ...
                bisecData.nRefinedElements'*bisecData.nInnerNodes;
            bisectedEdgeNodes = unique(mesh.edges(:,bisecData.bisectedEdges));
            obj.changedPatches{obj.nLevels+1} = false(nCNew, 1);
            idx = [unique(bisectedEdgeNodes); ((nCOld+1):nCNew)'];
            obj.changedPatches{obj.nLevels+1}(idx) = true;
        end
        
        function rho = p1LocalSmoothing(obj, k, res)
            idx = obj.changedPatches{k};
            rho = zeros(size(res));
            rho(idx,:) = obj.p1Smoother{k}(idx).*res(idx,:);
        end
        
        function rho = hoGlobalSmoothing(obj, res)
            rho = obj.patchwiseA \ res;
        end
    end
end

%% auxiliary functions

function inclusionMatrix = spaceProlongation(fromFes, toFes)
    %Create the prolongation matrix on the unit triangle
    unittriangle = Mesh.loadFromGeometry('unittriangle');
    fromFesUnitTriangle = FeSpace(unittriangle, ...
        HigherOrderH1Fe(fromFes.finiteElement.order));
    toFesUnitTriangle = FeSpace(unittriangle, ...
        HigherOrderH1Fe(toFes.finiteElement.order));
    Vertices = 1:getDofs(fromFesUnitTriangle).nDofs;
    unitTriangleInclusionMatrix = ...
        permute(interpolateData(eye(length(Vertices)),...
        fromFesUnitTriangle, toFesUnitTriangle), [2 1]);
    % Get local to global map in unit triangle
    
    unitTriangleInclusionMatrix = unitTriangleInclusionMatrix( ...
        getDofs(fromFesUnitTriangle).element2Dofs, ...
        getDofs(toFesUnitTriangle).element2Dofs);
    % Get dofs of the finite element spaces on the meshes
    fromFesDofs = getDofs(fromFes);
    toFesDofs = getDofs(toFes);
    
    % Copy the matrix from the unit triangle for each element
    mat = repmat(unitTriangleInclusionMatrix, ...
        [1, 1, fromFes.mesh.nElements]);
    
    % Create index sets for the accumarray
    temp = fromFesDofs.element2Dofs(:, :);
    temp2 = toFesDofs.element2Dofs(:, :);
    I = repmat(permute(temp, [1 3 2]), [1, size(temp2, 1), 1]);
    J = repmat(permute(temp2, [3 1 2]), [size(temp, 1), 1, 1]);
    ind = [I(:), J(:)];
    
    inclusionMatrix = ...
        accumarray(ind, mat(:), [], @mean, [], true);
end

function interpolatedData = interpolateData(data, fromFes, toFes)
    Dofs = 1:getDofs(toFes).nDofs;
    nComponents = size(data, 2);
    interpolatedData = zeros(numel(Dofs), nComponents);
    feFunctionWrapper = FeFunction(fromFes);
    for k = 1:nComponents
        feFunctionWrapper.setData(data(:,k));
        wholeData = nodalInterpolation(feFunctionWrapper, toFes);
        interpolatedData(:,k) = wholeData(Dofs);
    end
end

% error correction with optimal stepsize 
function [etaUpdate, sigmaUpdate] = computeOptimalUpdate(A, res, rho)
    rhoArho = scalarProduct(rho, A*rho);
    if max(abs(rho)) < eps
        lambda = 1; 
    else
        lambda = scalarProduct(res, rho) ./ rhoArho;
    end
    sigmaUpdate = lambda.*rho;
    etaUpdate = rhoArho.*(lambda.^2);

    if any(lambda > 3)
       warning('MG step-sizes no longer bound by d+1. Optimality of step size cannot be guaranteed!')
    end
end

function a = scalarProduct(x, y)
    a = sum(x.*y, 1);
end
