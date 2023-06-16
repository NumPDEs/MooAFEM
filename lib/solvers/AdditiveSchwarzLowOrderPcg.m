% AdditiveSchwartzPcg (subclass of PcgSolver) Solves linear equations
%   iteratively using the CG method with optimal multilevel additive
%   Schwartz preconditioner.

classdef AdditiveSchwarzLowOrderPcg < PcgSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end
    
    properties (Access=protected)
        blf
        hoFes
        loFes
        p1Acoarse
        p1Smoother
        P
        inclusionMatrix
        intergridMatrix
        freeVertices
        freeVerticesOld
        changedPatches
        patchwiseA
        patchwiseP1Matrix
    end

    properties (Access=private)
        listenerHandle
    end
    
    %% methods
    methods (Access=public)
        function obj = AdditiveSchwarzLowOrderPcg(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end
            obj = obj@PcgSolver();
            
            assert(fes.finiteElement.order > 1, ...
                'AdditiveSchwartzPcg only works for higher order finite elements.')
            assert(isempty(blf.b), ...
                'Additive Schwarz PCG solvers only tested for symmetric problems.')

            obj.nLevels = 0;

            mesh = fes.mesh;
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe, 'dirichlet', ':');
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

            p1Matrix = assemble(obj.blf, obj.loFes);
            p1Matrix = p1Matrix(obj.freeVertices, obj.freeVertices);
            
            if L == 1
                obj.p1Acoarse = p1Matrix;
            else
                obj.p1Smoother{L} = full(diag(p1Matrix)).^(-1);
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                obj.changedPatches{L} = obj.changedPatches{L}(obj.freeVertices);
                obj.patchwiseA = assemblePatchwise(obj.blf, obj.hoFes);
                obj.inclusionMatrix = spaceProlongation(obj.loFes, obj.hoFes);
                obj.inclusionMatrix = obj.inclusionMatrix(obj.freeVertices, freeVerticesHigher);

                % Global computation
                % obj.inclusionMatrix = interpolateFreeData(eye(length(obj.freeVertices)), obj.loFes, obj.hoFes);

                % obj.patchwiseP1Matrix{L} = assemblePatchwise(obj.blf, obj.loFes, obj.changedPatches{L});

                % DEBUG: exact on highest level
                % TMP = assemble(obj.blf, obj.hoFes);
                % obj.patchwiseA = TMP(getFreeDofs(obj.hoFes),getFreeDofs(obj.hoFes));
            end
            
            setupSystemMatrix@PcgSolver(obj, A);
        end
           
        function setupRhs(obj, b, varargin)
            setupRhs@PcgSolver(obj, b, varargin{:});
        end
        
        % preconditioner: inverse of diagonal on each level
        function Cx = preconditionAction(obj, residual)
            assert(~isempty(obj.nLevels), 'Data for multilevel iteration not given!')

            L = obj.nLevels;
            rho = cell(L, 1);

            if L == 1
                Cx = obj.A \ residual;
                return
            end

            % descending cascade
            rho{L} = hoGlobalSmoothing(obj, residual);

            % DEBUG: extra step size choice
            % [~, sUpd] = computeOptimalUpdate(obj.A, residual, rho{L});
            % rho{L} = sUpd;

            % Exact interpolation (NOT SYMMETRIC!)
            % residual = interpolateFreeData(residual, obj.hoFes, obj.loFes);

            residual = obj.inclusionMatrix * residual;

            % DEBUG: p1 smoothing only
            % rho{L} = p1LocalSmoothing(obj, L, residual);

            for k = L:-1:3
                residual = obj.intergridMatrix{k}'*residual;
                rho{k-1} = p1LocalSmoothing(obj, k-1, residual);
                % DEBUG: delete some corrections
                % if k > 5
                %     rho{k-1} = zeros(size(residual));
                % end
            end

            % exact solve on coarsest level
            residual = obj.intergridMatrix{2}'*residual;
            sigma = obj.p1Acoarse \ residual;
            sigma = obj.intergridMatrix{2}*sigma;

            % ascending cascade
            for k = 3:L
                sigma = sigma + rho{k-1};
                sigma = obj.intergridMatrix{k}*sigma;
            end

            % Correct interpolation (NOT SYMMETRIC!)
            % sigma = interpolateFreeData(sigma, obj.loFes, obj.hoFes);

            sigma = obj.inclusionMatrix' * sigma;
            sigma = sigma + rho{L};

            Cx = sigma;
        end
    end

    methods (Access=protected)
        % callback function to be executed before mesh refinement:
        % -> patches change for all vertices adjacent to refined edges and
        %    all new vertices
        function getChangedPatches(obj, mesh, bisecData)
            % Selects vertices of bisected edges and all newly created
            % vertices
            nCOld = mesh.nCoordinates;
            nCNew = mesh.nCoordinates + nnz(bisecData.bisectedEdges) + ...
                bisecData.nRefinedElements'*bisecData.nInnerNodes;
            bisectedEdgeNodes = unique(mesh.edges(:,bisecData.bisectedEdges));
            obj.changedPatches{obj.nLevels+1} = false(nCNew, 1);
            idx = [bisectedEdgeNodes; ((nCOld+1):nCNew)'];
            obj.changedPatches{obj.nLevels+1}(idx) = true;
        end

        function rho = p1LocalSmoothing(obj, k, res)
            idx = obj.changedPatches{k};
            rho = zeros(size(res));
            rho(idx,:) = obj.p1Smoother{k}(idx).*res(idx,:);
            % DEBUG: Local patchwise smoothing in P1
            % rho = obj.patchwiseP1Matrix{k} \ res;
        end
        
        function rho = hoGlobalSmoothing(obj, res)
            rho = obj.patchwiseA \ res;
        end
    end
end

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

 %% auxiliary functions
function interpolatedData = interpolateFreeData(data, fromFes, toFes)
    freeDofs = getFreeDofs(toFes);
    nComponents = size(data, 2);
    interpolatedData = zeros(numel(freeDofs), nComponents);
    feFunctionWrapper = FeFunction(fromFes);
    for k = 1:nComponents
        feFunctionWrapper.setFreeData(data(:,k));
        wholeData = nodalInterpolation(feFunctionWrapper, toFes);
        interpolatedData(:,k) = wholeData(freeDofs);
    end
end


%% DEBUG: COPIED FROM MG

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