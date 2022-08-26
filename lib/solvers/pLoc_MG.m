% pLoc_MG (subclass of MGSolver) Solves linear equations
%   iteratively using high order p optimal multilevel
%   local (additive Schwarz/block Jacobi) smoothing
%   one iteration = Vcycle(0,1) : only 1 post-smoothing

classdef pLoc_MG < MGSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end

    properties (Access=protected)
        p1Matrix
        p1Smoother
        highestOrderIsOne
        blf
        hoFes
        loFes
        intergridMatrix
        patchwiseChol
        patch
        changedPatches
        freeVertices
        freeVerticesOld
        P
    end
    
    %% event data
    properties (Access=private)
        listenerHandle
    end

    %% methods
    methods (Access=public)
        function obj = pLoc_MG(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end
            obj = obj@MGSolver();
            
            obj.hoFes = fes;
            obj.nLevels = 0;
            obj.patchwiseChol = {};
            obj.intergridMatrix = cell(0);
            obj.changedPatches = cell(1);
            obj.p1Matrix = {};
            obj.p1Smoother = {};
            obj.highestOrderIsOne = (obj.hoFes.finiteElement.order == 1);
            
            mesh = obj.hoFes.mesh;
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', ':');
            obj.P = LoFeProlongation(obj.loFes);
            obj.blf = blf;
            % TODO: check coefficients of blf for compatibility with theoretical results!
            
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end

        function setupLinearSystem(obj, A, b, x0)
            obj.nLevels = obj.nLevels + 1;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.loFes);
            obj.p1Matrix{obj.nLevels} = assemble(obj.blf, obj.loFes);
            obj.p1Matrix{obj.nLevels} = obj.p1Matrix{obj.nLevels}(obj.freeVertices, obj.freeVertices);
            obj.p1Smoother{obj.nLevels} = diag(obj.p1Matrix{obj.nLevels}).^(-1);
            
            if obj.nLevels >= 2
                obj.intergridMatrix{obj.nLevels} = ...
                    obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                
                localVer = obj.changedPatches{obj.nLevels}; %numbering on all vertices
                [~,localVer] = ismember(localVer,obj.freeVertices);
                obj.changedPatches{obj.nLevels} = localVer(localVer>0); %numbering on all inner vertices
            end

            freeDofs = getFreeDofs(obj.hoFes);
            global2freeDofs = zeros(getDofs(obj.hoFes).nDofs, 1);
            global2freeDofs(freeDofs) = 1:numel(freeDofs);
            
            if ~obj.highestOrderIsOne
                % store cholesky decomposition of local matrices associated
                % to patches (only upper triangle)
                obj.patch = assemblePatchDofs(obj.hoFes);
                obj.patch = cellfun(@(p) global2freeDofs(p), obj.patch, 'UniformOutput', false);
                obj.patchwiseChol = cellfun(@(p) full(chol(A(p,p))), obj.patch, 'UniformOutput', false);
            end

            setupLinearSystem@MGSolver(obj, A, b, x0);
        end

        % Geometric MultiGrid
        function [Cx, algError2] = Vcycle(obj, res)
            assert(isequal(size(res, 1), size(obj.A, 1)), ...
                'Setup for multilevel iteration not complete!')
            
            % if there is only one coarse level: exact solve
            if obj.nLevels == 1
                Cx = obj.A \ res;
                algError2 = Cx'*obj.A*Cx;
                return
            end

            L = obj.nLevels;
            algError2 = zeros(1, size(res, 2)); % built-in estimator of the algebraic error
            
            % descending cascade in P1: no smoothing
            residual{L} = projectFromPto1(obj, res);
            for k = L:-1:2
                residual{k-1}  = obj.intergridMatrix{k}'*residual{k};
            end

            % exact solve on coarsest level to compute accumulative lifting
            % of the residual (sigma)
            sigma = obj.p1Matrix{1} \ residual{1};
            algError2 = algError2 + sigma'*obj.p1Matrix{1}*sigma;

            % ascending cascade in P1 AND local smoothing
            for k = 2:(L-1)
                sigma = obj.intergridMatrix{k}*sigma;
                uptres = residual{k} - obj.p1Matrix{k}*sigma; %updated residual 
                rho = p1LocalSmoothing(obj, k, uptres);
                [aeUpd, sUpd] = computeEstimatorUpdate(obj.p1Matrix{k}, uptres, rho);
                algError2 = algError2 + aeUpd; sigma = sigma + sUpd;
            end
            
            % smooting on finest level dependent on p
            sigma = obj.intergridMatrix{L}*sigma;
            sigma = projectFrom1toP(obj, sigma);
            if obj.highestOrderIsOne
                uptres = residual{L} - obj.A*sigma;
                rho = p1LocalSmoothing(obj, L, uptres);
            else
                % finest level high-order patch-problems ONLY for p > 1
                uptres = res - obj.A*sigma;
                rho = hoGlobalSmoothing(obj, uptres);
            end
            [aeUpd, sUpd] = computeEstimatorUpdate(obj.A, uptres, rho);
            algError2 = algError2 + aeUpd; sigma = sigma + sUpd;

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
            obj.changedPatches{obj.nLevels+1} = ...
                [unique(bisectedEdgeNodes); ((nCOld+1):nCNew)'];
        end
        
        function rho = p1LocalSmoothing(obj, k, res)
            vershift = obj.changedPatches{k};
            rho = zeros(size(res)); 
            rho(vershift,:) = obj.p1Smoother{k}(vershift).*res(vershift,:);
        end
        
        function rho = hoGlobalSmoothing(obj, res)
            rho = zeros(size(res));
            lower = struct('UT', true, 'TRANSA', true); 
            upper = struct('UT', true);
            for ver = 1:obj.hoFes.mesh.nCoordinates
                U = obj.patchwiseChol{ver};
                dofshift = obj.patch{ver};

                % solve local patch problems (update is additive)
                rhoa = linsolve(U, linsolve(U, res(dofshift,:), lower), upper);
                rho(dofshift,:) = rho(dofshift,:) + rhoa;
            end
        end
        
        function y = projectFrom1toP(obj, x)
            if obj.highestOrderIsOne
                y = x;
            else
                y = interpolateFreeData(x, obj.loFes, obj.hoFes);
            end
        end
        
        function y = projectFromPto1(obj, x)
            if obj.highestOrderIsOne
                y = x;
            else
                y = interpolateFreeData(x, obj.hoFes, obj.loFes);
            end
        end
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

% error correction with optimal stepsize 
function [etaUpdate, sigmaUpdate] = computeEstimatorUpdate(A, res, rho)
    lambda = (res'*rho)/(rho'*A*rho);
    sigmaUpdate = lambda*rho;
    etaUpdate = (rho'*A*rho)*(lambda^2);

    if (lambda > 3)
       warning('MG step-sizes no longer bound by d+1. Optimality of step size cannot be guaranteed!')
    end
end