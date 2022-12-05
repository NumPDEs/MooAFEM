% AdditiveSchwartzPcg (subclass of PcgSolver) Solves linear equations
%   iteratively using the CG method with optimal multilevel additive
%   Schwartz preconditioner.

classdef AdditiveSchwartzPcg < PcgSolver
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
        freeVertices
        freeVerticesOld
        changedPatches
        highestOrderIsOne
        patchwiseA
    end

    properties (Access=private)
        listenerHandle
    end
    
    %% methods
    methods (Access=public)
        function obj = AdditiveSchwartzPcg(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end
            obj = obj@PcgSolver();
            
            obj.highestOrderIsOne = (fes.finiteElement.order == 1);
            obj.nLevels = 0;
            
            mesh = fes.mesh;
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', ':');
            obj.hoFes = fes;
            obj.P = LoFeProlongation(obj.loFes);
            obj.blf = blf; % TODO: check coefficients of blf for compatibility with theoretical results!
            
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
            obj.intergridMatrix = [];
        end
        
        function setupSystemMatrix(obj, A)
            obj.nLevels = obj.nLevels + 1;
            
            L = obj.nLevels;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.loFes);
            
            if obj.highestOrderIsOne
                obj.p1Matrix{L} = A;
            else
                obj.p1Matrix{L} = assemble(obj.blf, obj.loFes);
                obj.p1Matrix{L} = obj.p1Matrix{L}(obj.freeVertices, obj.freeVertices);
            end
            obj.p1Smoother{L} = full(diag(obj.p1Matrix{L})).^(-1);
            
            if L >= 2
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(obj.freeVertices));

                if ~obj.highestOrderIsOne
                    obj.patchwiseA = assemblePatchwise(obj.blf, obj.hoFes);
                end
            end
            
            setupSystemMatrix@PcgSolver(obj, A);
        end
           
        function setupRhs(obj, b, x0)
            setupRhs@PcgSolver(obj, b, x0);
        end
        
        % preconditioner: inverse of diagonal on each level
        function Cx = preconditionAction(obj, res)
            assert(~isempty(obj.nLevels), 'Data for multilevel iteration not given!')
            y = cell(obj.nLevels, 1);

            if obj.nLevels == 1
                Cx = obj.A \ res;
                return
            end
            
            L = obj.nLevels;
            
            % descending cascade
            residual{L} = projectFromPto1(obj, res);
            for k = L:-1:2
                residual{k-1} = obj.intergridMatrix{k}'*(obj.p1Smoother{k} .* residual{k});
            end
            
            % exact solve on coarsest level
            sigma = obj.p1Matrix{1} \ residual{1};
            
            % ascending cascade
            for k = 2:obj.nLevels-1
                sigma = obj.intergridMatrix{k}*sigma;
                uptres = residual{k} - obj.p1Matrix{k}*sigma; %updated residual 
                rho = p1LocalSmoothing(obj, k, uptres);
                sigma = sigma + rho;
            end
            
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
            sigma = sigma + rho;

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
function [etaUpdate, sigmaUpdate] = computeOptimalUpdate(A, res, rho)
    rhoArho = scalarProduct(rho, A*rho);
    lambda = scalarProduct(res, rho) ./ rhoArho;
    sigmaUpdate = lambda.*rho;
    etaUpdate = rhoArho.*(lambda.^2);

    if any(lambda > 3)
       warning('MG step-sizes no longer bound by d+1. Optimality of step size cannot be guaranteed!')
    end
end

function a = scalarProduct(x, y)
    a = sum(x.*y, 1);
end