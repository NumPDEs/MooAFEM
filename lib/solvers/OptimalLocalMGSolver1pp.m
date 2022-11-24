% OptimalLocalMGSolver (subclass of MGSolver) Solves linear equations
%   iteratively using high-order p-optimal multilevel local smoothing (additive
%   Schwarz/block Jacobi).
%   One iteration = Vcycle(0,1) (no pre-/one post-smoothing step).
%
% With the solver comes an 

classdef OptimalLocalMGSolver1pp < MGSolver
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
        Pinter
        intergridMatrix
        freeVertices
        freeVerticesOld
        freeDofs
        freeDofsOld
        changedPatches
        highestOrderIsOne
        patchwiseA
        coarsemesh
        freeVerticescoarse
        freeDofscoarse
    end
    
    %% event data
    properties (Access=private)
        listenerHandle
    end

    %% methods
    methods (Access=public)
        function obj = OptimalLocalMGSolver1pp(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end
            obj = obj@MGSolver();
            
            obj.highestOrderIsOne = (fes.finiteElement.order == 1);
            obj.nLevels = 0;
            
            mesh = fes.mesh;
            
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', ':');
            obj.hoFes = fes;
            obj.P = LoFeProlongation(obj.loFes);
            obj.blf = blf; % TODO: check coefficients of blf for compatibility with theoretical results!
            
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end

        function setupSystemMatrix(obj, A, Pinter)
            obj.nLevels = obj.nLevels + 1;
            
            L = obj.nLevels;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.loFes);

            if L == 1
             obj.freeVerticescoarse = getFreeDofs(obj.loFes);
             obj.freeDofscoarse = getFreeDofs(obj.hoFes);
            end

            obj.freeDofsOld = obj.freeDofs;
            obj.freeDofs = getFreeDofs(obj.hoFes);

            
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
                obj.Pinter{L} = Pinter(obj.freeDofs, obj.freeDofsOld);

                if ~obj.highestOrderIsOne
                  obj.patchwiseA{L} = assemblePatchwise(obj.blf, obj.hoFes, L, obj.changedPatches{L});
                end
            end
            
            setupSystemMatrix@MGSolver(obj, A);
        end

        function setupRhs(obj, b, x0)
            setupRhs@MGSolver(obj, b, x0);
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
            
            % descending cascade in Pp: no smoothing
            residual{L} = res;
            ppmatrix{L} = obj.A;
            for k = L:-1:2
                ppmatrix{k-1} = obj.Pinter{k}'*ppmatrix{k}*obj.Pinter{k};
                residual{k-1} = obj.Pinter{k}'*residual{k};
            end

            % coarse level in P1
            residual{1} = projectFromPto1(obj, residual{1});
            sigma = obj.p1Matrix{1} \ residual{1};
            algError2 = algError2 + scalarProduct(sigma, obj.p1Matrix{1}*sigma);

            
            % ascending cascade in Pp 
            sigma = projectFrom1toP(obj, sigma);
            sigma = obj.Pinter{2}*sigma;
            
            %level 2 : ALL patches are smoothed IF high-order
            uptres = residual{2} - ppmatrix{2}*sigma;
            if obj.highestOrderIsOne
                rho = p1LocalSmoothing(obj, L, uptres); %FIX LATER: p=1 case
            else 
                %getPatchDofs(obj.patchwiseA{2}, 1)
                rho = obj.patchwiseA{2} \ uptres;
            end


            [aeUpd, sUpd] = computeOptimalUpdate(ppmatrix{2}, uptres, rho);
            sigma = sigma + sUpd;
            algError2 = algError2 + aeUpd;
            
            % under construction: ONLY modified patches are smoothed
            for k = 3:L
                sigma = obj.Pinter{k}*sigma;
                uptres = residual{k} - ppmatrix{k}*sigma; %updated residual 
                rho = obj.patchwiseA{k} \ uptres; %p1LocalSmoothing(obj, k, uptres);
                [aeUpd, sUpd] = computeOptimalUpdate(ppmatrix{k}, uptres, rho);
                sigma = sigma + sUpd;
                algError2 = algError2 + aeUpd;
            end
            
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
                %y = interpolateFreeData(x, obj.loFes, obj.hoFes, obj.freeVerticesOld);
                y = interpolateFreeDatacoarse(x, obj.loFes, obj.hoFes, obj.freeDofscoarse);
            end
        end
        
        function y = projectFromPto1(obj, x)
            if obj.highestOrderIsOne
                y = x;
            else
                %y = interpolateFreeData(x, obj.hoFes, obj.loFes);
                y = interpolateFreeDatacoarse(x, obj.hoFes, obj.loFes, obj.freeVerticescoarse);
            end
        end
    end
end

%% auxiliary functions

function interpolatedData = interpolateFreeDatacoarse(data, fromFes, toFes, freeVer)
    coarsemesh = Mesh.loadFromGeometry('Lshape');  %should be general depending which problem is used
    freeDofs = freeVer;
    nComponents = size(data, 2);
    interpolatedData = zeros(numel(freeDofs), nComponents);
    

        fromFescoarse = FeSpace(coarsemesh, HigherOrderH1Fe(fromFes.finiteElement.order), 'dirichlet', ':');
        toFescoarse =  FeSpace(coarsemesh, HigherOrderH1Fe(toFes.finiteElement.order), 'dirichlet', ':');

     feFunctionWrapper = FeFunction(fromFescoarse);
    for k = 1:nComponents
        feFunctionWrapper.setFreeData(data(:,k));
        wholeData = nodalInterpolation(feFunctionWrapper, toFescoarse);
        interpolatedData(:,k) = wholeData(freeDofs);
    end
end

function interpolatedData = interpolateFreeData(data, fromFes, toFes)
    freeDofs = getFreeDofs(toFes);
    getFreeDofs(obj.loFes);
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
