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
        intergridMatrix
        freeDofs
        freeDofsOld
        changedPatches
        highestOrderIsOne
        patchwiseA
        Acoarse
        Afine
        coarsemesh
        coarseLoFes
        coarseHoFes
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
            
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', fes.bnd.dirichlet);
            obj.hoFes = fes;
            obj.blf = blf; % TODO: check coefficients of blf for compatibility with theoretical results!
            obj.Afine = cell(1,1);
            obj.Acoarse = assemble(obj.blf, obj.loFes);
            obj.Acoarse = obj.Acoarse(getFreeDofs(obj.loFes), getFreeDofs(obj.loFes));
            
            % set up FE spaces on coarsest level for projection
            if ~obj.highestOrderIsOne
                coarseMesh = clone(fes.mesh);
                obj.coarseLoFes = FeSpace(coarseMesh, LowestOrderH1Fe, 'dirichlet', fes.bnd.dirichlet);
                obj.coarseHoFes = FeSpace(coarseMesh, HigherOrderH1Fe(fes.finiteElement.order), 'dirichlet', fes.bnd.dirichlet);
            end
            
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end

        function setupSystemMatrix(obj, A, P)
            obj.nLevels = obj.nLevels + 1;
            
            L = obj.nLevels;
            freeVertices = getFreeDofs(obj.loFes);

            obj.freeDofsOld = obj.freeDofs;
            obj.freeDofs = getFreeDofs(obj.hoFes);
            
            if L >= 2
                obj.intergridMatrix{L} = P(obj.freeDofs, obj.freeDofsOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(freeVertices));

                if ~obj.highestOrderIsOne
                  obj.patchwiseA{L} = assemblePatchwise(obj.blf, obj.hoFes);
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
            obj.Afine{L} = obj.A;
            for k = L:-1:2
                residual{k-1} = obj.intergridMatrix{k}'*residual{k};
            end

            % coarse level in P1
            residual{1} = projectFromPto1(obj, residual{1});
            sigma = obj.Acoarse \ residual{1};
            algError2 = algError2 + scalarProduct(sigma, obj.Acoarse*sigma);

            
            % ascending cascade in Pp 
            sigma = projectFrom1toP(obj, sigma);
            sigma = obj.intergridMatrix{2}*sigma;
            
            %level 2 : ALL patches are smoothed IF high-order
            uptres = residual{2} - obj.Afine{2}*sigma;
            if obj.highestOrderIsOne
                rho = p1LocalSmoothing(obj, L, uptres); %FIX LATER: p=1 case
            else
                rho = obj.patchwiseA{2} \ uptres;
            end


            [aeUpd, sUpd] = computeOptimalUpdate(obj.Afine{2}, uptres, rho);
            sigma = sigma + sUpd;
            algError2 = algError2 + aeUpd;
            
            % under construction: ONLY modified patches are smoothed
            for k = 3:L
                sigma = obj.intergridMatrix{k}*sigma;
                uptres = residual{k} - obj.Afine{k}*sigma; %updated residual 
                rho = obj.patchwiseA{k} \ uptres;
                [aeUpd, sUpd] = computeOptimalUpdate(obj.Afine{k}, uptres, rho);
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
                y = interpolateCoarseFreeData(x, obj.coarseLoFes, obj.coarseHoFes);
            end
        end
        
        function y = projectFromPto1(obj, x)
            if obj.highestOrderIsOne
                y = x;
            else
                y = interpolateCoarseFreeData(x, obj.coarseHoFes, obj.coarseLoFes);
            end
        end
    end
end

%% auxiliary functions
function interpolatedData = interpolateCoarseFreeData(data, fromFes, toFes)
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
