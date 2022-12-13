% OptimalLocalMGSolver (subclass of MGSolver) Solves linear equations
%   iteratively using high-order p-optimal multilevel local smoothing (additive
%   Schwarz/block Jacobi).
%   One iteration = Vcycle(0,1) (no pre-/one post-smoothing step).
%
% With the solver comes an 

classdef LocalMgHighOrderVcycle < MGSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end

    properties (Access=protected)
        blf
        hoFes
        loFes
        P
        intergridMatrix
        freeDofs
        freeDofsOld
        changedPatches
        patchwiseA
        Acoarse
        Afine
        coarseLoFes
        coarseHoFes
    end
    
    %% event data
    properties (Access=private)
        listenerHandle
    end

    %% methods
    methods (Access=public)
        function obj = LocalMgHighOrderVcycle(fes, blf, P)
            arguments
                fes FeSpace
                blf BilinearForm
                P Prolongation
            end
            obj = obj@MGSolver();
            
            assert(fes.finiteElement.order > 1, ...
                'LocalMgLowOrderVcycle only works for higher order finite elements.')
            assert(isempty(blf.b), ...
                'Multigrid solvers only tested for symmetric problems.')
            
            obj.nLevels = 0;
            
            mesh = fes.mesh;
            obj.P = P;
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', fes.bnd.dirichlet);
            obj.hoFes = fes;
            obj.blf = blf;
            obj.Afine = cell(1,1);
            obj.Acoarse = assemble(obj.blf, obj.loFes);
            obj.Acoarse = obj.Acoarse(getFreeDofs(obj.loFes), getFreeDofs(obj.loFes));
            
            % set up FE spaces on coarsest level for projection
            coarseMesh = clone(fes.mesh);
            obj.coarseLoFes = FeSpace(coarseMesh, LowestOrderH1Fe, 'dirichlet', fes.bnd.dirichlet);
            obj.coarseHoFes = FeSpace(coarseMesh, HigherOrderH1Fe(fes.finiteElement.order), 'dirichlet', fes.bnd.dirichlet);
            
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end

        function setupSystemMatrix(obj, A)
            obj.nLevels = obj.nLevels + 1;
            
            L = obj.nLevels;
            freeVertices = getFreeDofs(obj.loFes);

            obj.freeDofsOld = obj.freeDofs;
            obj.freeDofs = getFreeDofs(obj.hoFes);
            
            if L >= 2
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeDofs, obj.freeDofsOld);
                obj.changedPatches{L} = freeVertices(obj.changedPatches{L}(freeVertices));

                if L == 2
                    obj.patchwiseA{L} = assemblePatchwise(obj.blf, obj.hoFes, ':');
                else
                    obj.patchwiseA{L} = assemblePatchwise(obj.blf, obj.hoFes, obj.changedPatches{L});
                end
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
            
            % descending cascade in p: no smoothing
            residual{L} = res;
            obj.Afine{L} = obj.A;
            for k = L:-1:2
                residual{k-1} = obj.intergridMatrix{k}'*residual{k};
            end

            % coarse level in P1
            residual{1} = interpolateCoarseFreeData(residual{1}, obj.coarseHoFes, obj.coarseLoFes);
            sigma = obj.Acoarse \ residual{1};
            algError2 = algError2 + scalarProduct(sigma, obj.Acoarse*sigma);
            
            % ascending cascade in Pp 
            sigma = interpolateCoarseFreeData(sigma, obj.coarseLoFes, obj.coarseHoFes);
            
            for k = 2:L
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
