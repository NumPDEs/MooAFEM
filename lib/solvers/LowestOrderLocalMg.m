% LowestOrderLocalMg (subclass of MGSolver) Solves linear equations
%   stemming from lowest order AFEM computations iteratively using a
%   multigrid solver where one iteration is characterized by Vcycle(0,1)
%   (no pre-/one post-smoothing step).
%   On every level, local Jacobi smoothing is employed on changed patches
%   and the smoothed update is weighted by an optimal step size.
%
%   See also [Innerberger, Miraçi, Praetorius, Streitberger; Optimal
%   computational costs of AFEM with optimal local hp-robust multigrid
%   solver].

classdef LowestOrderLocalMg < MGSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end

    properties (Access=protected)
        fes
        blf
        matrix
        smoother
        P
        intergridMatrix
        freeVertices
        freeVerticesOld
        changedPatches
    end
    
    %% event data
    properties (Access=private)
        listenerHandle
    end

    %% methods
    methods (Access=public)
        function obj = LowestOrderLocalMg(fes, blf, P)
            arguments
                fes FeSpace
                blf BilinearForm
                P Prolongation
            end
            obj = obj@MGSolver();
            
            assert(fes.finiteElement.order == 1, ...
                'LowestOrderLocalMg only works for lowest order finite elements.')
            assert(isempty(blf.b), ...
                'Multigrid solvers only tested for symmetric problems.')
            
            obj.nLevels = 0;
            
            mesh = fes.mesh;
            obj.fes = fes;
            obj.P = P;
            obj.blf = blf;
            
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end

        function setupSystemMatrix(obj, A)
            obj.nLevels = obj.nLevels + 1;
            
            L = obj.nLevels;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.fes);
            
            obj.matrix{L} = A;
            obj.smoother{L} = full(diag(obj.matrix{L})).^(-1);
            
            if L >= 2
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(obj.freeVertices));
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
            residual{L} = res;
            for k = L:-1:2
                residual{k-1}  = obj.intergridMatrix{k}'*residual{k};
            end

            % exact solve on coarsest level to compute accumulative lifting
            % of the residual (sigma)
            sigma = obj.matrix{1} \ residual{1};
            algError2 = algError2 + scalarProduct(sigma, obj.matrix{1}*sigma);

            % ascending cascade in P1 AND local smoothing
            for k = 2:L
                sigma = obj.intergridMatrix{k}*sigma;
                uptres = residual{k} - obj.matrix{k}*sigma; %updated residual 
                rho = localSmoothing(obj, k, uptres);
                [aeUpd, sUpd] = computeOptimalUpdate(obj.matrix{k}, uptres, rho);
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
        
        function rho = localSmoothing(obj, k, res)
            idx = obj.changedPatches{k};
            rho = zeros(size(res));
            rho(idx,:) = obj.smoother{k}(idx).*res(idx,:);
        end
    end
end

%% auxiliary functions
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
