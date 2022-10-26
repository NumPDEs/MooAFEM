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
        D
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
        
        function setupLinearSystem(obj, A, b, x0)
            obj.nLevels = obj.nLevels + 1;

            L = obj.nLevels;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.loFes);
            
            if obj.highestOrderIsOne
                obj.p1Matrix{L} = A;
                obj.D{L} = 1./full(diag(obj.p1Matrix{L}));
            else
                obj.p1Matrix{L} = assemble(obj.blf, obj.loFes);
                obj.p1Matrix{L} = obj.p1Matrix{L}(obj.freeVertices, obj.freeVertices);
                obj.D{L} = 1./full(diag(obj.p1Matrix{L}));
            end

            if L >= 2
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(obj.freeVertices));
                if ~obj.highestOrderIsOne
                    obj.patchwiseA = assemblePatchwise(obj.blf, obj.hoFes);
                end
            end
            
            setupLinearSystem@PcgSolver(obj, A, b, x0);
        end
        
        % preconditioner: inverse of diagonal on each level
        function Cx = preconditionAction(obj, x)
            assert(~isempty(obj.nLevels), 'Data for multilevel iteration not given!')
            y = cell(obj.nLevels, 1);

            if obj.nLevels == 1
                Cx = obj.A \ x;
                return
            end
            
            L = obj.nLevels;
            y = cell(L, 1);
            
            % descending cascade in P1: no smoothing
            residual{L} = projectFromPto1(obj, x);
            for k = L-1:-1:1
                residual{k} = obj.D{k+1} .* x;
                x  = obj.intergridMatrix{k+1}'*x;
            end

            % exact solve on coarsest level to compute accumulative lifting
            % of the residual (sigma)
            sigma = obj.p1Matrix{1} \ x;

            % ascending cascade in P1 AND local smoothing
            for k = 2:(L-1)
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