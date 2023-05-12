% LowestOrderAdditiveSchwartzPcg (subclass of PcgSolver) Solves linear
%   equations iteratively using the CG method with optimal multilevel
%   additive Schwartz preconditioner for lowest order finite elements.

classdef LowestOrderAdditiveSchwarzPcg < PcgSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end
    
    properties (Access=protected)
        blf
        fes
        matrix
        smoother
        P
        intergridMatrix
        freeDofs
        freeDofsOld
        changedPatches
        highestOrderIsOne
    end

    properties (Access=private)
        listenerHandle
    end
    
    %% methods
    methods (Access=public)
        function obj = LowestOrderAdditiveSchwarzPcg(fes, blf, P)
            arguments
                fes FeSpace
                blf BilinearForm
                P Prolongation
            end
            obj = obj@PcgSolver();
            
            assert(fes.finiteElement.order == 1, ...
                'LowestOrderAdditiveSchwartzPcg only works for lowest order finite elements.')
            assert(isempty(blf.b), ...
                'Additive Schwarz PCG solvers only tested for symmetric problems.')
            
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
            obj.freeDofsOld = obj.freeDofs;
            obj.freeDofs = getFreeDofs(obj.fes);
            
            obj.matrix{L} = A;
            obj.smoother{L} = full(diag(obj.matrix{L})).^(-1);
            
            if L >= 2
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeDofs, obj.freeDofsOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(obj.freeDofs));
            end
            
            setupSystemMatrix@PcgSolver(obj, A);
        end
           
        function setupRhs(obj, b, varargin)
            setupRhs@PcgSolver(obj, b, varargin{:});
        end
        
        % preconditioner: inverse of diagonal on each level
        function Cx = preconditionAction(obj, res)
            assert(~isempty(obj.nLevels), 'Data for multilevel iteration not given!')

            L = obj.nLevels;
            rho = cell(L, 1);

            if L == 1
                Cx = obj.A \ res;
                return
            end

            % descending cascade
            rho{L} = localSmoothing(obj, L, res);
            residual = res;
            for k = L-1:-1:2
                residual = obj.intergridMatrix{k+1}'*residual;
                rho{k} = localSmoothing(obj, k, residual);
            end
            
            % exact solve on coarsest level
            residual = obj.intergridMatrix{2}'*residual;
            sigma = obj.matrix{1} \ residual;
            
            % ascending cascade
            for k = 2:L-1
                sigma = obj.intergridMatrix{k}*sigma;
                sigma = sigma + rho{k};
            end
            sigma = obj.intergridMatrix{L}*sigma + rho{L};
           
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