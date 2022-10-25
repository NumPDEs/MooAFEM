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
        
        function setupLinearSystem(obj, A, b, x0)
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

            % descending cascade
            for k = L-1:-1:1
                y{k+1} = obj.D{k} .* x;
                x = obj.intergridMatrix{k}'*x;
            end
            
            % exact solve on coarsest level
            y{1} = obj.Acoarse \ x;
            
            % ascending cascade
            for k = 1:obj.nLevels-1
                y{k+1} = y{k+1} + obj.intergridMatrix{k}*y{k};
            end
            
            Cx = y{end};
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