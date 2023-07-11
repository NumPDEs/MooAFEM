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

            p1Matrix = assemble(obj.blf, obj.loFes);
            p1Matrix = p1Matrix(obj.freeVertices, obj.freeVertices);
            
            if L == 1
                obj.p1Acoarse = p1Matrix;
            else
                obj.p1Smoother{L} = full(diag(p1Matrix)).^(-1);
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                obj.changedPatches{L} = obj.changedPatches{L}(obj.freeVertices);
                obj.patchwiseA = assemblePatchwise(obj.blf, obj.hoFes);
                obj.inclusionMatrix = SpaceProlongation(obj.loFes, obj.hoFes) ...
                    .matrix(getFreeDofs(obj.loFes), getFreeDofs(obj.hoFes));
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

            residual = obj.inclusionMatrix * residual;

            for k = L:-1:3
                residual = obj.intergridMatrix{k}'*residual;
                rho{k-1} = p1LocalSmoothing(obj, k-1, residual);
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