% AdditiveSchwarzHighOrderPcg (subclass of PcgSolver) Solves linear equations
%   iteratively using the CG method with optimal multilevel additive
%   Schwartz preconditioner.

classdef AdditiveSchwarzHighOrderPcg < PcgSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end
    
    properties (Access=protected)
        blf
        hoFes
        loFes
        coarseLoFes
        coarseHoFes
        p1Acoarse
        smoother
        P
        intergridMatrix
        freeVertices
        freeVerticesOld
        freeDofs
        freeDofsOld
        changedPatches
        patchwiseA
    end

    properties (Access=private)
        listenerHandle
    end
    
    %% methods
    methods (Access=public)
        function obj = AdditiveSchwarzHighOrderPcg(fes, blf, P)
            arguments
                fes FeSpace
                blf BilinearForm
                P Prolongation
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
            obj.P = P;
            obj.blf = blf;

            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end
        
        function setupSystemMatrix(obj, A)
            obj.nLevels = obj.nLevels + 1;
            
            L = obj.nLevels;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.loFes);

            obj.freeDofsOld = obj.freeDofs;
            obj.freeDofs = getFreeDofs(obj.hoFes);

            ppMatrix = assemble(obj.blf, obj.hoFes);
            ppMatrix = ppMatrix(getFreeDofs(obj.hoFes), getFreeDofs(obj.hoFes));

            p1Matrix = assemble(obj.blf, obj.loFes);
            p1Matrix = p1Matrix(obj.freeVertices, obj.freeVertices);
            
            if L == 1
                obj.p1Acoarse = p1Matrix;
                fes = obj.hoFes;
                coarseMesh = clone(fes.mesh);
                obj.coarseLoFes = FeSpace(coarseMesh, LowestOrderH1Fe, 'dirichlet', fes.bnd.dirichlet);
                obj.coarseHoFes = FeSpace(coarseMesh, HigherOrderH1Fe(fes.finiteElement.order), 'dirichlet', fes.bnd.dirichlet);
            else
                obj.smoother{L} = full(diag(ppMatrix)).^(-1);
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeDofs, obj.freeDofsOld);
                obj.changedPatches{L} = obj.changedPatches{L}(obj.freeVertices);
                obj.patchwiseA{L} = assemblePatchwise(obj.blf, obj.hoFes, obj.changedPatches{L});
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
            for k = L:-1:2
                rho{k} = localSmoothing(obj, k, residual);
                residual = obj.intergridMatrix{k}'*residual;
            end

            % exact solve on coarsest level
            residual = interpolateFreeData(residual, obj.coarseHoFes, obj.coarseLoFes);
            sigma = obj.p1Acoarse \ residual;
            sigma = interpolateFreeData(sigma, obj.coarseLoFes, obj.coarseHoFes);
            
            % ascending cascade
            for k = 2:L
                sigma = obj.intergridMatrix{k}*sigma;
                sigma = sigma + rho{k};
            end
           
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
        
        function rho = localSmoothing(obj, k, res)
            rho = obj.patchwiseA{k} \ res;
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