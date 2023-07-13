% AdditiveSchwarzLowOrderCascade (subclass of Preconditioner) optimal multilevel
%   additive Schwarz preconditioner for arbitrary order finite elements: smooth
%   with diagonal of P1 matrix on changed patches on every level (local P1
%   smoothing) and with patchwise higher order matrix on finest level (global
%   higher order smoothing).
%
% See also: Preconditioner, PcgSolver

classdef AdditiveSchwarzLowOrderCascade < Preconditioner
    %% properties
    properties (GetAccess=public,SetAccess=protected)
    end
    
    properties (Access=protected)
        blf
        hoFes
        loFes
        P
        p1Acoarse
        hoAcoarse
        p1Smoother
        nLevels
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
        function obj = AdditiveSchwarzLowOrderCascade(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end
            
            assert(fes.finiteElement.order > 1, ...
                'Higher order additive Schwarz preconditioner only works for higher order finite elements.')
            assert(isempty(blf.b), ...
                'Additive Schwarz preconditioner only tested for symmetric problems.')

            obj.nLevels = 0;

            mesh = fes.mesh;
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe, 'dirichlet', fes.bnd.dirichlet);
            obj.hoFes = fes;
            obj.P = Prolongation.chooseFor(obj.loFes);
            obj.blf = blf;

            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end
        
        function setup(obj, A)
            obj.nLevels = obj.nLevels + 1;
            
            L = obj.nLevels;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.loFes);

            p1Matrix = assemble(obj.blf, obj.loFes);
            p1Matrix = p1Matrix(obj.freeVertices, obj.freeVertices);
            
            if L == 1
                obj.p1Acoarse = p1Matrix;
                obj.hoAcoarse = A;
            else
                obj.p1Smoother{L} = full(diag(p1Matrix)).^(-1);
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                obj.changedPatches{L} = obj.changedPatches{L}(obj.freeVertices);
                obj.patchwiseA = assemblePatchwise(obj.blf, obj.hoFes);
                obj.inclusionMatrix = SpaceProlongation(obj.loFes, obj.hoFes) ...
                    .matrix(getFreeDofs(obj.loFes), getFreeDofs(obj.hoFes));
            end
        end
        
        % preconditioner: inverse of diagonal on each level
        function Cx = apply(obj, x)
            assert(~isempty(obj.nLevels), 'Data for multilevel iteration not given!')

            L = obj.nLevels;
            rho = cell(L, 1);

            if L == 1
                Cx = obj.hoAcoarse \ x;
                return
            end

            % descending cascade
            rho{L} = hoGlobalSmoothing(obj, x);
            x = obj.inclusionMatrix * x;

            for k = L:-1:3
                x = obj.intergridMatrix{k}' * x;
                rho{k-1} = p1LocalSmoothing(obj, k-1, x);
            end

            % exact solve on coarsest level
            x = obj.intergridMatrix{2}' * x;
            Cx = obj.p1Acoarse \ x;
            Cx = obj.intergridMatrix{2}*Cx;

            % ascending cascade
            for k = 3:L
                Cx = Cx + rho{k-1};
                Cx = obj.intergridMatrix{k} * Cx;
            end

            Cx = obj.inclusionMatrix' * Cx;
            Cx = Cx + rho{L};
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
                bisecData.nRefinedElements' * bisecData.nInnerNodes;
            bisectedEdgeNodes = unique(mesh.edges(:,bisecData.bisectedEdges));
            obj.changedPatches{obj.nLevels+1} = false(nCNew, 1);
            idx = [bisectedEdgeNodes; ((nCOld+1):nCNew)'];
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
    end
end