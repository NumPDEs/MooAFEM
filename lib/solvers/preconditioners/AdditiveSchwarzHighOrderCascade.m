% AdditiveSchwarzHighOrderCascade (subclass of Preconditioner) optimal
%   multilevel additive Schwarz preconditioner for arbitrary order finite
%   elements: smooth with patchwise higher order matrix on changed patches on
%   all levels (local higher order smoothing) and globally with patchwise higher
%   order matrix on finest level (global higher order smoothing).
%
% See also: Preconditioner, PcgSolver

% TODO: this does not work as intended! Test and fix it!
classdef AdditiveSchwarzHighOrderCascade < Preconditioner
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end
    
    properties (Access=protected)
        blf
        hoFes
        loFes
        p1Acoarse
        hoAcoarse
        P
        inclusionMatrix
        intergridMatrix
        freeVertices
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
        function obj = AdditiveSchwarzHighOrderCascade(fes, blf, P)
            arguments
                fes FeSpace
                blf BilinearForm
                P Prolongation
            end
            
            assert(fes.finiteElement.order > 1, ...
                'Higher order additive Schwarz preconditioner only works for higher order finite elements.')
            assert(isempty(blf.b), ...
                'Additive Schwarz preconditioner only tested for symmetric problems.')

            obj.nLevels = 0;

            mesh = fes.mesh;
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe, 'dirichlet', fes.bnd.dirichlet);
            obj.hoFes = fes;
            obj.P = P;
            obj.blf = blf;

            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end
        
        function setup(obj, A)
            obj.nLevels = obj.nLevels + 1;
            
            L = obj.nLevels;
            obj.freeVertices = getFreeDofs(obj.loFes);

            obj.freeDofsOld = obj.freeDofs;
            obj.freeDofs = getFreeDofs(obj.hoFes);

            p1Matrix = assemble(obj.blf, obj.loFes);
            p1Matrix = p1Matrix(obj.freeVertices, obj.freeVertices);
            
            if L == 1
                obj.p1Acoarse = p1Matrix;
                obj.hoAcoarse = A;
                obj.inclusionMatrix = SpaceProlongation(obj.loFes, obj.hoFes) ...
                    .matrix(getFreeDofs(obj.loFes), getFreeDofs(obj.hoFes));
            else
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeDofs, obj.freeDofsOld);
                obj.changedPatches = obj.changedPatches(obj.freeVertices);
                if L == 2
                    % since order changes to 2nd mesh, we have to do global smoothing
                    obj.patchwiseA{L} = assemblePatchwise(obj.blf, obj.hoFes, ':');
                else
                    obj.patchwiseA{L} = assemblePatchwise(obj.blf, obj.hoFes, obj.changedPatches);
                end
            end
        end
        
        % preconditioner: patchwise higher order solve on changed patches
        function Cx = apply(obj, x)
            assert(~isempty(obj.nLevels), 'Data for multilevel iteration not given!')

            L = obj.nLevels;
            rho = cell(L, 1);

            if L == 1
                Cx = obj.hoAcoarse \ x;
                return
            end

            % descending cascade
            for k = L:-1:2
                rho{k} = localSmoothing(obj, k, x);
                x = obj.intergridMatrix{k}' * x;
            end

            % exact solve on coarsest level
            x = obj.inclusionMatrix * x;
            Cx = obj.p1Acoarse \ x;
            Cx = obj.inclusionMatrix' * Cx;
            
            % ascending cascade
            for k = 2:L
                Cx = obj.intergridMatrix{k} * Cx;
                Cx = Cx + rho{k};
            end
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
            obj.changedPatches = false(nCNew, 1);
            idx = [bisectedEdgeNodes; ((nCOld+1):nCNew)'];
            obj.changedPatches(idx) = true;
        end
        
        function rho = localSmoothing(obj, k, res)
            rho = obj.patchwiseA{k} \ res;
        end
    end
end