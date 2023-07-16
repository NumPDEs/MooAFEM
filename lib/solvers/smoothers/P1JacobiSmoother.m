% P1JacobiSmoother (subclass of MultilevelSmoother) multilevel Jacobi smoother
%   for lowest order finite elements: smooth with diagonal of stiffness matrix
%   on changed patches on every level.
%
% See also: MultilevelSmoother, OptimalVcycleMultigridSolver


classdef P1JacobiSmoother < MultilevelSmoother
    properties (Access=protected)
        fes
        blf
        inverseDiagonal
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
        function obj = P1JacobiSmoother(fes, blf, P)
            arguments
                fes FeSpace
                blf BilinearForm
                P Prolongation
            end
            obj = obj@MultilevelSmoother();
            
            assert(fes.finiteElement.order == 1, ...
                'P1JacobiSmoother only works for lowest order finite elements.')
            
            obj.fes = fes;
            obj.P = P;
            obj.blf = blf;
            
            mesh = fes.mesh;
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end

        function setup(obj, A)
            setup@MultilevelSmoother(obj, A);
            
            L = obj.nLevels;
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.fes);
            
            obj.inverseDiagonal{L} = full(diag(A)).^(-1);
            
            if L >= 2
                obj.intergridMatrix{L} = obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(obj.freeVertices));
                obj.inverseDiagonal{L} = obj.inverseDiagonal{L}(obj.changedPatches{L},:);
            end
        end

        function Cx = smooth(obj, x, k)
            idx = obj.changedPatches{k};
            Cx = zeros(size(x));
            Cx(idx,:) = obj.inverseDiagonal{k} .* x(idx,:);
        end

        function Px = prolongate(obj, x, k)
            Px = obj.intergridMatrix{k} * x;
        end

        function Px = restrict(obj, x, k)
            Px = obj.intergridMatrix{k}' * x;
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