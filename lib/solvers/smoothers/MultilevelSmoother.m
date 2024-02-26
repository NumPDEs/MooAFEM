% MultilevelSmoother (abstract handle class) Interface for smoothers to use in a
%   multigrid solver.
%
% Abstract methods to implement:
%
% smoother.setup(A) hook to set up the smoother for the matrix A.
%
% nonInvertedSmootherMatrix = smoother.setup(A) set up the smoother and also
%   return the non-inverted matrix, based on which the smoother is computed.
%
% Cx = smoother.smooth(x, k) applies the smoother to the vector x on level k.
%
% Px = smoother.prolongate(x, k) prolongate vector x from level k to level k+1.
%
% Px = smoother.restrict(x, k) restrict vector x from level k+1 to level k.
%
% See also: OptimalVcycleMultigridSolver


classdef MultilevelSmoother < handle
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        nLevels
    end

    properties (Access=protected)
        fes
        blf
        changedPatches
    end
    
    properties (Access=private)
        listenerHandle
    end

    %% methods
    methods (Access=public)
        function obj = MultilevelSmoother(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end

            assert(isempty(blf.b), ...
                'Multilevel smoothers only tested for symmetric problems.')

            obj.nLevels = 0;
            obj.fes = fes;
            obj.blf = blf;

            mesh = fes.mesh;
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end

        function nonInvertedSmootherMatrix = setup(obj, ~)
            obj.nLevels = obj.nLevels + 1;

            % should be implemented by subclasses
            nonInvertedSmootherMatrix = [];
        end
    end
    
    methods (Access=protected)
        % callback function to be executed before mesh refinement:
        % -> patches change for all vertices adjacent to refined edges and
        %    all new vertices
        function getChangedPatches(obj, mesh, bisecData)
            nCOld = mesh.nCoordinates;
            nCNew = mesh.nCoordinates + nnz(bisecData.bisectedEdges) + ...
                bisecData.nRefinedElements' * bisecData.nInnerNodes;
            bisectedEdgeNodes = unique(mesh.edges(:,bisecData.bisectedEdges));
            obj.changedPatches{obj.nLevels+1} = false(nCNew, 1);
            idx = [bisectedEdgeNodes; ((nCOld+1):nCNew)'];
            obj.changedPatches{obj.nLevels+1}(idx) = true;
        end
    end

    %% abstract methods
    methods (Abstract, Access=public)
        smooth(obj, x, k)
        prolongate(obj, x, k)
        restrict(obj, x, k)
    end
end