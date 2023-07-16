% MultilevelSmoother (abstract handle class) Interface for smoothers to use in a
%   multigrid solver.
%
% Abstract methods to implement:
%
% smoother.setup(A) hook to set up the smoother for the matrix A.
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

    %% methods
    methods (Access=public)
        function obj = MultilevelSmoother()
            obj.nLevels = 0;
        end

        function setup(obj, ~)
            obj.nLevels = obj.nLevels + 1;
        end
    end

    %% abstract methods
    methods (Abstract, Access=public)
        smooth(obj, x, k)
        prolongate(obj, x, k)
        restrict(obj, x, k)
    end
end