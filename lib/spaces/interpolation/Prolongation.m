% Prolongation (abstract handle class) Interface for prolongation of FeFunctions
%   from coarse to fine mesh.

classdef Prolongation < handle
    %% properties
    properties (Abstract, GetAccess='public', SetAccess='protected')
        interpolatedData
    end
    
    properties (Access='protected')
        u FeFunction
        listenerHandle
    end
    
    %% methods
    methods (Access='public')
        function obj = Prolongation(feFunction)
            obj.u = feFunction;
            mesh = feFunction.fes.mesh;
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.interpolate);
            obj.interpolatedData = [];
        end
    end
    
    methods (Abstract, Access='protected')
        interpolate(obj, src, event)
    end
end
