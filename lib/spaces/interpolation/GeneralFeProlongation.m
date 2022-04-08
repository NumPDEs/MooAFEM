% LinearFeProlongation (subclass of Prolongation) Provides piecewise linear
%   interpolant of FeFunction on a refined mesh.
%
%   Pu = LinearFeProlongation(u) returns the piecewise linear interpolation of u
%       that can be prolongated to a finer mesh.
%
%   Pu.interpolatedData is the data of the prolongation and is set automatically
%       at mesh refinement.

classdef LinearFeProlongation < Prolongation
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        interpolatedData
    end
    
    properties (Access='private')
        interpolant
    end
    
    %% methods
    methods (Access='public')
        function obj = LinearFeProlongation(feFunction)
            obj = obj@Prolongation(feFunction);
            obj.interpolant = [];
        end
    end
    
    methods
        function val = get.interpolatedData(obj)
            if ~isempty(obj.interpolant)
                obj.interpolatedData = nodalInterpolation(obj.interpolant, obj.u.fes);
                obj.interpolant = [];
            end
            val = obj.interpolatedData;
        end
    end
    
    methods (Access='protected')
        function interpolate(obj, ~, ~)
            obj.interpolant = LinearInterpolation(obj.u);
        end
    end
end
