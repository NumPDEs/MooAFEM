% Constant (subclass of Evaluable) Encapsulate constant functions.
%
%   f = Constant(mesh, value) represents the constant value on a given mesh.
%
%   See also: Evaluable

classdef Constant < Evaluable
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        mesh
        constantValue (:,1) double
    end
    
    %% methods
    methods (Access='public')
        function obj = Constant(mesh, value)
            obj.mesh = mesh;
            obj.constantValue = value;
        end
        
        function val = eval(obj, ~, ~)
            val = obj.constantValue;
        end
        
        function val = evalEdge(obj, ~, ~)
            val = obj.constantValue;
        end
    end
end
