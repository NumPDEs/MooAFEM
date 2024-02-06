% TangentialVector (subclass of Evaluable) Encapsulate edge tangential vector within the
%   uniform interface of Evaluable for use in integration routines.
%
%   f = TangentialVector(mesh, funcHandle) represents the tangential vector on the given
%       mesh.
%
%   See also: Evaluable

classdef TangentialVector < Evaluable
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        mesh
    end
    
    %% methods
    methods (Access=public)
        function obj = TangentialVector(mesh)
            obj.mesh = mesh;
        end
        
        function val = eval(obj, ~, ~)
            error('Tangential vector cannot be evaluated on elements.')
        end
        
        function val = evalEdge(obj, ~, idx)
            arguments
                obj
                ~
                idx {mustBeIndexVector} = ':'
            end
            
            trafo = getAffineTransformation(obj.mesh);
            val = [-1; 1] .* trafo.n([2;1],idx);
        end
    end
end
