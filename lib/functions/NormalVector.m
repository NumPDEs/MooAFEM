% NormalVector (subclass of Evaluable) Encapsulate edge normal vector within the
%   uniform interface of Evaluable for use in integration routines.
%
%   f = NormalVector(mesh, funcHandle) represents the normal vector on the given
%       mesh.
%
%   See also: Evaluable

classdef NormalVector < Evaluable
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        mesh
    end
    
    %% methods
    methods (Access='public')
        function obj = NormalVector(mesh)
            obj.mesh = mesh;
        end
        
        function val = eval(obj, ~, ~)
            error('Normal vector cannot be evaluated on elements.')
        end
        
        function val = evalEdge(obj, ~, idx)
            arguments
                obj
                ~
                idx {mustBeIndexVector} = ':'
            end
            
            trafo = getAffineTransformation(obj.mesh);
            val = trafo.n(:,idx);
        end
    end
end
