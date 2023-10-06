% LinearInterpolation (subclass of Evaluable) Provide piecewise linear
%   interpolant of function. This can be used for approximate prolongation.

classdef LinearInterpolation < Evaluable
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        mesh
        u
    end
    
    properties (Access=private)
        interpolant
    end
    
    %% methods
    methods (Access=public)
        function obj = LinearInterpolation(u)
            arguments
                u Evaluable
            end
            
            obj.mesh = u.mesh;
            obj.u = u;
            meshNodes = Barycentric2D(eye(3));
            vals = zeros(obj.mesh.nCoordinates, 1);
            vals(obj.mesh.elements') = shiftdim(eval(obj.u, meshNodes), 1);
            obj.interpolant = scatteredInterpolant(obj.mesh.coordinates', vals, 'linear');
        end
        
        function val = eval(obj, bary, idx)
            arguments
                obj
                bary Barycentric2D
                idx {mustBeIndexVector}= ':'
            end
            
            assert(~isempty(obj.interpolant), 'No data computed yet.')
            
            coordinates = elementwiseCoordinates(obj.mesh, bary, idx);
            dim = size(coordinates);
            dim(Dim.Vector) = 1;
            val = reshape(obj.interpolant(reshape(coordinates, 2, [])'), dim);
        end
    end
end
