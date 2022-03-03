% Barycentric Encapsulates a collection of barycentric coordinates.
%
%   bary = Barycentric(coordinates) constructs barycentric coordinates from the
%       (d+1)xN array coordinates, where N is the number of barycentric
%       coordinates and d is the number of spatial dimensions.
%
%   coordinates = arangeAlongDim(bary, dim) arranges the N barycentric
%       coordinates along the dimension given by dim > 2. The resulting array is
%       of dimension (d+1)x1x...x1xN, where N stands on the dim-th position.

classdef Barycentric
    %% properties
    properties (GetAccess = 'public', SetAccess = 'protected')
        coordinates
        nNodes
    end
    
    properties (Abstract, GetAccess = 'public', SetAccess = 'protected')
        dim
    end
    
    %% methods
    methods (Access = 'public')
        function obj = Barycentric(coordinates)
            arguments
                coordinates {Barycentric.verify}
            end
            
            obj.coordinates = coordinates;
            obj.nNodes = size(coordinates, Dim.Elements);
            assert(isequal(obj.dim, size(coordinates, Dim.Vector)-1), ...
                ['Barycentric coordinates must have dimension ', num2str(obj.dim)]);
        end
        
        function coordinates = arangeAlongDim(obj, dim)
            mustBeGreaterThanOrEqual(dim, 3)
            vectorizedSize = [obj.dim+1, ones(1,dim-2), obj.nNodes];
            coordinates = reshape(obj.coordinates, vectorizedSize);
        end
    end
    
    methods (Static)
        function verify(a)
            if ~isnumeric(a) ...
                    || any(a < 0, 'all') ...
                    || any(abs(sum(a, Dim.Vector) - 1) > 1e-8, 'all')
                eid = 'Barycentric:notBarycentric';
                msg = 'Barycentric coordinates must be positive and sum to 1.';
                throwAsCaller(MException(eid,msg))
            end
        end
    end
end
