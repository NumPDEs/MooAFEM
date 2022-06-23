% Barycentric2D (subclass of Barycentric) Specialization for Barycentric
%   coordinates on triangles.
%
%   bary2D = Barycentric2D.from1D(bary1D, edgeInd) constructs Barycentric2D from
%       Barycentric1D on edge given by edgeInd.

classdef Barycentric2D < Barycentric
    properties (GetAccess=public, SetAccess=protected)
        dim = 2;
    end
    
    methods (Static)
        function bary2D = from1D(bary1D, edgeInd)
            arguments
                bary1D (1,1) Barycentric1D
                edgeInd double {mustBeInRange(edgeInd, 1, 3)}
            end
            
            newNodes = zeros(3, bary1D.nNodes);
            s = [1,2; 2,3; 3,1];
            
            % assign 1D barycentric coordinates to edgeInd-th edge in triangle
            newNodes(s(edgeInd,:),:) = bary1D.coordinates;
            bary2D = Barycentric2D(newNodes);
        end
    end
end
