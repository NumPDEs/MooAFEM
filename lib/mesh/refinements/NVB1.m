% NVB1 (subclass of NVB) Class realizing the NVB1 refinement rule. Only the
%   initial marking step differs from generic NVB.

classdef NVB1 < NVB
    %% methods
    methods (Access=public)
        function obj = NVB1(mesh)
            obj = obj@NVB(mesh);
        end
    end
    
    methods (Access=protected)
        function markedEdges = markedElementsToEdges(obj, markedElements)
            markedEdges = false(1,obj.mesh.nEdges);
            markedEdges(obj.mesh.element2edges(1,markedElements)) = true;
        end
    end
end
