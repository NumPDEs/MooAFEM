% NVBEdge (subclass of NVB) Class realizing an edge driven version of the NVB1
%   refinement rule. This rule only differs from generic NVB in the closure
%   step.

classdef NVBEdge < NVB
    %% methods
    methods (Access='public')
        function obj = NVBEdge(mesh)
            obj = obj@NVB(mesh);
        end
    end
    
    %% internal methods
    methods (Access='protected')
        function markedEdges = markedElementsToEdges(obj, markedElements)
            markedEdges = false(1,obj.mesh.nEdges);
            markedEdges(markedElements) = true;
        end
    end
end
