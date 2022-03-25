% RGB (subclass of NVB) Class realizing the RGB refinement rule. This differs
%   from generic NVB only by the refinement of elements with three marked edges
%   and the choice of the refinement edge (longest edge).

classdef RGB < NVB
    %% methods
    methods (Access='public')
        function obj = RGB(mesh)
            obj = obj@NVB(mesh);
        end
    end
    
    %% methods
    methods (Access='protected')
        function markedEdges = markedElementsToEdges(obj, markedElements)
            % sort edges by length (while still retaining the orientation of the element)
            mesh = obj.mesh;
            trafo = getAffineTransformation(mesh);
            [~, maxIdx] = max(trafo.ds(mesh.element2edges), [], 1);
            mesh.changeRefinementEdge(maxIdx);
            
            markedEdges = false(mesh.nEdges,1);
            markedEdges(mesh.element2edges(:,markedElements)) = true;
        end
        
        function bisecGroups = groupElements(obj, markedEdges, markedElements)
            bisecGroups = groupElements@NVB(obj, markedEdges, markedElements);
            bisecGroups{4} = BisecRed(bisecGroups{4}.elementIdx) ;
        end
    end
end
