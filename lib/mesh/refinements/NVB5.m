% NVB5 (subclass of NVB) Class realizing the NVB5 refinement rule. This differs
%   from generic NVB by the refinement of the marked elements, where an interior
%   node is introduced.

classdef NVB5 < NVB
    %% methods
    methods (Access='public')
        function obj = NVB5(mesh)
            obj = obj@NVB(mesh);
        end
    end
    
    %% methods
    methods (Access='protected')
        function markedEdges = markedElementsToEdges(obj, markedElements)
            markedEdges = markedElementsToEdges@NVB(obj, markedElements);
        end
        
        function bisecGroups = groupElements(obj, markedEdges, markedElements)
            bisecGroups = groupElements@NVB(obj, markedEdges);
            
            % split elements with 3 marked edges into 'has interior node' and 'doesnt'
            noInteriorNode = false(1,obj.mesh.nElements);
            noInteriorNode(bisecGroups{4}.elementIdx) = true;
            noInteriorNode(markedElements) = false;
            bisecGroups{4} = Bisec123(find(noInteriorNode));
            bisecGroups{5} = Bisec5(markedElements);
        end
    end
end
