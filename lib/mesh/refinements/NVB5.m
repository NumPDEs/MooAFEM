% NVB5 (subclass of NVB) Class realizing the NVB5 refinement rule. This differs
%   from generic NVB by the refinement of the marked elements, where an interior
%   node is introduced.

classdef NVB5 < NVB
    %% methods
    methods (Access=public)
        function obj = NVB5(mesh)
            obj = obj@NVB(mesh);
        end
    end
    
    %% methods
    methods (Access=protected)
        function markedEdges = markedElementsToEdges(obj, markedElements)
            markedEdges = markedElementsToEdges@NVB(obj, markedElements);
        end
        
        function bisecData = groupElements(obj, markedEdges, markedElements)
            % split elements with 3 marked edges into 'has interior node' and 'doesnt'
            hasInteriorNode = false(1,obj.mesh.nElements);
            hasInteriorNode(markedElements) = true;
            
            edge = reshape(markedEdges(obj.mesh.element2edges), size(obj.mesh.element2edges));
            bisecData = BisectionEventData(obj.mesh.nElements, markedEdges, ...
                Bisec1,   find(edge(1,:) & ~edge(2,:) & ~edge(3,:)), ...
            	Bisec12,  find(edge(1,:) &  edge(2,:) & ~edge(3,:)), ...
            	Bisec13,  find(edge(1,:) & ~edge(2,:) &  edge(3,:)), ...
            	Bisec123, find(edge(1,:) &  edge(2,:) &  edge(3,:) & ~hasInteriorNode), ...
            	Bisec5,   find(edge(1,:) &  edge(2,:) &  edge(3,:) &  hasInteriorNode));
        end
    end
end
