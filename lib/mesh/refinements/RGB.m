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
            sortedEdges = mesh.element2edges;
            for k = 2:3
                idx = find(maxIdx == k);
                sortedEdges(:,idx) = circshift(sortedEdges(:,idx), 4-k, 1);
            end
            
            % initial marking (a bit of a hack to carry sortedEdges to closure step..)
            markedEdges = {false(mesh.nEdges,1), sortedEdges};
            markedEdges{1}(sortedEdges(:,markedElements)) = true;
        end
        
        function markedEdges = meshClosure(~, markedEdges)
            hasHangingNodes = 1;
            [markedEdges, sortedEdges] = deal(markedEdges{:}); % unpack data
            while nnz(hasHangingNodes) > 0
                edge = markedEdges(sortedEdges);
                hasHangingNodes = ( ~ edge(1,:) & (edge(2,:) | edge(3,:)) );
                markedEdges(sortedEdges(1,hasHangingNodes)) = true;
            end
        end
        
        function bisecGroups = groupElements(obj, markedEdges, markedElements)
            bisecGroups = groupElements@NVB(obj, markedEdges, markedElements);
            bisecGroups{4} = BisecRed(bisecGroups{4}.elementIdx) ;
        end
    end
end
