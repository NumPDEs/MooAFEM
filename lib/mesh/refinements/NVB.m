% NVB Class realizing the NVB3 refinement rule, which has the typical
%   computational flow for NVB-like refinement strategies (hence, the general
%   name 'NVB').
%
%   This implementation is taken from [Funken, Praetorius, Wissgott; 2011]
% 	>> Efficient Implementation of Adaptive P1-FEM in Matlab <<

classdef NVB
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        mesh
    end
    
    %% methods
    methods (Access='public')
        function obj = NVB(mesh)
            obj.mesh = mesh;
        end
        
        function [markedEdges, bisecGroups] = prepareRefinementData(obj, markedElements)
            markedEdges = obj.markedElementsToEdges(markedElements);
            markedEdges = obj.meshClosure(markedEdges);
            bisecGroups = obj.groupElements(markedEdges, markedElements);
        end
    end
    
    methods (Access='protected')
        function markedEdges = markedElementsToEdges(obj, markedElements)
            markedEdges = false(obj.mesh.nEdges,1);
            markedEdges(obj.mesh.element2edges(:,markedElements)) = true;
        end
        
        function markedEdges = meshClosure(obj, markedEdges)
            hasHangingNodes = 1;
            while nnz(hasHangingNodes) > 0
                edge = reshape(markedEdges(obj.mesh.element2edges), size(obj.mesh.element2edges));
                hasHangingNodes = ( ~ edge(1,:) & (edge(2,:) | edge(3,:)) );
                markedEdges(obj.mesh.element2edges(1,hasHangingNodes)) = true;
            end
        end
        
        function bisecGroups = groupElements(obj, markedEdges, ~)
            % group elements according to their refinement pattern
            edge = reshape(markedEdges(obj.mesh.element2edges), size(obj.mesh.element2edges));
            bisecGroups{1} = Bisec1  ( find( edge(1,:) & ~edge(2,:) & ~edge(3,:) ) );
            bisecGroups{2} = Bisec12 ( find( edge(1,:) &  edge(2,:) & ~edge(3,:) ) );
            bisecGroups{3} = Bisec13 ( find( edge(1,:) & ~edge(2,:) &  edge(3,:) ) );
            bisecGroups{4} = Bisec123( find( edge(1,:) &  edge(2,:) &  edge(3,:) ) );
        end
    end
end
