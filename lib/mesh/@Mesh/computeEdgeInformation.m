% computeEdgeInformation Compute edge information: edges, element2edges, and
%   boundary edges.
%
%   [edges, element2edges, flipEdges, boundary2edges] = ...
%       Mesh.computeEdgeInformation(elements) computes edge information from
%       elements only. The array boundary2edges is empty.
%
%   [edges, element2edges, flipEdges, boundary2edges] = ...
%       Mesh.computeEdgeInformation(elements, boundaries) computes edge
%       information and takes also boundary parts into account.
%
%   See also: Mesh

function [edges, element2edges, flipEdges, boundary2edges] = computeEdgeInformation(elements, boundaries)

arguments
    elements (3,:) double
    boundaries cell {mustBeBoundaryArrays} = {}
end

% get all edges (and boundary information if present)
nEdges = 3*size(elements, Dim.Elements);
edges = [reshape(elements([1,2,2,3,3,1],:),2,[]), horzcat(boundaries{:})]';

% find unique edges
reverse = edges(:,1) > edges(:,2);
edges(reverse,:) = fliplr(edges(reverse,:));
[edges, ia, ie] = unique(edges, 'rows');
edges = edges';

% use information to build arrays
element2edges = reshape(ie(1:nEdges),3,[]);
pointer = cumsum([nEdges, cellfun(@(x) size(x, Dim.Elements), boundaries)']);
leftIndices = num2cell( pointer(1:(end-1)) +  1 );
rightIndices = num2cell( pointer(2:end) );
boundary2edges = cellfun(@(l,r) ie(l:r), leftIndices, rightIndices, 'UniformOutput', false);

% preserve original orientation of edges (at least on boundary)
reverseUnique = reverse(ia);
edges(:,reverseUnique) = flipud(edges(:,reverseUnique));
reverse = xor(reverse, reverseUnique(ie));
flipEdges = reshape(reverse(1:nEdges), 3, []);

end


% nested function that hecks if input is cell array of boundaries
function mustBeBoundaryArrays(a)
    arguments
        a cell
    end

    % boundary part must be either of size (2,:) or empty
    if ~all(cellfun(@(x) size(x, Dim.Vector)==2 | isempty(x), a))
        eid = 'mustBeBoundaryArrays:BoundaryNotValid';
        msg = 'Boundary parts must be arrays of size (2,:).';
        throwAsCaller(MException(eid,msg))
    end
end