% refineLocally Refine a mesh locally based on marked elements/edges.
%
%   obj.refineLocally(marked) refines marked elements with NVB(3). Marked
%       edges/elements can be given by their index.
%
%   obj.refineLocally(marked, method) refines marked elements with one of the
%       following methods:
%       - NVB[1,3,5]: Element based Newest Vertex Bisecton with 1, 3, or 5
%           bisections for the marked triangles.
%       - NVBEdge: Edge based Newest Vertex Bisection.
%       - RGB: Element based Red-Green-Blue refinement.
%
%   See also Mesh

function refineLocally(obj, marked, method)

arguments
    obj (1,1) Mesh
    marked (:,1) double
    method (1,1) string = 'NVB'
end

%% prepare data for refinement and make it available to listeners
refinementAlgorithm = feval(method, obj);
[obj.markedEdges, obj.bisecGroups] = refinementAlgorithm.prepareRefinementData(marked);
obj.notify('IsAboutToRefine')

%%  do actual refinement and clear data of old mesh
obj.updateData();
obj.markedEdges = [];
obj.bisecGroups = [];
obj.trafo = [];
obj.notify('HasChanged');

end