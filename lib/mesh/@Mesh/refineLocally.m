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

function refineLocally(obj, marked, method, bool)



arguments
    obj (1,1) Mesh
    marked (:,1) double
    method (1,1) string = 'NVB'
    bool (1,1) = 1
end



%% prepare data for refinement and make it available to listeners
refinementAlgorithm = feval(method, obj);
bisecData = refinementAlgorithm.prepareRefinementData(marked);

% if bool == 1
%     bisecEdgescopy = bisecData.bisectedEdges;
%     bisecEdges = bisecEdgescopy;
%     bool = 0;
% end


obj.notify('IsAboutToRefine', bisecData);


% nNodes = length(obj.coordinates);
% nodesBisecEdges = obj.edges(:,bisecEdges);%which edges were refined

%%  do actual refinement and clear data of old mesh
obj.updateData(bisecData);

% %% compute intergrid matrix
% nupdatedNodes = length(obj.coordinates);
% newnodestwice = repmat(nNodes+1:nupdatedNodes,2,1);
% oldneigh = nodesBisecEdges; %edge2nodes(edge2newNode~=0,:);
% oldneigh = reshape(oldneigh', 2*length(oldneigh),1);
% Int_I = [[1:nNodes]'; newnodestwice(:)];
% Int_J = [[1:nNodes]'; nodesBisecEdges(:)];
% Int_val = [ones(nNodes,1); 1/2*ones(length(nodesBisecEdges(:)),1)];
% Int = sparse(Int_I, Int_J, Int_val);
% obj.intergrid{obj.level} = Int;


obj.trafo = [];
obj.notify('JustRefined');
obj.notify('RefineCompleted')

end