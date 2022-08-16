% intergridMatrix reusing refineLocally to
% construct intergrid operators for P1 only
% %Refine a mesh locally based on marked elements/edges.
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

function intergridMatrix(obj, marked, method)


arguments
    obj (1,1) Mesh
    marked (:,1) double
    method (1,1) string = 'NVB'
end


%% prepare data for refinement 
refinementAlgorithm = feval(method, obj);
bisecData = refinementAlgorithm.prepareRefinementData(marked);

nNodes = length(obj.coordinates);
nodesBisecEdges = obj.edges(:,bisecData.bisectedEdges);%which edges were refined

%%  do actual refinement and clear data of old mesh
obj.updateData(bisecData);

%% compute intergrid matrix
nupdatedNodes = length(obj.coordinates);
newnodestwice = repmat(nNodes+1:nupdatedNodes,2,1);

obj.locVert{end+1} = [unique(nodesBisecEdges);[nNodes+1:nupdatedNodes]'];
oldneigh = nodesBisecEdges; 
oldneigh = reshape(oldneigh', 2*length(oldneigh),1);

Int_I = [[1:nNodes]'; newnodestwice(:)];
Int_J = [[1:nNodes]'; nodesBisecEdges(:)];
Int_val = [ones(nNodes,1); 1/2*ones(length(nodesBisecEdges(:)),1)];
Int = sparse(Int_I, Int_J, Int_val);
obj.intergrid{end+1} = Int;


obj.trafo = [];
obj.notify('JustRefined');
obj.notify('RefineCompleted')

end