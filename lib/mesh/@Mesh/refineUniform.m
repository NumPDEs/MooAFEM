% refineUniform Refine a mesh uniformly n times.
%
%   obj.refineUniform() refines all elements with NVB(3)
%
%   obj.refineUniform(n) refines all elements n times with NVB(3)
%
%   obj.refineUniform(n, method) refines all elements/edges n times with one of
%       the following methods:
%       - NVB{1,3,5}: Element based Newest Vertex Bisecton with 1, 3, or 5
%           bisections or all triangles.
%       - NVBEdge: Edge based Newest Vertex Bisection.
%       - RGB: Element based Red-Green-Blue refinement.
%
%   See also Mesh


function refineUniform(obj, n, method)

arguments
    obj (1,1) Mesh
    n (1,1) double = 1
    method (1,1) string = 'NVB'
end

%% iterate refinements
for i=1:n
    % mark all elements/edges & refine
    if strcmp(method, 'NVBEdge') == true
        marked = 1:obj.nEdges;
    else
        marked = 1:obj.nElements;
    end
    obj.refineLocally(marked, method);
end

end
