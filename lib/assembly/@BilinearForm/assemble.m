% assemble Assembles matrix of bilinear form.
%
%   A = assemble(obj, fes) assembles the data of the given BilinearForm for a
%       FeSpace fes and returns a sparse matrix A.
%
%   See also: BilinearForm

function mat = assemble(obj, fes)
arguments
    obj
    fes FeSpace
end
%% construct sparse matrix
mat = computeVolumeData(obj, fes);
dofs = getDofs(fes);

if ~isequal(mat, 0)
    [I, J] = getLocalDofs(size(dofs.element2Dofs, Dim.Vector));
    mat = sparse(dofs.element2Dofs(I,:), dofs.element2Dofs(J,:), mat, ...
        dofs.nDofs, dofs.nDofs);
else
    mat = sparse(dofs.nDofs, dofs.nDofs);
end

%% integrate and add boundary data
% Robin data
if ~isempty(obj.robin)
    qrRobin = QuadratureRule.ifEmptyOfOrder(obj.qrRobin, 2*fes.finiteElement.order, "1D");
    phi = TestFunction(fes);
    f = CompositeFunction(@BilinearForm.robinPart, obj.robin, phi);
    idx = getCombinedBndEdges(fes.mesh, fes.bnd.robin);
    [I, J] = getLocalDofs(size(dofs.edge2Dofs, Dim.Vector));
    mat = mat + sparse(dofs.edge2Dofs(I,idx), dofs.edge2Dofs(J,idx), ...
        integrateEdge(f, qrRobin, idx), dofs.nDofs, dofs.nDofs);
end

end

function [localI, localJ] = getLocalDofs(n)
    localI = repmat(1:n, 1, n);
    localJ = repelem(1:n, 1, n);
end
