% assemble Assembles matrix of bilinear form.
%
%   A = assemble(obj) assembles the data of the given BilinearForm and
%   returns a sparse matrix A.
%
%   See also: BilinearForm

function mat = assemble(obj)
%% setup
if all([isempty(obj.a), isempty(obj.b), isempty(obj.c), isempty(obj.robin)])
    error('BilinearForm:NoCoefficients', 'No coefficients are given.')
end
fes = obj.fes;

%% integrate element data
mat = 0;
phi = TestFunction(fes);
Dphi = TestFunctionGradient(fes);

% diffusion
if ~isempty(obj.a)
    f = CompositeFunction(@BilinearForm.diffusionPart, obj.a, Dphi);
    mat = mat + integrateElement(f, obj.qra);
end

% convection
if ~isempty(obj.b)
    f = CompositeFunction(@BilinearForm.convectionPart, obj.b, Dphi, phi);
    mat = mat + integrateElement(f, obj.qrb);
end

% reaction
if ~isempty(obj.c)
    f = CompositeFunction(@BilinearForm.reactionPart, obj.c, phi);
    mat = mat + integrateElement(f, obj.qrc);
end

%% construct sparse matrix
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
    f = CompositeFunction(@BilinearForm.robinPart, obj.robin, phi);
    idx = getCombinedBndEdges(fes.mesh, obj.bndRobin);
    [I, J] = getLocalDofs(size(dofs.edge2Dofs, Dim.Vector));
    mat = mat + sparse(dofs.edge2Dofs(I,idx), dofs.edge2Dofs(J,idx), ...
        integrateEdge(f, obj.qrRobin, idx), dofs.nDofs, dofs.nDofs);
end

end

function [localI, localJ] = getLocalDofs(n)
    localI = repmat(1:n, 1, n);
    localJ = repelem(1:n, 1, n);
end
