% assemble Assembles vector of linear form.
%
%   F = assemble(obj, fes) assembles the data of the given LinearForm for a
%       FeSpace fes and returns a vector F.
%
%   See also: LinearForm

function vec = assemble(obj, fes)
arguments
    obj
    fes FeSpace
end
%% setup
if all([isempty(obj.f), isempty(obj.fvec), isempty(obj.neumann), isempty(obj.robin)])
    error('LinearForm:NoCoefficients', 'No coefficients are given.')
end

%% integrate element data
vec = 0;
phi = TestFunction(fes);
Dphi = TestFunctionGradient(fes);
p = fes.finiteElement.order;

% scalar rhs
if ~isempty(obj.f)
    qr = QuadratureRule.ifEmptyOfOrder(obj.qrf, p);
    f = CompositeFunction(@LinearForm.scalarPart, obj.f, phi);
    vec = vec + integrateElement(f, qr);
end

% vector rhs
if ~isempty(obj.fvec)
    qr = QuadratureRule.ifEmptyOfOrder(obj.qrfvec, p-1);
    f = CompositeFunction(@LinearForm.vectorPart, obj.fvec, Dphi);
    vec = vec + integrateElement(f, qr);
end

%% construct rhs vector from element data
dofs = getDofs(fes);
vec = accumarray(dofs.element2Dofs(:), vec(:), [dofs.nDofs, 1]);

%% integrate and add boundary data
% Robin data
if ~isempty(obj.robin)
    qr = QuadratureRule.ifEmptyOfOrder(obj.qrRobin, p, "1D");
    f = CompositeFunction(@LinearForm.robinPart, obj.robin, phi);
    idx = getCombinedBndEdges(fes.mesh, fes.bnd.robin);
    edgeData = integrateEdge(f, qr, idx);
    edge2Dofs = dofs.edge2Dofs(:,idx);
    vec = vec + accumarray(edge2Dofs(:), edgeData(:), [dofs.nDofs, 1]);
end

% Neumann data
if ~isempty(obj.neumann)
    qr = QuadratureRule.ifEmptyOfOrder(obj.qrNeumann, p-1, "1D");
    f = CompositeFunction(@LinearForm.neumannPart, obj.neumann, phi);
    idx = getCombinedBndEdges(fes.mesh, fes.bnd.neumann);
    edgeData = integrateEdge(f, qr, idx);
    edge2Dofs = dofs.edge2Dofs(:,idx);
    vec = vec + accumarray(edge2Dofs(:), edgeData(:), [dofs.nDofs, 1]);
end

end
