% assemble Assembles vector of linear form.
%
%   F = assemble(obj) assembles the data of the given LinearForm and
%   returns a vector F.
%
%   See also: LinearForm

function vec = assemble(obj)
%% setup
if all([isempty(obj.f), isempty(obj.fvec), isempty(obj.neumann), isempty(obj.robin)])
    error('LinearForm:NoCoefficients', 'No coefficients are given.')
end
fes = obj.fes;

%% integrate element data
vec = 0;
phi = TestFunction(fes);
Dphi = TestFunctionGradient(fes);

% scalar rhs
if ~isempty(obj.f)
    f = CompositeFunction(@LinearForm.scalarPart, obj.f, phi);
    vec = vec + integrateElement(f, obj.qrf);
end

% vector rhs
if ~isempty(obj.fvec)
    f = CompositeFunction(@LinearForm.vectorPart, obj.fvec, Dphi);
    vec = vec + integrateElement(f, obj.qrfvec);
end

%% construct rhs vector from element data
dofs = getDofs(fes);
vec = accumarray(dofs.element2Dofs(:), vec(:), [dofs.nDofs, 1]);

%% integrate and add boundary data
% Robin data
if ~isempty(obj.robin)
    f = CompositeFunction(@LinearForm.robinPart, obj.robin, phi);
    idx = getCombinedBndEdges(fes.mesh, fes.bnd.robin);
    edgeData = integrateEdge(f, obj.qrRobin, idx);
    edge2Dofs = dofs.edge2Dofs(:,idx);
    vec = vec + accumarray(edge2Dofs(:), edgeData(:), [dofs.nDofs, 1]);
end

% Neumann data
if ~isempty(obj.neumann)
    f = CompositeFunction(@LinearForm.neumannPart, obj.neumann, phi);
    idx = getCombinedBndEdges(fes.mesh, fes.bnd.neumann);
    edgeData = integrateEdge(f, obj.qrNeumann, idx);
    edge2Dofs = dofs.edge2Dofs(:,idx);
    vec = vec + accumarray(edge2Dofs(:), edgeData(:), [dofs.nDofs, 1]);
end

end
