% nodalInterpolation Compute nodal interpolant of given function.
%
%   val = nodalInterpolation(f, fes) compute data for nodal interpolation of f
%       into FeSpace fes.

function val = nodalInterpolation(f, fes)

arguments
    f Evaluable
    fes FeSpace
end

fe = fes.finiteElement;

assert(isa(fe, 'NodalFiniteElement'), ...
    'Nodal interpolation defined only for nodal finite elements.')

bary = getDofLocations(fe);
localMatrix = squeeze(evalShapeFunctions(fe, bary));
feval = eval(f, bary);

assert(isequal(size(feval, Dim.Vector), 1), ...
    'Only scalar function can be interpolated.')

feval = repeatAlongDimension(feval, Dim.QuadratureNodes, bary.nNodes);
elementValues = (shiftdim(feval,1) / localMatrix)';
elementValues = repeatAlongDimension(elementValues, Dim.Elements, fes.mesh.nElements);

dofs = getDofs(fes);
val = zeros(dofs.nDofs, 1);
val(dofs.element2Dofs) = elementValues;

end

function val = repeatAlongDimension(val, dim, n)
    if size(val, dim) ~= n
        repetitions = ones(1, dim);
        repetitions(dim) = n;
        val = repmat(val, repetitions);
    end
end
