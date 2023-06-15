mesh = Mesh.loadFromGeometry('unitsquare');
mesh.refineUniform();
trafo = getAffineTransformation(mesh);

fes1 = FeSpace(mesh, HigherOrderH1Fe(2), 'dirichlet', ':');
fes2 = FeSpace(mesh, HigherOrderH1Fe(4), 'dirichlet', ':');

testMatrix = ones(getDofs(fes1).nDofs, 1);
inclusionMatrix = FeProlongation(fes1, fes2);

function inclusionMatrix = FeProlongation(fromFes, toFes)
    unittriangle = Mesh.loadFromGeometry('unittriangle');
    fromFesUnitTriangle = FeSpace(unittriangle, HigherOrderH1Fe(fromFes.finiteElement.order));
    toFesUnitTriangle = FeSpace(unittriangle, HigherOrderH1Fe(toFes.finiteElement.order));
    Vertices = 1:getDofs(fromFesUnitTriangle).nDofs;
    unitTriangleInclusionMatrix = interpolateData(eye(length(Vertices)), ...
                                   fromFesUnitTriangle, toFesUnitTriangle);

    fromFesDofs = getDofs(fromFes);
    toFesDofs = getDofs(toFes);

    [row, ~] = getLocalDofs(size(fromFesDofs.element2Dofs, Dim.Vector));
    [~, column] = getLocalDofs(size(toFesDofs.element2Dofs, Dim.Vector));

    I = fromFesDofs.element2Dofs(row, :);
    J = toFesDofs.element2Dofs(:, column);
    sub2ind(size(I), I(:, :));
    
    inclusion = accumarray([I, J], unitTriangleInclusionMatrix, [], @mean, [], true);
end

function [localI, localJ] = getLocalDofs(n)
    localI = repmat(1:n, 1, n);
    localJ = repelem(1:n, 1, n);
end

function interpolatedData = interpolateData(data, fromFes, toFes)
    Dofs = 1:getDofs(toFes).nDofs;
    nComponents = size(data, 2);
    interpolatedData = zeros(numel(Dofs), nComponents);
    feFunctionWrapper = FeFunction(fromFes);
    for k = 1:nComponents
        feFunctionWrapper.setData(data(:,k));
        wholeData = nodalInterpolation(feFunctionWrapper, toFes);
        interpolatedData(:,k) = wholeData(Dofs);
    end
end
