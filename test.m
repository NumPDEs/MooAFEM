mesh = Mesh.loadFromGeometry('unitsquare');
mesh.refineUniform();
trafo = getAffineTransformation(mesh);

fes1 = FeSpace(mesh, HigherOrderH1Fe(4), 'dirichlet', ':');
fes2 = FeSpace(mesh, HigherOrderH1Fe(2), 'dirichlet', ':');

testMatrix = ones(getDofs(fes2).nDofs, 1);
inclusionMatrix = FeProlongation(fes1, fes2);
inclusionMatrix = full(inclusionMatrix);
Vertices = 1:getDofs(fes1).nDofs;
inclusionMatrix2 = interpolateData(eye(length(Vertices)), fes1, fes2)';


function inclusionMatrix = FeProlongation(fromFes, toFes)
    %Create the prolongation matrix on the unit triangle
    unittriangle = Mesh.loadFromGeometry('unittriangle');
    fromFesUnitTriangle = FeSpace(unittriangle, ...
                             HigherOrderH1Fe(fromFes.finiteElement.order));
    toFesUnitTriangle = FeSpace(unittriangle, ...
                               HigherOrderH1Fe(toFes.finiteElement.order));
    Vertices = 1:getDofs(fromFesUnitTriangle).nDofs;
    unitTriangleInclusionMatrix = ...
        permute(interpolateData(eye(length(Vertices)),...
                           fromFesUnitTriangle, toFesUnitTriangle), [2 1]);

    % Get dofs of the finite element spaces on the meshes
    fromFesDofs = getDofs(fromFes);
    toFesDofs = getDofs(toFes);

    % Copy the matrix from the unit triangle for each element
    mat = repmat(unitTriangleInclusionMatrix, ...
                 [1, 1, fromFes.mesh.nElements]);

    % Create index sets for the accumarray
    temp = fromFesDofs.element2Dofs(:, :);
    temp2 = toFesDofs.element2Dofs(:, :);
    I = repmat(permute(temp, [1 3 2]), [1, size(temp2, 1), 1]);
    J = repmat(permute(temp2, [3 1 2]), [size(temp, 1), 1, 1]);
    ind = [I(:), J(:)];
    
    inclusionMatrix = ...
      accumarray(ind, mat(:), [], @mean, [], true);
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


