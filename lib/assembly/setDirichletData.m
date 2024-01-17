% setDirichletData Includes inhomogeneous Dirichlet data into linear system
%   by nodal interpolation in all fixed boundary nodes
%
%   [rhs, uDirData] = setDirichletData(uDir, fes, A, F) evaluates the given
%       function handle uDir in all fixed dofs of the FeSpace fes and
%       corrects computes the right-hand side F of the linear system
%
%   See also: nodalInterpolation

function [rhs, uDirData] = setDirichletData(uDir, fes, A, F)

arguments
    uDir function_handle
    fes (1,1) FeSpace
    A double
    F double
end

% interpolate Dirichlet data on boundary elements(!)
uDirData = nodalInterpolation(MeshFunction(fes.mesh, uDir), fes, ...
                                           fes.mesh.boundaryElems);

% vanish coefficients inside of the domain
freeDofs = getFreeDofs(fes);
uDirData(freeDofs) = 0;

% update right-hand side by lifting
rhs = F - A * uDirData;

% restrict Dirichlet-data interpolation to fixed dofs
uDirData(freeDofs) = [];

end