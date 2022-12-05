% buildLocalMatrixFactorization Patchwise assembly of matrix and cholesky
%   factorization.

function buildLocalMatrixFactorization(obj, data)

% since this operates only on free dofs, store global numbering
freeDofs = getFreeDofs(obj.fes);
dofs = getDofs(obj.fes);
global2freeDofs = zeros(dofs.nDofs, 1);
global2freeDofs(freeDofs) = 1:numel(freeDofs);
obj.patchesAsFreeDofs = cellfun(@(x) global2freeDofs(x), obj.patchDofs, 'UniformOutput', false);
nDofsPerElement = size(dofs.element2Dofs, 1);

% store cholesky decomposition of local matrices associated to patches
obj.patchwiseChol = cell(numel(obj.patchDofs), 1);
for k = 1:numel(obj.patchDofs)
    % get free dofs on patch and all dofs on patch
    pDofs = obj.patchDofs{k};
    n = numel(pDofs);
    pElements = obj.patchElements{k};
    allPatchDofs = dofs.element2Dofs(:,pElements);
    
    % find local numbering of free dofs and relation between all and free dofs
    [isNeeded, localNumbering] = ismember(allPatchDofs, pDofs);
    isNeeded = reshape(isNeeded, size(isNeeded,1), 1, size(isNeeded,2));
    idx = find(and(isNeeded, pagetranspose(isNeeded)));
    
    % assemble arrays for sparse matrix with all dofs on patch ...
    I = repmat(localNumbering, nDofsPerElement, 1);
    J = repelem(localNumbering, nDofsPerElement, 1);
    localData = data(:,pElements);
    
    % ... and assemble/decompose only the parts that correspond to free dofs
    obj.patchwiseChol{k} = ...
        chol(accumarray([I(idx), J(idx)], localData(idx), [n, n]));
end

end