% buildLocalMatrixFactorization Patchwise assembly of matrix and cholesky
%   factorization.

function buildLocalMatrixFactorization(obj, data)

% since this operates only on free dofs, store global numbering
freeDofs = getFreeDofs(obj.fes);
obj.global2freeDofs = zeros(getDofs(obj.fes).nDofs, 1);
obj.global2freeDofs(freeDofs) = 1:numel(freeDofs);

% store cholesky decomposition of local matrices associated to patches
% TODO: this is the most time intensive line in the code due
% to sparse matrix indexing. Is there a better way?
obj.patchwiseChol = cellfun(@(p) chol(full(data(p,p))), obj.patchDofs, 'UniformOutput', false);

end