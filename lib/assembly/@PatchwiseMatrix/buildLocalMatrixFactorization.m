% buildLocalMatrixFactorization Patchwise assembly of matrix and cholesky
%   factorization.

function buildLocalMatrixFactorization(obj, data)

freeDofs = getFreeDofs(obj.fes);
global2freeDofs = zeros(getDofs(obj.fes).nDofs, 1);
global2freeDofs(freeDofs) = 1:numel(freeDofs);
% store cholesky decomposition of local matrices associated
% to patches (only upper triangle)
obj.patchDofs = cellfun(@(p) global2freeDofs(p), obj.patchDofs, 'UniformOutput', false);

% TODO: this is the most time intensive line in the code due
% to sparse matrix indexing. Is there a better way?
obj.patchwiseChol = cellfun(@(p) full(chol(data(p,p))), obj.patchDofs, 'UniformOutput', false);

end