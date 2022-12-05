% mldivide Solve linear System A * x = y for patchwise matrix A. This means
%   solving local problems on patches associated to vertices 
%   (depending on level) and add the patchwise solutions
%   together. mldivide(A, y) can also be called by A \ y.

function x = mldivide(obj, y)

x = zeros(size(y));
lower = struct('UT', true, 'TRANSA', true); 
upper = struct('UT', true);

for k = obj.activePatches
    U = obj.patchwiseChol{k};
    idx = obj.patchesAsFreeDofs{k};

    % solve local patch problems (update is additive)
    update = linsolve(U, linsolve(U, y(idx,:), lower), upper);
    x(idx,:) = x(idx,:) + update;
end

end
