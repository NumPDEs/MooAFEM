% markGoafemBET GOAFEM marking from [Becker, Estecahandy, Trujillo; 2011].
%
%   [marked, total] = markGoafemBET(eta, zeta, theta) marks elements based on
%       primal indicators eta, dual indicators zeta, and marking parameter theta.

function [marked, total] = markGoafemBET(eta, zeta, theta)

arguments
    eta (:,1) double
    zeta (:,1) double
    theta (1,1) double {mustBeInRange(theta, 0, 1)}
end

rho = eta*sum(zeta) + sum(eta)*zeta;
[marked, total] = markDoerflerSorting(rho, theta);

end