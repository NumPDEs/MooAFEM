% markGoafemFPZ GOAFEM marking from [Feischl, Praetorius, van der Zee; 2016].
%
%   [marked, total] = markGoafemFPZ(eta, zeta, theta) marks elements based on
%       primal indicators eta, dual indicators zeta, and marking parameter
%       theta. Marks equally many elements for eta and zeta.
%
%   [marked, total] = markGoafemFPZ(eta, zeta, theta, C) Mark up to factor
%       C > 0 more elements for larger set.

function [marked, total] = markGoafemFPZ(eta, zeta, theta, C)

arguments
    eta (:,1) double
    zeta (:,1) double
    theta (1,1) double {mustBeInRange(theta, 0, 1)}
    C (1,1) double {mustBeGreaterThan(C, 0)} = 1
end

[markedEta, totalEta] = markDoerflerSorting(eta, theta);
[markedZeta, totalZeta] = markDoerflerSorting(zeta, theta);

total = totalEta*totalZeta;
N = min(numel(markedEta), numel(markedZeta));
M = ceil(C * N);
if numel(markedEta) > numel(markedZeta)
    marked = union(markedZeta(1:N), markedEta(1:min(M, end)));
else
    marked = union(markedEta(1:N), markedZeta(1:min(M, end)));
end

end