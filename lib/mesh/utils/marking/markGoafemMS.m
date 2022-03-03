% markGoafemMS GOAFEM marking from [Mommer, Stevenson; 2009].
%
%   [marked, total] = markGoafemMS(eta, zeta, theta) marks elements based on
%       primal indicators eta, dual indicators zeta, and marking parameter
%       theta.

function [marked, total] = markGoafemMS(eta, zeta, theta)

arguments
    eta (:,1) double
    zeta (:,1) double
    theta (1,1) double {mustBeInRange(theta, 0, 1)}
end

[markedEta, totalEta] = markDoerflerSorting(eta, theta);
[markedZeta, totalZeta] = markDoerflerSorting(zeta, theta);

total = totalEta*totalZeta;
if numel(markedEta) > numel(markedZeta)
    marked = markedZeta;
else
    marked = markedEta;
end

end