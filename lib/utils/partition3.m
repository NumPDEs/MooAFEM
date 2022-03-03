% partition3 Compute all partitions of given integer into 3 non-negative numbers.
%
%   p = partition3(n) returns all sets [k1, k2, k3] of non-negative numbers such
%       that n = sum([k1, k2, k3]). The sets are the columns of p.

function res = partition3(n)

arguments
    n (1,1) double
end

[a, b] = find(flipud(triu(ones(n+1))));
res = [(n+1)-b, (n+1)-a, a+b-(n+2)]';

end