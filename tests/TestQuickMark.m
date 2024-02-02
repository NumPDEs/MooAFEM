
N = 1e6;

% rng(10);
eta = abs(fix(sqrt(N))*randn(N, 1));

theta = 0.3;


% Doerfler marking by sorting
[etaSorting, Ind] = sort(eta, 'descend');
cumsumEta = cumsum(etaSorting);
nMarked = find(cumsumEta >= theta * cumsumEta(end), 1, 'first');
eMarked = Ind(1:nMarked);

% Doerfler marking by QuickMark
eQuickMarked = markDoerflerQuick(eta, theta);

assert(all(sort(eMarked) == sort(eQuickMarked)));