
N = 1e2;

eta = rand(N, 1);

theta = 0.3;


% Doerfler marking by sorting
[eta, Ind] = sort(eta, 'descend');
cumsumEta = cumsum(eta);
nMarked = find(cumsumEta >= theta * cumsumEta(end), 1, 'first');
eMarked = Ind(1:nMarked);

% Doerfler marking by QuickMark
eQuickMarked = markDoerflerQuick(eta, theta);

assert(all(eMarked == eQuickMarked));