
N = 10;

% rng(10);
eta = abs(fix(sqrt(N))*randn(N, 1));

theta = 0.3;


% Doerfler marking by sorting
eMarked = markDoerflerSorting(eta, theta);

eNotSoOptimalMarked = markDoerflerBinning(eta, theta);

% Doerfler marking by QuickMark
% eQuickMarked = markDoerflerQuick(eta, theta);
eVeryQuickMarked = quickMark(eta, theta);

eEvenQuickerMarked = quickMarkOptimized(eta, theta);

assert(all(sort(eMarked) == sort(eVeryQuickMarked)));