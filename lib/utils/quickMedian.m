function [m, position] = quickMedian(x, varargin)
    % if nargin < 2
    %     varargin = {};
    % end
    % Compute median (MATLAB built-in with linear complexity)
    % m = median(x, varargin{:});
    m = matlab.internal.math.columnmedian(x, false);
    % Choose larger one of two values closest to the median
    position = find(abs(abs(x - m) - min(abs(x - m))) < 1e-12);
    [~, idx] = min(x(position));
    position = position(idx);

end