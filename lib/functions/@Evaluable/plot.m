% plot Plot Evaluable object.
%
%   plot(obj) Plot function with predefined options.
%
%   plot(obj, 'SubdivisionFactor', sf) Subdivides all triangles that have
%       area larger than 1/sf of area of whole domain. This is used to
%       accurately plot non-linear functions over coarse meshes.
%
%   plot(__, options) Plot function with some additional options given as
%       name-value pairs. These options can include any of the options for
%       the patch function. Note that not all of these patch related
%       options are tested for compatibility with the plot method and the
%       options are discarded if a 'SubdivisionFactor' is given!
%
%   See also: trisurf, patch

function plot(obj, customOptions, plotOptions)

arguments
    obj (1,1) Evaluable
    customOptions.SubdivisionFactor (1,1) double = 1;
    plotOptions.?matlab.graphics.primitive.Patch
end

%% preallocate axes handles
data = eval(obj, Barycentric2D([1;0;0]), 1);
nComponents = size(data, Dim.Vector);
ax = cell(nComponents, 1);
for i = 1:nComponents
    ax{i} = axes(gcf);
end

%% determine large triangles
mesh = obj.mesh;
trafo = mesh.getAffineTransformation();
subdivision = ceil(sqrt(customOptions.SubdivisionFactor * trafo.area ./ ...
                        sum(trafo.area, Dim.Elements)));
[smin, smax] = bounds(subdivision);

%% plot data for each subdivision
plotOptions = namedargs2cell(plotOptions);
if isequal(smin, smax, 1)
    % for no subdivision, patch can be also used to plot lines
    %   -> options can be piped through    
    bary = Barycentric2D(partition3(1));
    coordinates = mesh.elementwiseCoordinates(bary);
    Xel = splitFirstDimension(coordinates, [1,2,3]);
    data = guaranteeSize(eval(obj, bary), mesh.nElements, bary.nNodes);
    Zel = splitFirstDimension(data, [1,2,3]);

    for i = 1:nComponents
        patch(ax{i}, Xel{1}, Xel{2}, Zel{i}, Zel{i}, plotOptions{:})
        hold on
    end
else
    % for subdivisions, patch/line must be invoked separately
    %   -> discard custom options
    if ~isempty(plotOptions)
        warning('Custom plot options discarded due to subdivisions.')
    end
    for s = smin:smax
        currentElems = find(subdivision == s);
        if ~isempty(currentElems)
            % subtriangles on one element
            bary = Barycentric2D(partition3(s)/s);
            [elements, outerEdges] = getSubtriangles(s);
            
            % compute plot data
            coordinates = mesh.elementwiseCoordinates(bary, currentElems);
            Xel = splitFirstDimension(coordinates, elements);
            Xlin = splitFirstDimension(coordinates, outerEdges);
            data = guaranteeSize(eval(obj, bary, currentElems), numel(currentElems), bary.nNodes);
            Zel = splitFirstDimension(data, elements);
            Zlin = splitFirstDimension(data, outerEdges);
            
            for i = 1:nComponents
                patch(ax{i}, Xel{1}, Xel{2}, Zel{i}, Zel{i}, 'LineStyle', 'none')
                hold on
                line(ax{i}, Xlin{1}, Xlin{2}, Zlin{i}, 'Color', 'k');
            end
        end
    end
end

%% make output look nice
variableName = inputname(1);
for i = 1:nComponents
    titleString = sprintf('Plot of %s', variableName); 
    if nComponents > 1
        titleString = sprintf([titleString, ' (component %d)'], i);
    end
    title(ax{i}, titleString)
    grid(ax{i}, 'on')
    view(ax{i}, 3)
end

end

%% local functions
% to understand how this (admittedly somewhat arcane) piece of code works,
% make a sketch of subtriangles and number the nodes according to the order
% given by partition3
function [subtriangles, outerEdges] = getSubtriangles(n)
    N = n*(n+1)/2;
    st1 = cumsum([1:N; repelem(1:n, 1:n); ones(1,N)], 1);
    k = n-1;
    K = k*(k+1)/2;
    st2 = cumsum([setdiff(1:N, cumsum([1, 1:k])); -ones(1,K); repelem(2+(1:k), 1:k)], 1);
    subtriangles = [st1, st2]';
    
    oe = [cumsum([1, 1:k]), cumsum(1:n), N + (1:n)];
    outerEdges = cumsum([oe; 1:n, (1:n)+1, ones(1,n)], 1)';
end

% compute components of data fed to plot commands
function splitArray = splitFirstDimension(array, idx)
    n = size(array, Dim.Vector);
    k = size(idx, 2);
    splitArray = cell(n,1);
    nodalData = array(:,:,idx);
    for i = 1:n
        splitArray{i} = reshape(nodalData(i,:,:), [], k)';
    end
end

% 
function data = guaranteeSize(data, n, k)
    if size(data, Dim.Elements) == 1
        data = data.*ones(1, n, k);
    end 
end
