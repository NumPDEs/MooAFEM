% plot Plot Mesh object.
%
%   plot(obj) Plot mesh with some predefined options
%
%   plot(obj, options) plots mesh with some additional options given as
%       name-value pairs. These options can include any of the options for the
%       patch function. Note that not all of these patch related options are
%       tested for compatibility with the plot method!
%
%       Additionally, the following options can be given:
%       - 'labelCoordinates', true: Print coordinate labels near coordinates.
%       - 'labelElements', true: Print element labels at element midpoints.
%       - 'labelEdges', true: Print edge labels at edge midpoints.
%       - 'colorBoundaries', true: Color boundary edges for different boundary
%                                  parts differently.
%
%   See also Mesh, trisurf, patch

function plot(obj, labelOptions, plotOptions)
arguments
    obj (1,1) Mesh
    labelOptions.labelCoordinates (1,1) logical = false
    labelOptions.labelElements (1,1) logical = false
    labelOptions.labelEdges (1,1) logical = false
    labelOptions.colorBoundaries (1,1) logical = false
    plotOptions.?matlab.graphics.primitive.Patch
end
    
%% plot mesh and make plot look '2D'
plotOptions = namedargs2cell(plotOptions);
X = reshape(obj.coordinates(1,obj.elements), 3, []);
Y = reshape(obj.coordinates(2,obj.elements), 3, []);
patch(X, Y, 'black', 'facealpha', 0.1, plotOptions{:})
grid on
axis equal

%% add custom plot options
hold on
if labelOptions.colorBoundaries
    cmap = colormap('lines');
    for i = 1:obj.nBoundaries
        color = cmap(i,:);
        bndNodes = obj.edges(:,obj.boundaries{i});
        x = reshape(obj.coordinates(1,bndNodes), 2, []);
        y = reshape(obj.coordinates(2,bndNodes), 2, []);
        line(x, y, 'LineWidth', 2, 'color', color)
    end
end
if labelOptions.labelCoordinates
    text(obj.coordinates(1,:), obj.coordinates(2,:), num2str((1:obj.nCoordinates)'), 'HorizontalAlignment', 'center')
end
if labelOptions.labelElements
    elementMidpoint = Barycentric2D([1;1;1]/3);
    locations = elementwiseCoordinates(obj, elementMidpoint);
    text(locations(1,:), locations(2,:), num2str((1:obj.nElements)'), 'HorizontalAlignment', 'center')
end
if labelOptions.labelEdges
    edgeMidpoint = Barycentric1D([1;1]/2);
    locations = edgewiseCoordinates(obj, edgeMidpoint);
    text(locations(1,:), locations(2,:), num2str((1:obj.nEdges)'), 'HorizontalAlignment', 'center')
end

hold off

end
