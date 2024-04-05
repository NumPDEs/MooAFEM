function plotCoordinates(obj, idx)
    
arguments
    obj (1,1) Mesh
    idx (:,1) {mustBeIndexVector} = ':'
end

selectedCoordinates = obj.coordinates(:,idx);
scatter(selectedCoordinates(1,:), selectedCoordinates(2,:), 'filled', ...
        'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [0.8 0 0]);

end
