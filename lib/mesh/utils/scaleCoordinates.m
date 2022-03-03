% scaleCoordinates Scale coordinates by given factor which can be a scalar
%   (isotrope scaling), or a 2-vector (scaling in x- and y-direction).
%
%   scale(coordinates, scaleFactor) scales coordinates starting at origin.
%
%   scale(coordinates, scaleFactor, refPoint) scales coordinates starting at a
%       reference point.

function newCoordinates = scaleCoordinates(coordinates, scaleFactor, refPoint)

arguments
    coordinates (2,:) double
    scaleFactor (:,1) double {mustBePositive(scaleFactor)}
    refPoint (2,1) double = [0;0]
end

assert(size(scaleFactor, Dim.Elements) <= 2, 'Scaling factor must be 1- or 2-dimensional.')

newCoordinates = scaleFactor .* (coordinates - refPoint) + refPoint;

end
