% rotateCoordinates Rotate coordinates by given angle in radiants.
%
%   rotate(coordinates, rotationAngle) rotates coordinates around origin.
%
%   rotate(coordinates, rotationAngle, refPoint) rotates coordinates around
%       reference point.

function newCoordinates = rotateCoordinates(coordinates, rotationAngle, refPoint)

arguments
    coordinates (2,:) double
    rotationAngle (1,1) double
    refPoint (2,1) double = [0;0]
end

A = [cos(rotationAngle) -sin(rotationAngle); sin(rotationAngle) cos(rotationAngle)];
newCoordinates = A * (coordinates - refPoint);
newCoordinates = newCoordinates + refPoint;

end