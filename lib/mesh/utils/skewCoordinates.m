% skewCoordinates Skew coordinates by given angle in radiants along x- or
%   y-axis. The use of this method should be well-considered, because skewing
%   can lead to degenerate triangles.
%
%   skew(coordinates, skewAngle) skews coordinates around origin along x-axis.
%
%   skew(coordinates, skewAngle, XorY) skews coordinates around origin and
%       choose axis ('x', or 'y').
%
%   skew(coordinates, skewAngle, XorY, refPoint) skews coordinates around
%       reference point and choose axis.

function newCoordinates = skewCoordinates(coordinates, skewAngle, XorY, refPoint)

arguments
    coordinates (2,:) double
    skewAngle (1,1) double {mustBeInRange(skewAngle, 0, 1.57)}
    XorY (1,1) char {mustBeMember(XorY, {'x', 'y'})} = 'x'
    refPoint (2,1) double = [0;0]
end

if XorY == 'x', A = [1 -tan(skewAngle); 0 1]; end
if XorY == 'y', A = [1 0; tan(skewAngle) 1]; end
newCoordinates = A * (coordinates - refPoint);
newCoordinates = newCoordinates + refPoint;

end