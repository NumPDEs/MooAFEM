% translateCoordinates Translate coordinates by given vector.
%
%  translate(coordinates, shiftVector) translates coordinates by shiftVector.

function newCoordinates = translateCoordinates(coordinates, shiftVector)

arguments
    coordinates (2,:) double
    shiftVector (2,1) double
end

newCoordinates = coordinates + shiftVector;

end