classdef TestCoordinateTransformations < matlab.unittest.TestCase
    
properties
    coordinates
end

methods (TestMethodSetup)
    function createUnitTriangleCoordinates(testCase)
        mesh = Mesh.unitTriangle();
        testCase.coordinates = mesh.coordinates;
        testCase.fatalAssertNumElements(testCase.coordinates, 6);
    end
end

methods (Test)
    function validateTranslation(testCase)
        newCoordinates = translateCoordinates(testCase.coordinates, [1,2]);
        testCase.verifyEqual(newCoordinates, [1,2;2,2;1,3]', 'AbsTol', eps)
    end

    function validateRotation(testCase)
        newCoordinates = rotateCoordinates(testCase.coordinates, pi/2, [1,1]/2);
        testCase.verifyEqual(newCoordinates, [1,0;1,1;0,0]', 'AbsTol', eps)
    end

    function validateScaling(testCase)
        newCoordinates = scaleCoordinates(testCase.coordinates, 2, [1,1]);
        testCase.verifyEqual(newCoordinates, [-1,-1;1,-1;-1,1]', 'AbsTol', eps)
    end

    function validateSkewX(testCase)
        newCoordinates = skewCoordinates(testCase.coordinates, pi/4, 'x', [1,0]);
        testCase.verifyEqual(newCoordinates, [0,0;1,0;-1,1]', 'AbsTol', eps)
    end

    function validateSkewY(testCase)
        newCoordinates = skewCoordinates(testCase.coordinates, pi/4, 'y', [0,2]);
        testCase.verifyEqual(newCoordinates, [0,0;1,1;0,1]', 'AbsTol', eps)
    end
end
    
end