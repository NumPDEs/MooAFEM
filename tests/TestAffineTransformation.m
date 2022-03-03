classdef TestAffineTransformation < matlab.unittest.TestCase
    
properties
    mesh
    trafo
end

methods (TestMethodSetup)
    function createDataForCustomMesh(testCase)
        coordinates = [0,0;1,0;0,1;0.5,0.5;0.5,0.0]';
        elements = [1,4,3;1,5,4;5,2,4]';
        boundary = {[1,5;5,2;2,4;4,3;3,1]'};
        testCase.mesh = Mesh(coordinates, elements, boundary);
        testCase.trafo = getAffineTransformation(testCase.mesh);
    end
end

methods (Test)
    function obj = elementDeterminantsSumToOneHalf(testCase)
        testCase.verifyEqual(sum(testCase.trafo.area), 1/2, 'AbsTol', eps)
    end

    function obj = totalEdgeLengthIsCorrect(testCase)
        expected = (5 + 3*sqrt(2)) / 2;
        testCase.verifyEqual(sum(testCase.trafo.ds), expected, 'AbsTol', eps)
    end

    function inverseTrafoTransformsToReferenceTriangle(testCase)
        z1 = testCase.mesh.coordinates(:,testCase.mesh.elements(1,:));
        z2 = testCase.mesh.coordinates(:,testCase.mesh.elements(2,:));
        z3 = testCase.mesh.coordinates(:,testCase.mesh.elements(3,:));
        vectors = [z2 - z1; z3 - z1];
        actual = vectorProduct(testCase.trafo.DFinv, vectors, [2,2], [2,2]);
        threeTimesTref = repmat([1;0;0;1], 1, 3);
        testCase.verifyEqual(actual, threeTimesTref, 'AbsTol', eps)
    end

    function checkUnitNormalProperties(testCase)
        totalLength = sum(sqrt(vectorProduct(testCase.trafo.n, testCase.trafo.n)));
        testCase.verifyEqual(totalLength, testCase.mesh.nEdges, 'AbsTol', eps)
        z1 = testCase.mesh.coordinates(:,testCase.mesh.edges(1,:));
        z2 = testCase.mesh.coordinates(:,testCase.mesh.edges(2,:));
        edgeVectors = z2 - z1;
        totalOrthogonality = sum(sqrt(vectorProduct(testCase.trafo.n, edgeVectors)));
        testCase.verifyEqual(totalOrthogonality, 0, 'AbsTol', eps)  
    end
end

end
