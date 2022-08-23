classdef TestMesh < matlab.unittest.TestCase
    
methods (Test)
    function predefinedMeshesAreValidMeshes(testCase)
        rootDir = fileparts(fileparts(mfilename('fullpath')));
        geometryDir = fullfile(rootDir, 'lib', 'mesh', '@Mesh', 'geometries');
        folders = dir(geometryDir);
        folders = {folders([folders(:).isdir]).name};
        folders = {folders{~startsWith(folders, '.')}};
        for i = 1:numel(folders)
            testCase.verifyClass(Mesh.loadFromGeometry(folders{i}), 'Mesh')
        end
    end
    
    function cloneMethodMakesClone(testCase)
        mesh = Mesh.loadFromGeometry('Lshape');
        clonedMesh = clone(mesh);
        testCase.verifyEqual(mesh.coordinates, clonedMesh.coordinates);
        testCase.verifyEqual(mesh.elements, clonedMesh.elements);
        testCase.verifyEqual(mesh.edges, clonedMesh.edges);
        testCase.verifyEqual(mesh.element2edges, clonedMesh.element2edges);
        testCase.verifyEqual(mesh.flipEdges, clonedMesh.flipEdges);
        testCase.verifyEqual(mesh.boundaries, clonedMesh.boundaries);
    end

    function LShapeDetailsAreCorrect(testCase)
        % this is tailored to the Lshape to ensure correctness of arrays
        mesh = Mesh.loadFromGeometry('Lshape');
        testCase.verifyEqual(mesh.nCoordinates, 11);
        testCase.verifyEqual(mesh.nElements, 12);
        testCase.verifyEqual(mesh.nEdges, 22);
        testCase.verifyEqual(mesh.nBoundaries, 2);
        testCase.verifyNumElements(mesh.boundaries{1}, 2)
        testCase.verifyNumElements(mesh.boundaries{2}, 6)
    end

    function edgeMidpointsAreSymmetricInSquare(testCase)
        bary = Barycentric1D([1;1]/2);
        mesh = Mesh.loadFromGeometry('unitsquare');
        centerOfMass = mean(edgewiseCoordinates(mesh, bary), Dim.Elements);
        testCase.verifyEqual(centerOfMass, [0.5;0.5], 'AbsTol', eps)
        centerOfMass = mean(edgewiseCoordinates(mesh, bary, [1,2,4,6]), Dim.Elements);
        testCase.verifyEqual(centerOfMass, [0.5;0.5], 'AbsTol', eps)
    end

    function elementMidpointsAreSymmetricInSquare(testCase) 
        bary = Barycentric2D([1;1;1]/3);
        mesh = Mesh.loadFromGeometry('unitsquare');
        centerOfMass = mean(elementwiseCoordinates(mesh, bary), Dim.Elements);
        testCase.verifyEqual(centerOfMass, [0.5;0.5], 'AbsTol', eps)
        centerOfMass = mean(elementwiseCoordinates(mesh, bary, [1,3]), Dim.Elements);
        testCase.verifyEqual(centerOfMass, [0.5;0.5], 'AbsTol', eps)
        centerOfMass = mean(elementwiseCoordinates(mesh, bary, logical([0,1,0,1])), Dim.Elements);
        testCase.verifyEqual(centerOfMass, [0.5;0.5], 'AbsTol', eps)
    end
    
    function patchAreasAreThreeTimesDomainSize(testCase)
        mesh = Mesh.loadFromGeometry('unitsquare');
        mesh.refineUniform(3);
        patchAreas = computePatchAreas(mesh);
        testCase.verifyEqual(sum(patchAreas), 3, 'AbsTol', eps);
    end
end
    
end
