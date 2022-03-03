classdef TestRefinement < matlab.unittest.TestCase
    
properties
    mesh
end
    
properties (TestParameter)
    method = {'NVB', 'NVB1', 'NVB3', 'NVB5', 'NVBEdge', 'RGB'};
end

methods (TestMethodSetup)
    function createMesh(testCase)
        testCase.mesh = Mesh.loadFromGeometry('unitsquare');
    end
end

methods (Test)
    function predefinedMeshesAreValid(testCase)
        testCase.verifyEdgeConnectivityIntegrity(testCase.mesh);
    end
    
    function localRefinementIsValid(testCase, method)
        testCase.mesh.refineLocally([1], method);
        testCase.verifyEdgeConnectivityIntegrity(testCase.mesh);
    end
    
    function uniformRefinementIsValid(testCase, method)
        testCase.mesh.refineUniform(1, method);
        testCase.verifyEdgeConnectivityIntegrity(testCase.mesh);
    end
end

methods (Access='protected')
    function verifyEdgeConnectivityIntegrity(testCase, mesh)
        [edgesRef, e2eRef, ~, bndRef] = Mesh.computeEdgeInformation( ...
            mesh.elements, cellfun(@(bnd) mesh.edges(:,bnd), mesh.boundaries, 'UniformOutput', false));

        % edges are the same and properly connected
        [edgesA, ~, ia] = unique(sort(mesh.edges, 1)', 'rows');
        testCase.verifyEqual(edgesA', sort(edgesRef, 1))
        testCase.verifyEqual(ia(mesh.element2edges), e2eRef)
        for i = 1:numel(bndRef)
            testCase.verifyEqual(sort(ia(mesh.boundaries{i})), sort(bndRef{i}))
        end

        % edge orientation is correct: edge vectors of an element must sum to zero
        edgeVec = mesh.coordinates(:,mesh.edges(2,:)) - mesh.coordinates(:,mesh.edges(1,:));
        edgeVec = edgeVec(:,mesh.element2edges);
        edgeVec(:,mesh.flipEdges) = -edgeVec(:,mesh.flipEdges);
        edgeVec = reshape(edgeVec, 2, 3, []);
        summedEdgeVectors = sum(abs(sum(edgeVec, 2)), 3);
        testCase.verifyEqual(summedEdgeVectors, [0;0], 'AbsTol', eps)
    end
end

end