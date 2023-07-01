% SpaceProlongation (temporary value class) assembles the prolongation
% operation from S^p to S^q with q > p

classdef SpaceProlongation
    
   properties (GetAccess=public, SetAccess=protected)
        matrix
    end
    
    properties (Access=protected)
        listenerHandle
    end
    
    %% methods
    methods (Access=public)
        function obj = SpaceProlongation(sourceFes, targetFes)
            unitTriangleInclusionMatrix = SpaceProlongation.getUnitTriangleInclusionMatrix(sourceFes.finiteElement, targetFes.finiteElement);
            
            % Copy the matrix from the unit triangle for each element
            mat = repmat(unitTriangleInclusionMatrix(:), [1, sourceFes.mesh.nElements]);
            
            % Create index sets for the accumarray
            fromDofs = getDofs(sourceFes);
            toDofs = getDofs(targetFes);
            I = repmat(fromDofs.element2Dofs, size(toDofs.element2Dofs, 1), 1);
            J = repelem(toDofs.element2Dofs, size(fromDofs.element2Dofs, 1), 1);
            
            % assemble matrix and account for multiply computed dofs
            obj.matrix = sparse(I(:), J(:), mat(:), fromDofs.nDofs, toDofs.nDofs);
            dofCounts = accumarray(toDofs.element2Dofs(:), 1, [toDofs.nDofs, 1])';
            obj.matrix = obj.matrix ./ dofCounts;
        end
    end

    methods (Static, Access=protected)
        function matrix = getUnitTriangleInclusionMatrix(sourceFE, targetFE)
            unittriangle = Mesh.unitTriangle();
            fromFes = FeSpace(unittriangle, sourceFE);
            toFes = FeSpace(unittriangle, targetFE);
            Vertices = 1:getDofs(fromFes).nDofs;
            matrix = permute(SpaceProlongation.interpolateData(eye(length(Vertices)), fromFes, toFes), [2 1]);

            % Get local to global map in unit triangle
            matrix = matrix(getDofs(fromFes).element2Dofs, getDofs(toFes).element2Dofs);
        end

        function interpolatedData = interpolateData(data, fromFes, toFes)
            Dofs = 1:getDofs(toFes).nDofs;
            nComponents = size(data, 2);
            interpolatedData = zeros(numel(Dofs), nComponents);
            feFunctionWrapper = FeFunction(fromFes);
            for k = 1:nComponents
                feFunctionWrapper.setData(data(:,k));
                wholeData = nodalInterpolation(feFunctionWrapper, toFes);
                interpolatedData(:,k) = wholeData(Dofs);
            end
        end
    end
end

