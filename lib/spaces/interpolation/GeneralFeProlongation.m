% GeneralFeProlongation (subclass of Prolongation) Provides prolongation
%   operator for general FeFunction to a refined mesh.
%
%   P = GeneralFeProlongation(fes) returns a handle to the prolongation object
%       associated to the finite element space fes. The prolongation matrix
%       P.matrix is set automatically at mesh refinement.
%
%   prolongate(P, u) returns the prolongated data of FeFunction u.

classdef GeneralFeProlongation < Prolongation
    %% properties
    properties (Access='private')
        postRefineListener
    end
    
    %% methods
    methods (Access='public')
        function obj = GeneralFeProlongation(fes)
            obj = obj@Prolongation(fes);
            obj.postRefineListener = fes.mesh.listener('RefineCompleted', @obj.connectDofs);
            
            assert(isa(obj.fes.finiteElement, 'NodalFiniteElement'), ...
                'Prolongation defined only for nodal finite elements.')
        end
    end
    
    methods (Access='protected')
        function setupMatrix(obj, mesh, data)
            % general idea: compute *all* dofs for each new element and store
            % them consecutively
            % those are connected to the actual new dofs when their numbering
            % per element is available (-> obj.connectDofs)
            
            % de-allocate for memory-efficiency
            obj.matrix = [];
            
            % get dof information
            oldDofs = getDofs(obj.fes);
            fe = obj.fes.finiteElement;
            dofLocations = getDofLocations(fe);
            nLocalDofs = dofLocations.nNodes;
            
            % pre-allocate
            element2newDof = cumsum([1; getNChildrenPerElement(data)*nLocalDofs]);
            n = data.nNonRefinedElements*nLocalDofs + sum(data.nRefinedElements)*nLocalDofs^2;
            [I, J, V] = deal(zeros(n, 1));
            
            % for non-refined elements, just transfer the dofs
            n = data.nNonRefinedElements*nLocalDofs;
            if n ~= 0
                idx = (1:n)';
                nonRefinedElems = getRefinedElementIdx(data, 0);
                I(idx) = getConsecutiveIndices(element2newDof(nonRefinedElems), nLocalDofs);
                J(idx) = reshape(oldDofs.element2Dofs(:,nonRefinedElems), [], 1);
                V(idx) = 1;
            end
            ptr = n;
            
            % handle refined elements (nodal interpolation on new elements)
            localMatrix = squeeze(evalShapeFunctions(fe, dofLocations));
            for k = find(data.nRefinedElements ~= 0)'
                unitTriangle = getBisectedUnitTriangle(class(data.bisection{k}));
                newDofLocations = getNewLocationsElementwise(unitTriangle, dofLocations);
                newDofWeights = divideElementwise(...
                    reshape(evalShapeFunctions(fe, newDofLocations), nLocalDofs, nLocalDofs, []), localMatrix);

                nNewLocs = newDofLocations.nNodes;
                elems = getRefinedElementIdx(data, k);
                n = data.nRefinedElements(k) * numel(newDofWeights);
                idx = ptr + (1:n)';
                I(idx) = repelem(getConsecutiveIndices(element2newDof(elems), nNewLocs), nLocalDofs);
                J(idx) = reshape(repelem(oldDofs.element2Dofs(:,elems), 1, nNewLocs), [], 1);
                V(idx) = reshape(repmat(newDofWeights, 1, data.nRefinedElements(k)), [], 1);
                ptr = ptr + n;
            end
            
            nNewDofs = element2newDof(end)-1;
            obj.matrix = sparse(I, J, V, nNewDofs, oldDofs.nDofs);
        end
    
        function connectDofs(obj, ~, ~)
            dofs = getDofs(obj.fes);
            nOutput = size(obj.matrix, 1);
            idx = zeros(dofs.nDofs, 1);
            
            idx(dofs.element2Dofs(:)) = 1:nOutput;
            obj.matrix = obj.matrix(idx,:);
        end
    end
end

%% local functions
function mesh = getBisectedUnitTriangle(bisecMethod)
    mesh = Mesh.loadFromGeometry('unittriangle');

    % bisect unit triangle according to given bisection rule
    % TODO: this should be handled by bisection method itself?
    switch bisecMethod
        case 'Bisec1'
            mesh.refineLocally(3, 'NVBEdge');
        case 'Bisec12'
            mesh.refineLocally([2,3], 'NVBEdge');
        case 'Bisec13'
            mesh.refineLocally([1,3], 'NVBEdge');
        case 'Bisec123'
            mesh.refineUniform(1, 'NVB');
        case 'Bisec5'
            mesh.refineUniform(1, 'NVB5');
        case 'BisecRed'
            mesh.refineUniform(1, 'RGB');
        otherwise
            error(['Bisection method ', bisecMethod, ' not related to refinement of unit triangle.'])
    end
end

function elementwiseLocations = getNewLocationsElementwise(mesh, locations)
    % collect dofs such that all locations of one child appear before those of other children
    % -> permute dimensions: 2 = children, 3 = barycentric coordinates
    elementwiseLocations = permute(elementwiseCoordinates(mesh, locations), [1,3,2]);
    elementwiseLocations = reshape(elementwiseLocations, 2, [], 1);
    elementwiseLocations = Barycentric2D(...
        cat(Dim.Vector, elementwiseLocations, 1-sum(elementwiseLocations, Dim.Vector)));
end

function C = divideElementwise(A, B)
    % TODO: from 2022a on, this can be replaced by pagemrdivide(A, B)
    n = size(A,1);
    C = B' \ reshape(permute(A, [2,1,3]), n, []);
    C = reshape(permute(reshape(C, n, n, []), [2,1,3]), n, []);
end

function idx = getConsecutiveIndices(parents, nIndices)
    idx = reshape(((0:(nIndices-1)) + parents)', [], 1);
end