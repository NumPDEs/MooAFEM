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
            fe = obj.fes.finiteElement;
            dofLocations = getDofLocations(fe);
            nLocalDofs = dofLocations.nNodes;
            oldDofs = getDofs(obj.fes);
            
            % get old-elements-to-new-dofs connectivity
            nChildElems = ones(mesh.nElements, 1);
            for k = 1:data.nBisecGroups
                nChildElems(data.bisecGroups{k}.elementIdx) = data.bisecGroups{k}.nDescendants;
            end
            element2newDof = cumsum([1; nChildElems*nLocalDofs]);
            
            % pre-allocate
            n = nnz(nChildElems == 1)*nLocalDofs + sum(nChildElems(nChildElems > 1))*nLocalDofs^2;
            [I, J, V] = deal(zeros(n, 1));
            
            % for non-refined elements, just transfer the dofs
            elems = find(nChildElems == 1);
            n = numel(elems)*nLocalDofs;
            idx = (1:n)';
            I(idx) = getConsecutiveIndices(element2newDof(elems), nLocalDofs);
            J(idx) = reshape(oldDofs.element2Dofs(:,elems), [], 1);
            V(idx) = 1;
            ptr = n;
            
            % handle refined elements
            localMatrix = squeeze(evalShapeFunctions(fe, dofLocations));
            for k = 1:data.nBisecGroups
                if data.bisecGroups{k}.nMembers ~= 0
                    unitTriangle = getBisectedUnitTriangle(class(data.bisecGroups{k}));
                    newDofLocations = getElementwiseLocationsAsBarycentric(unitTriangle, dofLocations);
                    newDofWeights = localMatrix \ squeeze(evalShapeFunctions(fe, newDofLocations));

                    elems = data.bisecGroups{k}.elementIdx;
                    n = data.bisecGroups{k}.nMembers * numel(newDofWeights);
                    idx = ptr + (1:n)';
                    I(idx) = repelem(getConsecutiveIndices(element2newDof(elems), newDofLocations.nNodes), nLocalDofs);
                    J(idx) = reshape(repelem(oldDofs.element2Dofs(:,elems), 1, newDofLocations.nNodes), [], 1);
                    V(idx) = reshape(repmat(newDofWeights, 1, data.bisecGroups{k}.nMembers), [], 1);
                    ptr = ptr + n;
                end
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
            mesh.refineLocally(1, 'NVBEdge');
        case 'Bisec12'
            mesh.refineLocally([1,2], 'NVBEdge');
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

function elementwiseLocations = getElementwiseLocationsAsBarycentric(mesh, locations)
    elementwiseLocations = reshape(elementwiseCoordinates(mesh, locations), 2, [], 1);
    elementwiseLocations = Barycentric2D(...
        cat(Dim.Vector, elementwiseLocations, 1-sum(elementwiseLocations, Dim.Vector)));
end

function idx = getConsecutiveIndices(parents, nIndices)
    idx = reshape((0:(nIndices-1)) + parents, [], 1);
end