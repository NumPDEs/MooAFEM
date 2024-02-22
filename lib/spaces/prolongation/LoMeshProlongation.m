% LoMeshProlongation (subclass of MeshProlongation) Prolongate lowest order L^2/H^1
%   conforming finite element function to refined mesh.
%
%   P = LoMeshProlongation(fes) returns a handle to the prolongation object
%       associated to the finite element space fes. The prolongation matrix
%       P.matrix is set automatically at mesh refinement.
%
% See also Prolongation

classdef LoMeshProlongation < MeshProlongation
    %% properties
    properties (Access=protected)
        feType
    end

    %% methods
    methods (Access=public)
        function obj = LoMeshProlongation(fes)
            obj = obj@MeshProlongation(fes);
            switch class(fes.finiteElement)
                case 'LowestOrderL2Fe'
                    obj.feType = 'L2';
                case 'LowestOrderH1Fe'
                    obj.feType = 'H1';
                case 'LowestOrderCRFe'
                    obj.feType = 'CR';
                otherwise
                    eid = 'LoMeshProlongation:wrongFeType';
                    msg = 'LoMeshProlongation needs a lowest order L2, H1, or CR finite element space.';
                    throwAsCaller(MException(eid, msg));
            end
        end
    end
    
    methods (Access=protected)
        function setupMatrix(obj, mesh, data)
            if isequal(obj.feType, 'L2')
                setupMatrixL2(obj, mesh, data);
            elseif isequal(obj.feType, 'H1')
                setupMatrixH1(obj, mesh, data);
            elseif isequal(obj.feType, 'CR')
                setupMatrixCR(obj, mesh, data);
            end
        end
        
        function setupMatrixL2(obj, mesh, data)
            % use that for lowest order L2 elements the dofs correspond to
            % elements and indices of new elements is known
            nChildren = getNChildrenPerElement(data);
            nNewElements = sum(nChildren);
            
            obj.matrix = sparse(1:nNewElements, repelem(1:mesh.nElements, nChildren), ...
                ones(nNewElements,1), nNewElements, mesh.nElements);
        end


        function setupMatrixCR(obj, mesh, data)
            % use that for lowest order CR elements the dofs correspond to
            % coordinates and new coordinates reside on edges or on inner nodes
            %
            % fixed boundary dofs remain zero 
            dofs = getDofs(obj.fes);
            nEntries = 9 * mesh.nElements + nnz(data.bisectedEdges) + ...
                3 * sum(data.nInnerNodes .* data.nRefinedElements);
            [I, J, V] = deal(zeros(nEntries, 1));

            conFES = FeSpace(obj.fes.mesh, LowestOrderH1Fe());
            conDofs = getDofs(conFES);
            conFreeDofs = getFreeDofs(conFES);
            conFixedDofs = getFixedDofs(conFES);
            
            angles = mesh.computeElementAngles();

            % old nodal S1 dofs on the coarse level get averaged values of CR function
            CRtoS1local = [1 -1 1; 1 1 -1; -1 1 1];
            % TODO: vectorize this!
            for j = 1:mesh.nElements
                idx = (j-1)*9 + (1:9);
                I(idx) = repmat(conDofs.element2Dofs(:,j), 3, 1);
                J(idx) = repelem(dofs.element2Dofs(:,j), 3, 1);
                V(idx) = reshape(CRtoS1local .* angles(:,j) ./ 2 ./ pi, [], 1);
            end
            idxNr = 9 * mesh.nElements;
            dofNr = mesh.nCoordinates;
            
            % new nodal S1 dofs are copied from unique dofs of old CR function
            % NB: This uses that each edge is bisected at most once
            nNewEdgeNodes = nnz(data.bisectedEdges);
            idx = idxNr + (1:nNewEdgeNodes);
            I(idx) = dofNr + transpose(1:nNewEdgeNodes);
            J(idx) = transpose(dofs.edge2Dofs(:,data.bisectedEdges));
            V(idx) = ones(nNewEdgeNodes, 1);
            dofNr = dofNr + nNewEdgeNodes;
            idxNr = idxNr + nNewEdgeNodes;
            
            % inner dofs are weighted sum of element dofs
            for k = find(data.nRefinedElements)'
                nNewInnerNodes = data.nRefinedElements(k)*data.nInnerNodes(k);
                idx = idxNr + (1:3*nNewInnerNodes);
                I(idx) = repelem(dofNr + (1:nNewInnerNodes)', 3);
                J(idx) = reshape(repelem(dofs.element2Dofs, [1;1;1], data.nInnerNodes(k)), [], 1);
                evalCR = transpose(CRtoS1local) * data.bisection{k}.innerNodes;
                V(idx) = reshape(repmat(evalCR, 1, data.nInnerNodes(k)), [], 1);
                idxNr = idxNr + 3*nNewInnerNodes;
                % dofNr = dofNr + nNewInnerNodes;
            end
            
            nNewDofs = mesh.nCoordinates + nnz(data.bisectedEdges) + ...
                sum(data.nInnerNodes.*data.nRefinedElements);
            obj.matrix = sparse(I, J, V, nNewDofs, dofs.nDofs);
            obj.matrix(conFixedDofs, :) = 0;
        end
        
        function setupMatrixH1(obj, mesh, data)
            % use that for lowest order H1 elements the dofs correspond to
            % coordinates and new coordinates reside on edges or on inner nodes
            dofs = getDofs(obj.fes);
            nEntries = mesh.nCoordinates + 2*nnz(data.bisectedEdges) + ...
                3*sum(data.nInnerNodes.*data.nRefinedElements);
            [I, J, V] = deal(zeros(nEntries, 1));
            
            % node dofs stay the same
            idx = 1:mesh.nCoordinates; 
            I(idx) = idx';
            J(idx) = idx';
            V(idx) = ones(idx(end), 1);
            dofNr = mesh.nCoordinates;
            
            % edge dofs are mean of edge dofs
            n = nnz(data.bisectedEdges);
            idx = idx(end) + (1:2*n);
            I(idx) = repelem(dofNr + (1:n)', 2);
            J(idx) = reshape(dofs.edge2Dofs(:,data.bisectedEdges), [], 1);
            V(idx) = ones(2*n, 1)/2;
            dofNr = idx(end) + n;
            
            % inner dofs are weighted sum of element dofs
            idxEnd = idx(end);
            for k = find(data.nRefinedElements)'
                n = data.nRefinedElements(k)*data.nInnerNodes(k);
                idx = idxEnd + (1:3*n);
                I(idx) = repelem(dofNr + (1:n)', 3);
                J(idx) = reshape(repelem(dofs.element2Dofs, [1;1;1], data.nInnerNodes(k)), [], 1);
                V(idx) = reshape(repmat(data.bisection{k}.innerNodes, 1, data.nInnerNodes(k)), [], 1);
                dofNr = dofNr + n;
                idxEnd = idxEnd + 3*n;
            end
            
            nNewDofs = mesh.nCoordinates + nnz(data.bisectedEdges) + ...
                sum(data.nInnerNodes.*data.nRefinedElements);
            obj.matrix = sparse(I, J, V, nNewDofs, dofs.nDofs);
        end

        % no-op override to avoid wrong call to superclass method
        function connectDofs(~, ~, ~)
        end
    end
end
