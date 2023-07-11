% LoFeProlongation (subclass of FeProlongation) Prolongate lowest order L^2/H^1
%   conforming finite element function to refined mesh.
%
%   P = LoFeProlongation(fes) returns a handle to the prolongation object
%       associated to the finite element space fes. The prolongation matrix
%       P.matrix is set automatically at mesh refinement.
%
% See also Prolongation

classdef LoFeProlongation < FeProlongation
    %% properties
    properties (Access=protected)
        feType
    end

    %% methods
    methods (Access=public)
        function obj = LoFeProlongation(fes)
            obj = obj@FeProlongation(fes);
            switch class(fes.finiteElement)
                case 'LowestOrderL2Fe'
                    obj.feType = 'L2';
                case 'LowestOrderH1Fe'
                    obj.feType = 'H1';
                otherwise
                    eid = 'LoFeProlongation:wrongFeType';
                    msg = 'LoFeProlongation needs a lowest order L2 or H1 finite element space.';
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
