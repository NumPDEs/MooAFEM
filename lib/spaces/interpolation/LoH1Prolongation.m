% LoH1Prolongation (subclass of Prolongation) Prolongate lowest order H^1
%   conforming finite element function to refined mesh.
%
%   P = LoH1Prolongation(fes) returns a handle to the prolongation object
%       associated to the finite element space fes. The prolongation matrix
%       P.matrix is set automatically at mesh refinement.
%
%   prolongate(P, u) returns the prolongated data of FeFunction u.

classdef LoH1Prolongation < Prolongation
    %% methods
    methods (Access='public')
        function obj = LoH1Prolongation(fes)
            assert(isa(fes.finiteElement, 'LowestOrderH1Fe'), ...
                'LoH1Prolongation needs a lowest order H1 finite element space.')
            obj = obj@Prolongation(fes);
        end
    end
    
    methods (Access='protected')
        function setupMatrix(obj, mesh, data)
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
                V(idx) = reshape(repmat(data.bisecMethods{k}.innerNodes, 1, data.nInnerNodes(k)), [], 1);
                dofNr = dofNr + n;
                idxEnd = idxEnd + 3*n;
            end
            
            nNewDofs = mesh.nCoordinates + nnz(data.bisectedEdges) + ...
                sum(data.nInnerNodes.*data.nRefinedElements);
            obj.matrix = sparse(I, J, V, nNewDofs, dofs.nDofs);
        end
    end
end
