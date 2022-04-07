% LoH1Prolongation (subclass of Prolongation) Prolongate lowest order H^1
%   conforming finite element function to refined mesh.
%
%   Pu = LoH1Prolongation(u) returns a handle to the prolongation object.
%
%   Pu.interpolatedData is the data of the prolongation and is set automatically
%       at mesh refinement.

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
            
            nNewDofs = mesh.nCoordinates + nnz(data.bisectedEdges) + ...
                sum(cellfun(@(x) x.nInnerNodes*x.nMembers, data.bisecGroups));
            dofs = getDofs(obj.fes);
            
            [I, J, V] = deal(zeros(nNewDofs, 1));
            
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
            for k = 1:data.nBisecGroups
                nInnerDofs = data.bisecGroups{k}.nInnerNodes;
                n = data.bisecGroups{k}.nMembers*nInnerDofs;
                idx = idxEnd + (1:3*n);
                I(idx) = repelem(dofNr + (1:n)', 3);
                J(idx) = reshape(repelem(dofs.element2Dofs, [1;1;1], nInnerDofs), [], 1);
                V(idx) = reshape(repmat(data.bisecGroups{k}.innerNodes, 1, nInnerDofs), [], 1);
                dofNr = dofNr + n;
                idxEnd = idxEnd + 3*n;
            end
            
            obj.matrix = sparse(I, J, V, nNewDofs, dofs.nDofs);
        end
    end
end
