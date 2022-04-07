% LoL2Prolongation (subclass of Prolongation) Prolongate lowest order L^2 finite
%   element function to refined mesh.
%
%   Pu = LoL2Prolongation(u) returns a handle to the prolongation object.
%
%   Pu.interpolatedData is the data of the prolongation and is set automatically
%       at mesh refinement.

classdef LoL2Prolongation < Prolongation
    %% methods
    methods (Access='public')
        function obj = LoL2Prolongation(fes)
            assert(isa(fes.finiteElement, 'LowestOrderL2Fe'), ...
                'LoL2Prolongation needs a lowest order L2 finite element space.')
            obj = obj@Prolongation(fes);
        end
    end
    
    methods (Access='protected')
        function setupMatrix(obj, mesh, data)
            % use that for lowest order L2 elements the dofs correspond to elements
            % and indices of new elements is known
            nChildren = ones(mesh.nElements, 1);
            for k = 1:data.nBisecGroups
                nChildren(data.bisecGroups{k}.elementIdx) = data.bisecGroups{k}.nDescendants;
            end
            nNewElements = sum(nChildren);
            
            obj.matrix = sparse(1:nNewElements, repelem(1:mesh.nElements, nChildren), ...
                ones(nNewElements,1), nNewElements, mesh.nElements);
        end
    end
end
