% LoL2Prolongation (subclass of Prolongation) Prolongate lowest order L^2 finite
%   element function to refined mesh.
%
%   Pu = LoL2Prolongation(u) returns a handle to the prolongation object.
%
%   Pu.interpolatedData is the data of the prolongation and is set automatically
%       at mesh refinement.

classdef LoL2Prolongation < Prolongation
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        interpolatedData
    end
    
    %% methods
    methods (Access='public')
        function obj = LoL2Prolongation(feFunction)
            assert(isa(feFunction.fes.finiteElement, 'LowestOrderL2Fe'), ...
                'Prolongation only implemented for lowest order L2 finite elements.')
            obj = obj@Prolongation(feFunction);
        end
    end
    
    methods (Access='protected')
        function interpolate(obj, src, ~)
            % use that for lowest order L2 elements the dofs correspond to elements
            % and indices of new elements is known
            nChildElems = ones(src.nElements, 1);
            for k = 1:length(src.bisecGroups)
                nChildElems(src.bisecGroups{k}.elementIdx) = src.bisecGroups{k}.nDescendants;
            end
            obj.interpolatedData = repelem(obj.u.data, nChildElems);
        end
    end
end
