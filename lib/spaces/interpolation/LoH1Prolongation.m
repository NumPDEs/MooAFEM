% LoH1Prolongation (subclass of Prolongation) Prolongate lowest order H^1
%   conforming finite element function to refined mesh.
%
%   Pu = LoH1Prolongation(u) returns a handle to the prolongation object.
%
%   Pu.interpolatedData is the data of the prolongation and is set automatically
%       at mesh refinement.

classdef LoH1Prolongation < Prolongation
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        interpolatedData
    end
    
    %% methods
    methods (Access='public')
        function obj = LoH1Prolongation(feFunction)
            assert(isa(feFunction.fes.finiteElement, 'LowestOrderH1Fe'), ...
                'Prolongation only implemented for lowest order H1 conforming finite elements.')
            obj = obj@Prolongation(feFunction);
        end
    end
    
    methods (Access='protected')
        function interpolate(obj, src, ~)
            % use that for lowest order H1 elements the dofs correspond to
            % coordinates and new coordinates reside on edges or on inner nodes
            uInnerNodes = cell(numel(src.bisecGroups), 1);
            for k = 1:numel(src.bisecGroups)
                if ~isempty(src.bisecGroups{k}.innerNodes)
                    uInnerNodes{k} = eval(obj.u, src.bisecGroups{k}.innerNodes, ...
                        src.bisecGroups{k}.elementIdx);
                end
            end
            midpoint = Barycentric1D([1;1]/2);
            obj.interpolatedData = [obj.u.data, ...
                evalEdge(obj.u, midpoint, src.markedEdges), ...
                horzcat(uInnerNodes{:})];
        end
    end
end
