% Prolongation (abstract handle class) Interface for prolongation of FeFunctions
%   from coarse to fine mesh.
%
%   prolongate(P, u) returns the prolongated data of FeFunction u.
%
%   restrict(P, u) returns the restricted data of FeFunction u.

classdef Prolongation < handle
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        matrix {mustBeSparse}
    end
    
    
    %% methods
    methods (Access=public)
        function obj = Prolongation(fes) %#ok<INUSD>
            obj.matrix = [];
        end
        
        function newFeData = prolongate(obj, u)
            arguments
                obj
                u FeFunction
            end
            
            if isempty(obj.matrix)
                error('Prolongation matrix not set. Are you sure there was mesh refinement?')
            end
            newFeData = obj.matrix * u.data';
        end

        function newFeData = restrict(obj, u)
            arguments
                obj
                u FeFunction
            end
            
            if isempty(obj.matrix)
                error('Prolongation matrix not set. Are you sure there was mesh refinement?')
            end
            newFeData = (u.data * obj.matrix)';
        end
    end
    
    methods (Static, Access=public)
        function P = chooseFor(fes)
            if isa(fes.finiteElement, 'LowestOrderH1Fe') ...
                    || isa(fes.finiteElement, 'LowestOrderL2Fe')
                P = LoMeshProlongation(fes);
            else
                P = MeshProlongation(fes);
            end
        end
    end
end

function mustBeSparse(A)
    if ~(isempty(A) || issparse(A))
        eidType = 'Prolongation:matrixNotSparse';
        msgType = 'Property ''matrix'' must be a sparse matrix.';
        throwAsCaller(MException(eidType,msgType))
    end
end
