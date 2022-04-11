% Prolongation (abstract handle class) Interface for prolongation of FeFunctions
%   from coarse to fine mesh.

classdef Prolongation < handle
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        matrix {mustBeSparse}
    end
    
    properties (Access='protected')
        fes FeSpace
        listenerHandle
    end
    
    %% methods
    methods (Access='public')
        function obj = Prolongation(fes)
            obj.fes = fes;
            mesh = fes.mesh;
            obj.matrix = [];
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.setupMatrix);
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
    end
    
    methods (Abstract, Access='protected')
        setupMatrix(obj, src, event)
    end
end

function mustBeSparse(A)
    if ~(isempty(A) || issparse(A))
        eidType = 'Prolongation:matrixNotSparse';
        msgType = 'Property ''matrix'' must be a sparse matrix.';
        throwAsCaller(MException(eidType,msgType))
    end
end