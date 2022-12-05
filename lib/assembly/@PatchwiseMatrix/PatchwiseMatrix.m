% PatchwiseMatrix (handle class) Patchwise matrix of bilinear form.

classdef PatchwiseMatrix < handle
    properties (GetAccess=public, SetAccess=protected)
        fes
        activePatches
    end
    
    properties (Access=protected)
        patchDofs
        patchElements
        patchesAsFreeDofs
        patchwiseChol
    end

    properties (Dependent)
        nPatches
    end

    methods (Access=public)
        function obj = PatchwiseMatrix(fes, data, patches)
            arguments
                fes FeSpace
                data
                patches {mustBeIndexVector}
            end
            
            obj.fes = fes;
            [obj.patchDofs, obj.patchElements] = assemblePatchDofs(obj);
            idx = 1:obj.nPatches;
            obj.activePatches = reshape(idx(patches), 1, []);
            obj.buildLocalMatrixFactorization(data);
        end
        
        function dofs = getPatchDofs(obj, idx)
            arguments
                obj
                idx {mustBeIndexVector}
            end
            dofs = vertcat(obj.patchDofs{idx});
        end
        
        function elements = getPatchElements(obj, idx)
            arguments
                obj
                idx {mustBeIndexVector}
            end
            elements = vertcat(obj.patchElements{idx});
        end
    end

    methods
        function nPatches = get.nPatches(obj)
            nPatches = numel(obj.patchDofs);
        end
    end

    methods (Access=protected)
        [patchDofs, patchElements] = assemblePatchDofs(obj)
        buildLocalMatrixFactorization(obj, data)
    end
end

