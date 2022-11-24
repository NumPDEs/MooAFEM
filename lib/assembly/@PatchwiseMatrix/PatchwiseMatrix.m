% PatchwiseMatrix (handle class) Patchwise matrix of bilinear form.

classdef PatchwiseMatrix < handle
    properties (GetAccess=public, SetAccess=protected)
        fes
        patchwiseChol
        L
        V
    end
    
    properties (Access=protected)
        patchDofs
        patchElements
        global2freeDofs
    end

    properties (Dependent)
        nPatches
    end

    methods (Access=public)
        function obj = PatchwiseMatrix(fes, data, L, V)
            obj.fes = fes;
            [obj.patchDofs, obj.patchElements] = assemblePatchDofs(obj);
            obj.buildLocalMatrixFactorization(data);
            obj.L = L;
            obj.V = V;
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
            nPatches = numel(obj.patchwiseChol);
        end
    end

    methods (Access=protected)
        [patchDofs, patchElements] = assemblePatchDofs(obj)
        buildLocalMatrixFactorization(obj, data)
    end
end

