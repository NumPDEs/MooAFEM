% AdditiveSchwartzPcg (subclass of PcgSolver) Solves linear equations
%   iteratively using the CG method with optimal multilevel additive
%   Schwartz preconditioner.

classdef AdditiveSchwartzPcg < PcgSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end
    
    properties (Access=protected)
        Acoarse
        fine2coarse
        transferMatrix
        oldFreeDofs
        oldPatchAreas
        D
        P
    end
    
    %% methods
    methods (Access=public)
        function obj = AdditiveSchwartzPcg(P)
            assert(isa(P, 'Prolongation'), 'P must be an instance of Prolongation')
            assert(isa(P.fes.finiteElement, 'LowestOrderH1Fe'), 'Only lowest-order H1 finite elements are supported')
            
            obj = obj@PcgSolver();
            obj.nLevels = 0;
            obj.fine2coarse = [];
            obj.D = {};
            obj.Acoarse = [];
            obj.transferMatrix = [];
            obj.P = P;
        end
        
        function setupLinearSystem(obj, A, b, x0)
            mesh = obj.P.fes.mesh;
            freeDofs = getFreeDofs(obj.P.fes);
            patchAreas = computePatchAreas(mesh);
            
            if obj.nLevels == 0
                % on coarsest level store whole matrix
                obj.Acoarse = A;
            else    
                % find nodes where patch has changed to ensure that no node gets
                % considered twice (new nodes are just appended)
                obj.oldPatchAreas(mesh.nCoordinates) = 0;
                stableNodes = abs(obj.oldPatchAreas - patchAreas) < 2*eps;
                obj.transferMatrix{obj.nLevels} = obj.P.matrix(freeDofs, obj.oldFreeDofs);
                
                % on every other level only store inverse of diagonal
                % (only on free dofs and dofs where refinement happened;
                % we heavily use that dofs and nodes coincide, here)
                obj.D{obj.nLevels} = 1./full(diag(A));
                obj.D{obj.nLevels}(stableNodes(freeDofs)) = 0;
            end
            
            obj.oldFreeDofs = freeDofs;
            obj.oldPatchAreas = patchAreas;
            obj.nLevels = obj.nLevels + 1;
            
            setupLinearSystem@PcgSolver(obj, A, b, x0);
        end
        
        % preconditioner: inverse of diagonal on each level
        function Cx = preconditionAction(obj, x)
            assert(~isempty(obj.nLevels), 'Data for multilevel iteration not given!')
            y = cell(obj.nLevels, 1);
            
            % descending cascade
            for k = obj.nLevels-1:-1:1
                y{k+1} = obj.D{k} .* x;
                x = obj.transferMatrix{k}'*x;
            end
            
            % exact solve on coarsest level
            y{1} = obj.Acoarse \ x;
            
            % ascending cascade
            for k = 1:obj.nLevels-1
                y{k+1} = y{k+1} + obj.transferMatrix{k}*y{k};
            end
            
            Cx = y{end};
        end
    end
end