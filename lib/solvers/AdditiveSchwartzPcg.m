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
        D
    end
    
    %% methods
    methods (Access=public)
        function obj = AdditiveSchwartzPcg()
            obj = obj@PcgSolver();
            obj.nLevels = 0;
            obj.fine2coarse = [];
            obj.D = {};
            obj.Acoarse = [];
            obj.transferMatrix = [];
        end
        
        function setup(obj, A, b, x0, P)
            assert(isa(P, 'Prolongation'), 'P must be an instance of Prolongation')
            setup@PcgSolver(obj, A, b, x0);
            
            mesh = P.fes.mesh;
            
            
            
            if isempty(obj.nLevels)
                % on coarsest level store whole matrix
                obj.nLevels = 0;
                obj.Acoarse = obj.A;
            else
                % on every other level only store diagonal
                obj.nLevels = obj.nLevels + 1;
                
                % Find nodes of newly created nodes + neigbours where patch has changed
                % this is to ensure that no node gets considered twice, since fine2coarse
                % also contains node-relations in areas where no refinement happens
                % TODO: is there a better way to get this information?
                [~,oldNodes] = find(obj.fine2coarse == 1);
                patchAreaOldProlongated = zeros(mesh.nCoords(),1);
                patchAreaOldProlongated(oldNodes) = obj.patchAreaOld;
                nonRefinedArea = abs(patchAreaOldProlongated-obj.patchArea) < 2*eps;
                obj.transferMatrix{obj.nLevels} = obj.fine2coarse(obj.freeDofsOld, obj.freeDofs);
                
                % invert diagonal (only on free dofs and dofs where refinement happened)
                obj.D{obj.nLevels} = 1./full(diag(A));
                obj.D{obj.nLevels}(nonRefinedArea) = 0;
                obj.D{obj.nLevels} = obj.D{obj.nLevels}(obj.freeDofs);
            end
        end
        
        % add transfer between old and new coordinates based on marked elements
        function registerRefinedElements(obj, marked)
            % TODO: get rid of refineMeshFB, there is a lot of unnecessary stuff going on
            % obtain full mapping coarse mesh to fine mesh
            mesh = obj.fes.mesh;
            
            % refineMeshFB cannot deal with empty neumann boundary
            if isempty(mesh.neumann)
                [~,~,~,~,~,obj.fine2coarse] = refineMeshFB(mesh.coordinates, mesh.elements, mesh.dirichlet, marked, [], 1);
            elseif isempty(mesh.dirichlet)
                [~,~,~,~,~,obj.fine2coarse] = refineMeshFB(mesh.coordinates, mesh.elements, mesh.neumann, marked, [], 1);
            else
                [~,~,~,~,~,~,obj.fine2coarse] = refineMeshFB(mesh.coordinates, mesh.elements, mesh.dirichlet, mesh.neumann, marked, [], [], 1);
            end
        end
        
        % preconditioner: inverse of diagonal on each level
        function Cx = preconditionAction(obj, x)
            assert(~isempty(obj.nLevels), 'Data for multilevel iteration not given!')
            y = cell(obj.nLevels+1, 1);
            
            % descending cascade
            for k = obj.nLevels:-1:1
                y{k+1} = obj.D{k} .* x;
                x = obj.transferMatrix{k}*x;
            end
            
            % exact solve on coarsest level
            y{1} = obj.Acoarse \ x;
            
            % ascending cascade
            for k = 1:obj.nLevels
                y{k+1} = y{k+1} + obj.transferMatrix{k}'*y{k};
            end
            
            Cx = y{end};
        end
    end
end