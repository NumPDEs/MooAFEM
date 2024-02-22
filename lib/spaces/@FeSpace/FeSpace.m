% FeSpace (handle class) Contains a finite element and handles the global
%   connectivity for degrees of freedom (DOFs). Frequently used data are cached
%   for efficient use.

classdef FeSpace < handle
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        mesh
        finiteElement
        bnd
    end
    
    properties (Access=private)
        dofs
        freeDofs
        fixedDofs
    end
    
    %% events
    properties (Access=private)
        listenerHandle
    end
    
    %% methods
    methods (Access=public)
        function obj = FeSpace(mesh, fe, bnd)
            % Construct FeSpace from given mesh and finite element.
            %
            %   fes = FeSpace(mesh, fe) constructs FeSpace from FiniteElement fe
            %       on mesh, whole boundary is dirichlet boundary.
            %
            %   fes = FeSpace(mesh, fe, bndName, bndIdx, ...) Boundaries can be
            %       given by tuples consisting of the name ('dirichlet',
            %       'neumann', 'robin') and an index vector.
            
            arguments
                mesh (1,1) Mesh
                fe (1,1) FiniteElement
                bnd.dirichlet {mustBeIndexVector} = []
                bnd.neumann {mustBeIndexVector} = []
                bnd.robin {mustBeIndexVector} = []
            end
            
            % register fespace as an observer of mesh
            obj.listenerHandle = mesh.listener('JustRefined', @obj.handleChangedMesh);
            
            % initialize lazily evaluated properties
            obj.dofs = [];
            obj.freeDofs = [];
            
            obj.mesh = mesh;
            obj.finiteElement = fe;
            
            % store indices of boundaries and do sanity checks
            obj.bnd = bnd;
            if all(structfun(@isempty, bnd))
                obj.bnd.dirichlet = ':';
            end
            
            bndNames = {'dirichlet', 'neumann', 'robin'};
            referenceCount = zeros(1, mesh.nBoundaries);
            for k=1:numel(bndNames)
                idx = obj.bnd.(bndNames{k});
                if ~ischar(idx)
                    mustBeInRange(idx, 1, mesh.nBoundaries)
                end
                referenceCount(idx) = referenceCount(idx) + 1;
            end
            if any(referenceCount ~= 1)
                error('FeSpace:invalidBoundary', ...
                    'Every boundary part must be connected to exactly one boundary.')
            end
        end
        
        function dofs = getDofs(obj)
            if isempty(obj.dofs)
                obj.dofs = obj.assembleDofs();
            end
            dofs = obj.dofs;
        end
        
        function freeDofs = getFreeDofs(obj)
            if isempty(obj.freeDofs)
                obj.freeDofs = obj.assembleFreeDofs();
            end
            freeDofs = obj.freeDofs;
        end
        
        function fixedDofs = getFixedDofs(obj)
            if isempty(obj.fixedDofs)
                obj.fixedDofs = ...
                    setdiff(1:obj.getDofs().nDofs, obj.getFreeDofs(), 'stable');
            end
            fixedDofs = obj.fixedDofs;
        end
        
        dofs = createVertexDofs(obj, idx)
        dofs = createInnerEdgeDofs(obj, idx)
        dofs = createInnerElementDofs(obj, idx)
    end
    
    methods (Access=protected)
        dofs = assembleDofs(obj)
        freeDofs = assembleFreeDofs(obj)
    end
    
    methods (Access=private)
        handleChangedMesh(obj, src, event)
    end
end
