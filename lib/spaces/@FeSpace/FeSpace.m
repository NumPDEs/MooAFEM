% FeSpace (handle class) Contains a finite element and handles the global
%   connectivity for degrees of freedom (DOFs). Frequently used data are cached
%   for efficient use.

classdef FeSpace < handle
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        mesh
        finiteElement
        dirichlet
    end
    
    properties (Access=private)
        dofs
        freeDofs
    end
    
    %% events
    properties (Access=private)
        listenerHandle
    end
    
    %% methods
    methods (Access=public)
        function obj = FeSpace(mesh, fe, opt)
            % Construct FeSpace from given mesh and finite element.
            %
            %   fes = FeSpace(mesh, fe) constructs FeSpace from FiniteElement fe
            %       on mesh, whole boundary is dirichlet boundary.
            %
            %   fes = FeSpace(mesh, fe, 'dirichlet', dir) Dirichlet boundary is
            %       given by dir, which is either an array of boundary indices,
            %       or 'all'.
            
            arguments
                mesh (1,1) Mesh
                fe (1,1) FiniteElement
                opt.dirichlet = 'all'
            end
            
            % register fespace as an observer of mesh
            obj.listenerHandle = mesh.listener('JustRefined', @obj.handleChangedMesh);
            
            % initialize lazily evaluated properties
            obj.dofs = [];
            obj.freeDofs = [];
            
            obj.mesh = mesh;
            obj.finiteElement = fe;
            
            % store indices of dirichlet boundary
            if isa(opt.dirichlet, 'char')
                if strcmp(opt.dirichlet, 'all')
                    obj.dirichlet = 1:mesh.nBoundaries;
                else
                    error('FeSpace:invalidDirichletBoundary', 'Choice of dirichlet boundary not valid!')
                end
            elseif isa(opt.dirichlet, 'double')
                if any(opt.dirichlet < 1 | opt.dirichlet > mesh.nBoundaries, 'all')
                    error('FeSpace:invalidDirichletBoundary', 'Choice of dirichlet boundary not valid!')
                else
                    obj.dirichlet = unique(opt.dirichlet);
                end
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
    end
    
    methods (Access=protected)
        dofs = assembleDofs(obj)
        freeDofs = assembleFreeDofs(obj)
    end
    
    methods (Access=private)
        handleChangedMesh(obj, src, event)
    end
end
