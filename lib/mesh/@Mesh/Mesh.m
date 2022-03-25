% Mesh (handle class) Encapsulates geometric data necessary for a mesh:
%   coordinates, edges, elements, and boundary parts.
% 
%   Transformation information is computed on access and stored in an internal
%   field until the mesh changes. All data fields are read only.

classdef Mesh < handle
    %% properties
    properties (GetAccess = 'public', SetAccess = 'private')
        coordinates (2,:) double
        elements (3,:) double
        edges (2,:) double
        element2edges (3,:) double
        flipEdges (3,:) logical
        boundaries (:,1) cell
    end
    
    properties (GetAccess = 'public', SetAccess = 'private', Hidden = true)
        markedEdges (:,1) logical
        bisecGroups (:,1) cell
    end
    
    %% dependent properties (computed from data)
    properties (Dependent)
        nCoordinates
        nElements
        nEdges
        nBoundaries
    end
    
    %% private properties for caching purposes
    properties (Access = 'private')
        trafo
    end
    
    %% events
    events
        IsAboutToRefine
        HasChanged
    end
    
    %% public methods
    methods (Access = 'public')
        function obj = Mesh(coordinates, elements, boundaries)
            % Construct Mesh object from given coordinate, element, and boundary
            % arrays.
            %
            %   Mesh(coordinates, elements, boundaries)
            
            obj.coordinates = coordinates;
            obj.elements = elements;
            [obj.edges, obj.element2edges, obj.flipEdges, obj.boundaries] ...
                = obj.computeEdgeInformation(obj.elements, boundaries);
            obj.trafo = [];
        end
        
        % get methods for cached data
        function trafo = getAffineTransformation(obj)
            if isempty(obj.trafo)
                obj.trafo = AffineTransformation(obj);
            end
            trafo = obj.trafo;
        end
        
        % further methods
        points = elementwiseCoordinates(obj, bary, idx)
        points = edgewiseCoordinates(obj, bary, idx)
        patchAreas = computePatchAreas(obj)
        plot(obj, varargin)
        saveTikzConforming(obj, varargin)
        refineLocally(obj, marked, method)
        refineUniform(obj, n, method)
        changeRefinementEdge(obj, newRefinementEdge)
        edges = getCombinedBndEdges(obj, idx)
    end
    
    %% protected methods
    methods (Access='protected')
        updateData(obj)
    end
    
    %% get methods for dependent data
    methods
        function val = get.nCoordinates(obj)
            val = size(obj.coordinates, Dim.Elements);
        end
        
        function val = get.nElements(obj)
            val = size(obj.elements, 2);
        end
        
        function val = get.nEdges(obj)
            val = size(obj.edges, 2);
        end
        
        function val = get.nBoundaries(obj)
            val = numel(obj.boundaries);
        end
    end
    
    %% static methods
    methods (Static)
        obj = loadFromGeometry(geometryIdentifier)
        [edges, element2edges, flipEdges, boundary2edges] = computeEdgeInformation(elements, boundaries)
    end
end
