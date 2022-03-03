% AffineTransformation (handle class) Encapsulates data about the affine
%   transformation from a reference triangle to the mesh.
% 
%   The reference triangle is given as {(0,0), (1,0), (0,1)}
%   Here, F: \hat{T} -> T is the affine transformation and DF its Jacobian.
%   This is used in the transformation formula for integrals. For
%   transforming derivatives of integrands, the (transposed) inverse Jacobian
%   DF^(-T) = DFinv is used for \nabla u = DF^(-T) * \nabla \hat{u}.

classdef AffineTransformation < handle
    %% properties
    properties (GetAccess = 'public', SetAccess = 'private')
        DFinv (4,:) double
        area (1,:) double
        ds (1,:) double
        n (2,:) double
    end
    
    %% methods
    methods
        function obj = AffineTransformation(mesh)
            % compute edge Vectors
            z1 = mesh.coordinates(:,mesh.elements(1,:));
            edgeVectors21 = mesh.coordinates(:,mesh.elements(2,:)) - z1;
            edgeVectors31 = mesh.coordinates(:,mesh.elements(3,:)) - z1;
            
            % compute element areas (= 1/2 * element-wise determinant of edge vectors)
            obj.area = (edgeVectors21(1,:).*edgeVectors31(2,:) - edgeVectors21(2,:).*edgeVectors31(1,:)) / 2;
            
            % compute DFinv from edge vectors (column-wise)
            obj.DFinv = [edgeVectors31(2,:); -edgeVectors21(2,:); ...
                -edgeVectors31(1,:);  edgeVectors21(1,:)] ./ (2*obj.area);
            
            % compute edge transformations and normals
            edgeVectors = mesh.coordinates(:,mesh.edges(2,:)) - mesh.coordinates(:,mesh.edges(1,:));
            % ds (= edge length)
            obj.ds = sqrt(vectorProduct(edgeVectors, edgeVectors));
            % unit normal vector (flipped right 90 degrees relative to edge)
            obj.n = [edgeVectors(2,:); -edgeVectors(1,:)] ./ obj.ds;
        end
    end
end
