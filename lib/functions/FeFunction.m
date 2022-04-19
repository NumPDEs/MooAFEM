% FeFunction (subclass of Evaluable) Function on a finite element space. The
%   data of the function is flagged invalid if the underlying mesh changes.
%
%   f = FeFunction(fes) initializes a function on the finite element space fes.
%
%   f.setData(val) sets the data to the constant val.
%
%   f.setData(data) explicitely sets the data of the finite element function.
%       Knowledge of the underlying data structure is required.
%
%   f.setFreeData(val) sets the data to the constant val on all dofs that are
%   flagged as 'free' in the finite element space.
%
%   f.setData(data) explicitely sets the data on all dofs that are flagged as
%       'free' in the finite element space. Knowledge of the underlying data
%       structure is required.
%
%   See also: Evaluable

classdef FeFunction < Evaluable
    %% properties
    properties (SetAccess=protected, GetAccess=public)
        mesh
        fes
        data (1,:) double
        dataIsValid (1,1) logical
    end
    
    properties (Access = private)
        listenerHandle
    end
    
    %% methods
    methods (Access=public)
        function obj = FeFunction(fes)
            arguments
                fes (1,1) FeSpace
            end
            
            obj.fes = fes;
            obj.mesh = fes.mesh;
            obj.data = [];
            obj.dataIsValid = true;
            
            % register FeFunction as an observer of mesh
            obj.listenerHandle = fes.mesh.listener('JustRefined', @obj.handleChangedMesh);
        end
        
        function setData(obj, data)
            dofs = getDofs(obj.fes);
            assertAdmissibleDataSize(obj, data, dofs.nDofs)
            obj.allocateZeroDataIfWrongFormat();
            obj.data(:) = data(:);
            obj.dataIsValid = true;
        end
        
        function setFreeData(obj, data)
            freeDofs = getFreeDofs(obj.fes);
            assertAdmissibleDataSize(obj, data, length(freeDofs))
            obj.allocateZeroDataIfWrongFormat();
            obj.data(freeDofs) = data(:);
            obj.dataIsValid = true;
        end
        
        function val = eval(obj, bary, idx)
            arguments
                obj
                bary Barycentric2D
                idx {mustBeIndexVector}= ':'
            end
            
            assertDataValidity(obj);
            dofs = getDofs(obj.fes);
            phi = evalShapeFunctions(obj.fes.finiteElement, bary);
            I = dofs.element2Dofs(:,idx);
            val = vectorProduct(reshape(obj.data(I), size(I)), phi);
        end
        
        function val = evalEdge(obj, bary, idx)
            arguments
                obj
                bary Barycentric1D
                idx {mustBeIndexVector}= ':'
            end
            
            assertDataValidity(obj);
            dofs = getDofs(obj.fes);
            phi = evalEdgeShapeFunctions(obj.fes.finiteElement, bary);
            I = dofs.edge2Dofs(:,idx);
            val = vectorProduct(reshape(obj.data(I), size(I)), phi);
        end
    end
    
    methods (Access = private)
        function handleChangedMesh(obj, ~, ~)
            obj.dataIsValid = false;
        end
        
        function allocateZeroDataIfWrongFormat(obj)
            if isempty(obj.data) || ~obj.dataIsValid
                dofs = getDofs(obj.fes);
                obj.data = zeros(dofs.nDofs, 1);
            end
        end
        
        function assertAdmissibleDataSize(~, data, size)
            if ~(isscalar(data) || (numel(data) == size))
                me = MException('FeFunction:wrongDataDimension', ...
                    'Data must correspond to free dofs of underlying FeSpace.');
                throwAsCaller(me);
            end
        end
        
        function assertDataValidity(obj)
            if ~obj.dataIsValid
                me = MException('FeFunction:dataNotValid', ...
                    'Data cannot be evaluated, since underlying Mesh / FeSpace has changed.');
                throwAsCaller(me);
            end
        end
    end
end
