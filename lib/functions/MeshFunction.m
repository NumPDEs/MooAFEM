% MeshFunction (subclass of Evaluable) Encapsulate function handles within the
%   uniform interface of Evaluable to provide element/edge-wise  evaluation.
%
%   f = MeshFunction(mesh, funcHandle) encapsulates the function represented by
%       funcHandle to enable element/edge-wise evaluation on the given mesh.
%
%   See also: Evaluable

classdef MeshFunction < Evaluable
    %% properties
    properties (GetAccess='public', SetAccess='protected')
        mesh
        functionHandle
    end
    
    %% methods
    methods (Access='public')
        function obj = MeshFunction(mesh, functionHandle)
            arguments
                mesh (1,1) Mesh
                functionHandle (1,1) function_handle
            end
            
            obj.mesh = mesh;
            obj.functionHandle = functionHandle;
        end
        
        function val = eval(obj, bary, idx)
            arguments
                obj
                bary Barycentric2D
                idx {mustBeIndexVector} = ':'
            end
            
            val = obj.functionHandle(elementwiseCoordinates(obj.mesh, bary, idx));
        end
        
        function val = evalEdge(obj, bary, idx)
            arguments
                obj
                bary Barycentric1D
                idx {mustBeIndexVector} = ':'
            end
            
            val = obj.functionHandle(edgewiseCoordinates(obj.mesh, bary, idx));
        end
    end
end
