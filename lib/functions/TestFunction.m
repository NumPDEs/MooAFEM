% TestFunction (subclass of Evaluable) Vector of test functions on a finite
%   element space.
%
%   f = TestFunction(fes) initializes a test function on the finite element
%       space fes. The evaluation of this function is a vector of the values of
%       all shape functions at the given barycentric coordinate.
%
%   See also: Evaluable

classdef TestFunction < Evaluable
    %% properties
    properties (SetAccess=protected, GetAccess=public)
        fes
        mesh
    end
    
    %% methods
    methods (Access=public)
        function obj = TestFunction(fes)
            arguments
                fes (1,1) FeSpace
            end
            
            obj.fes = fes;
            obj.mesh = fes.mesh;
        end
        
        function val = eval(obj, bary, ~)
            arguments
                obj
                bary Barycentric2D
                ~
            end
            
            val = evalShapeFunctions(obj.fes.finiteElement, bary);
        end
        
        function val = evalEdge(obj, bary, ~)
            arguments
                obj
                bary Barycentric1D
                ~
            end
            
            val = evalEdgeShapeFunctions(obj.fes.finiteElement, bary);
        end
    end
end
