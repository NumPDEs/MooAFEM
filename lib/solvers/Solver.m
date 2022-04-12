% Solver (abstract handle class) Interface for solver class.

classdef Solver < handle
    methods (Abstract)
        solve(obj, A, b, varargin)
    end
    
    %% hooks to overwrite for additional steps
    methods
        function obj = Solver(varargin)
        end
        
        function setup(obj, varargin)
        end
    end
    
    %% validation functions to use within this class
    methods (Static, Access=protected)
        function mustBeQuadratic(A)
            if ~(size(A,1) == size(A,2))
                eidType = 'Solver:matrixNotQuadratic';
                msgType = 'Matrix must be quadratic.';
                throwAsCaller(MException(eidType,msgType))
            end
        end
        
        function mustBeCompatible(A, b)
            if ~(size(A,2) == size(b,1))
                eidType = 'Solver:dataNotCompatible';
                msgType = 'Vector must have size compatible with matrix.';
                throwAsCaller(MException(eidType,msgType))
            end
        end
    end
end