% CompositeFunction (subclass of Evaluable) Combine multiple Evaluables to a new function.
%
%   f = CompositeFunction(funcHandle, funcArg1, ..., funcArgN) combine
%       Evaluables funcArg1, ..., funcArgN by means of the function handle
%       funcHandle. funcHandle must have exactly N input arguments, representing
%       each of the Evaluables. In this way, complicated functions can be
%       created, e.g., for coefficients of (bi-)linear forms. In particular,
%       also FeFunctions on different FeSpaces can be combined.
%       For edge evaluation, funcArg1, ..., funcArgN must implement the evalEdge
%       method.
%
%   See also: Evaluable

classdef CompositeFunction < Evaluable
    %% properties
    properties (SetAccess='protected', GetAccess='public')
        mesh
        functionHandle
        functionArguments
    end
    
    %% methods
    methods (Access='public')
        function obj = CompositeFunction(functionHandle, functionArguments)
            arguments
                functionHandle (1,1) function_handle
            end
            
            arguments (Repeating)
                functionArguments (1,1) Evaluable
            end
            
            obj.functionHandle = functionHandle;
            obj.functionArguments = functionArguments;
            obj.mesh = functionArguments{1}.mesh;
        end
        
        function val = eval(obj, varargin)
            val = compositeEval(obj, @eval, varargin{:});
        end
        
        function val = evalEdge(obj, varargin)
            val = compositeEval(obj, @evalEdge, varargin{:});
        end
    end
    
    methods (Access='private')
        function val = compositeEval(obj, evalMethod, bary, idx)
            arguments
                obj
                evalMethod
                bary Barycentric
                idx {mustBeIndexVector} = ':'
            end
            
            evaluatedArgs = cellfun(@(x) evalMethod(x, bary, idx), ...
                obj.functionArguments, 'UniformOutput', false);
            val = obj.functionHandle(evaluatedArgs{:});
        end
    end
end
