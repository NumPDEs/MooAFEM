% Evaluable (abstract handle class) Abstract base class for function classes.
%
%   Subclases must implement the eval method:
%       eval(obj, bary, idx)
%   where bary is a Barycentric2D and idx is an index vector. The method must
%   return the evaluation of the Evaluable on all barycentric coordinates in
%   bary on all triangles given by idx.

classdef Evaluable < handle
    %% properties
    properties (Abstract, GetAccess=public, SetAccess=protected)
        mesh
    end
    
    %% methods
    methods (Abstract, Access=public)
        eval(obj, bary, idx)
    end
    
    methods (Access=public)
        plot(obj, plotOptions)
        evalEdge(obj, bary, idx)
    end
end
