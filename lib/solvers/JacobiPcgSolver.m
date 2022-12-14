% JacobiPcgSolver (subclass of PcgSolver) Solves linear equations
%   iteratively using the CG method with diagonal preconditioner.

classdef JacobiPcgSolver < PcgSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        fes
    end

    properties (Access=private)
        C
    end
    
    %% methods
    methods (Access=public)
        % preconditioner action: inverse of diagonal
        function obj = JacobiPcgSolver(fes)
            arguments
                fes FeSpace
            end
            
            obj = obj@PcgSolver();
            obj.fes = fes;
        end

        function setupSystemMatrix(obj, A)
            obj.C = full(diag(A)).^(-1);
            setupSystemMatrix@PcgSolver(obj, A);
        end
           
        function setupRhs(obj, b, varargin)
            setupRhs@PcgSolver(obj, b, varargin{:});
        end
        
        function Cx = preconditionAction(obj, x)
            Cx = obj.C .* x;
        end
    end
end