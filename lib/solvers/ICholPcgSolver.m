% ICholPcgSolver (subclass of PcgSolver) Solves linear equations iteratively
%   using the CG method with incomplete Cholesky decomposition preconditioner.

classdef ICholPcgSolver < PcgSolver
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
        function obj = ICholPcgSolver(fes)
            arguments
                fes FeSpace
            end
            
            obj = obj@PcgSolver();
            obj.fes = fes;
        end

        function setupSystemMatrix(obj, A)
            obj.C = ichol(A);
            setupSystemMatrix@PcgSolver(obj, A);
        end
           
        function setupRhs(obj, b, varargin)
            setupRhs@PcgSolver(obj, b, varargin{:});
        end
        
        function Cx = preconditionAction(obj, x)
            Cx = obj.C' \ (obj.C \ x);
        end
    end
end