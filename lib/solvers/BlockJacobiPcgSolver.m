% BlockJacobiPcgSolver (subclass of PcgSolver) Solves linear equations
%   iteratively using the CG method with block-diagonal preconditioner.

classdef BlockJacobiPcgSolver < PcgSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        fes
    end

    properties (Access=private)
        C
        blf
    end
    
    %% methods
    methods (Access=public)
        % preconditioner action: patchwise inverse of diagonal
        function obj = BlockJacobiPcgSolver(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end
            
            obj = obj@PcgSolver();
            obj.fes = fes;
            obj.blf = blf;
        end

        function setupSystemMatrix(obj, A)
            obj.C = assemblePatchwise(obj.blf, obj.fes, ':');
            setupSystemMatrix@PcgSolver(obj, A);
        end
           
        function setupRhs(obj, b, varargin)
            setupRhs@PcgSolver(obj, b, varargin{:});
        end
        
        function Cx = preconditionAction(obj, x)
            Cx = obj.C \ x;
        end
    end
end