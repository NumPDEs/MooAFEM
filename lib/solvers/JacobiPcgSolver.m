% JacobiPcgSolver (subclass of PcgSolver) linear equations iteratively using
% the CG method with diagonal preconditioner.

classdef JacobiPcgSolver < PcgSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        C
        hoFes
        blf
    end

    properties (Access=private)
        setupPreconditioner
        apply
    end
    
    %% methods
    methods (Access=public)
        % preconditioner action: inverse of diagonal
        function obj = JacobiPcgSolver(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end
            obj = obj@PcgSolver();
            obj.hoFes = fes;
            obj.blf = blf;

            if fes.finiteElement.order == 1
                obj.setupPreconditioner = @(A, blf, hoFes) full(diag(A)).^(-1);
                obj.apply = @(C, x) C .* x;
            else
                obj.setupPreconditioner = @(A, blf, hoFes) assemblePatchwise(blf, hoFes, ':');
                obj.apply = @(C, x) C \ x;
            end
        end

        function setupSystemMatrix(obj, A)
            obj.C = obj.setupPreconditioner(A, obj.blf, obj.hoFes);
            setupSystemMatrix@PcgSolver(obj, A);
        end
           
        function setupRhs(obj, b, varargin)
            setupRhs@PcgSolver(obj, b, varargin{:});
        end
        
        function Cx = preconditionAction(obj, x)
            Cx = obj.apply(obj.C, x);
        end
    end
end