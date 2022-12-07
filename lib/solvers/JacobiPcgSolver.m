% JacobiPcgSolver (subclass of PcgSolver) linear equations iteratively using
% the CG method with diagonal preconditioner.

classdef JacobiPcgSolver < PcgSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        C
        hoFes
        blf

        highestOrderIsOne
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

            obj.highestOrderIsOne = (fes.finiteElement.order == 1);

            obj.hoFes = fes;
            obj.blf = blf; % TODO: check coefficients of blf for compatibility with theoretical results!
            
            %obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end
        function setupSystemMatrix(obj, A)
            
            if obj.highestOrderIsOne
                obj.C = full(diag(A)).^(-1);
            else
                obj.C = assemblePatchwise(obj.blf, obj.hoFes, ':');
            end
            
            setupSystemMatrix@PcgSolver(obj, A);
        end
           
        function setupRhs(obj, b, x0)
            setupRhs@PcgSolver(obj, b, x0);
        end
        
        function Cx = preconditionAction(obj, x)
            Cx = obj.C \ x;
        end
    end
end