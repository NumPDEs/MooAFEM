% PcgSolver (abstract subclass of IterativeSolver) Solves linear equations
%   iteratively by a preconditioned conjugate gradient (PCG) method.
%
% See also: IterativeSolver

classdef PcgSolver < IterativeSolver
    %% properties
    properties (SetAccess=protected, GetAccess=public)
        residualCNorm
    end
    
    properties (Access=protected)
        C
        residual
        Cresidual
        searchDirection
        normb
    end

    %% constructor
    methods (Access=public)
        function obj = PcgSolver(preconditioner)
            arguments
                preconditioner (1,1) Preconditioner
            end

            obj = obj@IterativeSolver();
            obj.C = preconditioner;
        end
    end
    
    %% extend superclass methods
    methods (Access=public)
        function setupSystemMatrix(obj, A)
            setupSystemMatrix@IterativeSolver(obj, A);
            obj.C.setup(A);
        end
        
        function setupRhs(obj, b, varargin)
            setupRhs@IterativeSolver(obj, b, varargin{:});
            
            % initialize residual & search direction
            obj.residual = b - obj.A * obj.x;
            obj.Cresidual = obj.C.apply(obj.residual);
            obj.searchDirection = obj.Cresidual;
            obj.residualCNorm = sqrt(dot(obj.residual, obj.Cresidual, 1));
            obj.normb = vecnorm(b, 2, 1);
        end
    end
    
    %% implement abstract superclass methods
    methods (Access=public)
        function tf = isConverged(obj)
            tf = ((obj.iterationCount >= obj.maxIter) ...
                | (obj.residualCNorm < obj.tol) ...
                | (obj.residualCNorm ./ obj.normb < obj.tol));
        end
    end
    
    methods (Access=protected)
        function computeUpdate(obj)
            % only iterate on active components
            idx = find(obj.activeComponents);
            
            % update solution
            AsearchDirection = obj.A * obj.searchDirection(:,idx);
            dAd = dot(obj.searchDirection(:,idx), AsearchDirection, 1);
            if dAd < eps
                alpha = 1;
            else
                alpha = obj.residualCNorm(:,idx).^2 ./ dAd;
            end
            obj.x(:,idx) = obj.x(:,idx) + alpha .* obj.searchDirection(:,idx);

            % update residual
            obj.residual(:,idx) = obj.residual(:,idx) - alpha .* AsearchDirection;
            obj.Cresidual(:,idx) = obj.C.apply(obj.residual(:,idx));
            residualCNormOld = obj.residualCNorm(:,idx);
            obj.residualCNorm(:,idx) = sqrt(dot(obj.residual(:,idx), obj.Cresidual(:,idx), 1));
            
            % update search direction
            beta = (obj.residualCNorm(:,idx) ./ residualCNormOld).^2;
            obj.searchDirection(:,idx) = obj.Cresidual(:,idx) + beta .* obj.searchDirection(:,idx);
        end
    end
end