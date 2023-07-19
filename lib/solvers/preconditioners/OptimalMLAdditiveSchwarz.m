% OptimalMLAdditiveSchwarz (subclass of Preconditioner) optimal multilevel
%   additive Schwarz preconditioner with given smoother.
%
% See also: Preconditioner, PcgSolver, MultilevelSmoother

classdef OptimalMLAdditiveSchwarz < Preconditioner
    %% properties
    properties (Access=protected)
        smoother
        Acoarse
        projectedAcoarse
    end

    properties (Dependent)
        nLevels
    end

    properties (Access=private)
        listenerHandle
    end
    
    %% methods
    methods (Access=public)
        function obj = OptimalMLAdditiveSchwarz(smoother)
            arguments
                smoother MultilevelSmoother
            end
            
            obj = obj@Preconditioner();
            obj.smoother = smoother;
        end
        
        function setup(obj, A)
            if obj.nLevels == 0
                obj.Acoarse = A;
                obj.projectedAcoarse = obj.smoother.setup(A);
            else
                obj.smoother.setup(A);
            end
        end
           
        function Cx = apply(obj, x)
            assert(~isempty(obj.nLevels), 'Data for multilevel iteration not given!')

            L = obj.nLevels;
            rho = cell(L, 1);

            if L == 1
                Cx = obj.Acoarse \ x;
                return
            end

            % descending cascade
            residual = x;
            for k = L:-1:2
                rho{k} = obj.smoother.smooth(residual, k);
                residual = obj.smoother.restrict(residual, k);
            end
            
            % exact solve on coarsest level
            Cx = obj.projectedAcoarse \ residual;
            
            % ascending cascade
            for k = 2:L
                Cx = obj.smoother.prolongate(Cx, k);
                Cx = Cx + rho{k};
            end
        end
    end

    methods
        function n = get.nLevels(obj)
            n = obj.smoother.nLevels;
        end
    end
end