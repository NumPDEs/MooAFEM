% LocalAdditiveMultilevelCRPreconditioner 

classdef LocalAdditiveMultilevelCRPreconditioner < Preconditioner
    %% properties
    properties (Access=protected)
        smoother
        Acoarse
    end

    properties (Dependent)
        nLevels
    end

    %% methods
    methods
        function obj = LocalAdditiveMultilevelCRPreconditioner(smoother)
            obj = obj@Preconditioner();
            obj.smoother = smoother;
        end
        
        function setup(obj, A)
            if obj.nLevels == 0
                obj.Acoarse = A;
            end
            obj.smoother.setup(A);
        end
        
        function Cx = apply(obj, x)
            Cx = Vcycle(obj, x);
        end

        function n = get.nLevels(obj)
            n = obj.smoother.nLevels;
        end

    end
        
    methods (Access=protected)
        function Cx = Vcycle(obj, x)
            % one step of the Vcycle
            
            % if there is only one coarse level: exact solve
            L = obj.nLevels;
            if L == 1
                Cx = obj.Acoarse \ x;
                % algEstimator = sqrt(dot(Cx, obj.A * Cx, Dim.Vector));
                return
            end

            res = x;

            % descending cascade: no smoothing
            resPerLevel = cell(L, 1);
            resPerLevel{L} = res;
            for k = L:-1:3
                resPerLevel{k-1} = obj.smoother.averageRestrict(res, k);
                res = obj.smoother.restrict(res, k);
            end
            resPerLevel{1} = obj.smoother.averageRestrict(res, 2);

            % exact solution step on level 1
            Cx = obj.Acoarse \ resPerLevel{1};
            Cx = obj.smoother.average(Cx, 2);

            % ascending cascade: no smoothing
            for k = 3:L
                update = obj.smoother.smooth(resPerLevel{k-1}, k-1);
                Cx = obj.smoother.prolongate(Cx, k) + obj.smoother.average(update, k);
            end

            % final coarse smooth
            Cx = obj.smoother.smooth(Cx, L);
        end
    end
end