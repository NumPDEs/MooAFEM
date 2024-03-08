% LocalAdditiveMultilevelCRSmoother 

classdef LocalAdditiveMultilevelCRSmoother < CRJacobiCascade
    %% methods
    methods
        function obj = LocalAdditiveMultilevelCRSmoother(fes, blf, Av, P)
            obj = obj@CRJacobiCascade(fes, blf, Av, P);
        end
    end
        
    methods (Access=protected)

        % Take from super-class
        % prolongate(obj, x, k)
        % restrict(obj, x, k)
        
        function Cx = smooth(obj, x, k)
        end

        function Cx = cycle(obj, x)
            assert(isequal(size(x, 1), size(obj.A, 1)), ... 
                'Setup for multilevel iteration not complete!')
            
            % if there is only one coarse level: exact solve
            L = obj.nLevels;
            if L == 1
                Cx = obj.A \ x;
                algEstimator = sqrt(dot(Cx, obj.A * Cx, Dim.Vector));
                return
            end

            % sequential multiplicative V-cycle update
            algEstimator = 0;
            res = x;
            v = obj.x;
            for ell = 1:L
                % Apply V-cycle down to level ell
                xUpdate = obj.Vcycle(res, ell);
                % A priori step size
                if ell == L
                    rho = 1;
                else
                    rho = 0.8;
                end
                v = v + rho * xUpdate;
                % Compute optimal step size
                % updatedResidual = x - obj.A * Cx;
                % rho = obj.smoother.smooth(updatedResidual, ell+1);
                % [etaUpdate2, xUpdate] = MultigridSolver.computeOptimalUpdate(obj.A, updatedResidual, rho);
                % Cx = Cx + xUpdate;
                % Update estimator
                % algEstimator = algEstimator + etaUpdate2;
                % Update residual
                res = res - obj.A * rho * xUpdate;
            end
            Cx = v;
        end

        function Cx = Vcycle(obj, x, ell) 
            % one step of the Vcycle
            
            L = obj.nLevels;

            if ell == L
                Cx = obj.smoother.smooth(x, ell);
                return
            end

            % descending cascade: no smoothing
            res = x;
            for k = L:-1:ell+2
                res = obj.smoother.restrict(res, k);
            end
            res = obj.smoother.averageRestrict(res, ell+1);

            % smoothing step on level ell
            if ell == 1
                Cx = obj.projectedMatrix{1} \ res;
            else
                Cx = obj.smoother.smooth(res, ell);
            end
            Cx = obj.smoother.average(Cx, ell+1);

            % ascending cascade: no smoothing
            for k = ell+2:L
                Cx = obj.smoother.prolongate(Cx, k);
            end
        end
    end
end