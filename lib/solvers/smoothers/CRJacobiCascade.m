% CRJacobiCascade (subclass of MultilevelSmoother) multilevel Jacobi smoother
%   for lowest order nonconforming Crouzeix-Raviart finite elements: smooth 
%   with diagonal of stiffness matrix on changed patches on every level.
%
% See also: MultilevelSmoother, OptimalVcycleMultigridSolver


classdef CRJacobiCascade < MultilevelSmoother
    properties (Access=protected)
        Av
        P
        confFes
        patchwiseA
        intergridAveraging  % averaging of CR on fine level
        intergridProlongation  % inclusion of conforming S1 given by CR coefficients
        freeDofs
        freeDofsOld
    end
    
    %% event data
    properties (Access=private)
        listenerHandle
    end

    %% methods
    methods (Access=public)
        function obj = CRJacobiCascade(fes, blf, Av, P)
            arguments
                fes
                blf
                Av Prolongation
                P Prolongation
            end

            assert(isa(fes.finiteElement, 'LowestOrderCRFe'), ...
                ['CRJacobiCascade only works for lowest order nonconforming', ... 
                'Crouzeix-Raviart finite elements.'])
            assert(fes.finiteElement.order == 1, ...
                'CRJacobiCascade only works for lowest order finite elements.')
            
            obj = obj@MultilevelSmoother(fes, blf);
            obj.Av = Av;
            obj.P = P;
            obj.confFes = FeSpace(fes.mesh, LowestOrderH1Fe);
        end

        function nonInvertedSmootherMatrix = setup(obj, A)
            setup@MultilevelSmoother(obj, A);
            
            L = obj.nLevels;
            obj.freeDofsOld = obj.freeDofs;
            obj.freeDofs = getFreeDofs(obj.fes);
            freeVertices = getFreeDofs(obj.confFes);
            nonInvertedSmootherMatrix = A;
            
            if L > 1
                obj.intergridAveraging{L} = obj.Av.matrix(obj.freeDofs, obj.freeDofsOld);
                obj.intergridProlongation{L} = obj.P.matrix(obj.freeDofs, obj.freeDofsOld);
                obj.changedPatches{L} = find(obj.changedPatches{L}(freeVertices));
                obj.patchwiseA{L} = assemblePatchwise(obj.blf, obj.fes, obj.changedPatches{L});
            end
        end

        function Cx = smooth(obj, x, k, nSteps)
            if nargin < 4
                nSteps = 1;
            end
            Cx = obj.patchwiseA{k} \ x;
            for j = 2:nSteps
                Cx = obj.patchwiseA{k} \ Cx;
            end
        end

        function Px = prolongate(obj, x, k)
            Px = obj.intergridProlongation{k} * x;
        end

        function Px = restrict(obj, x, k)
            Px = obj.intergridProlongation{k}' * x;
        end

        function Px = average(obj, x, k)
            Px = obj.intergridAveraging{k} * x;
        end

        function Px = averageRestrict(obj, x, k)
            Px = obj.intergridAveraging{k}' * x;
        end

        % Debugging
        function plot(obj, k)
            plot(obj.fes.mesh)
        end
    end
end
