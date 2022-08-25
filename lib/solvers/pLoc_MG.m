% pLoc_MG (subclass of MGSolver) Solves linear equations
%   iteratively using high order p optimal multilevel
%   local (additive Schwarz/block Jacobi) smoothing
%   one iteration = Vcycle(0,1) : only 1 post-smoothing

classdef pLoc_MG < MGSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end

    properties (Access=protected)
        Acoarse
        feslow
        blflow
        lowfinest
        D
        locAhigh
        patch
        P
    end

    %% methods
    methods (Access=public)
        function obj = pLoc_MG(P)
            obj = obj@MGSolver();
            
            assert(isa(P, 'Prolongation'), 'P must be an instance of Prolongation')
            obj.P = P;
            obj.nLevels = 0;
            obj.D = {};
            obj.Acoarse = [];
            obj.locAhigh = {};
            obj.lowfinest = [];
            
            mesh = obj.P.fes.mesh;
            obj.feslow = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', ':');
            obj.blflow = BilinearForm(obj.feslow);
            obj.blflow.a = Constant(obj.feslow.mesh, 1);
            obj.blflow.qra = QuadratureRule.ofOrder(max(2*1-2, 1));
        end

        function setupLinearSystem(obj, A, b, x0)
            polDeg = obj.P.fes.finiteElement.order;
            obj.lowfinest = assemble(obj.blflow);

            %if high order, assemble patch Dofs
            if (polDeg>1)
                %adding the patch dof indices
                obj.patch = assemblePatchDofs(obj.P.fes);

                %storing cholesky decomposition of local matrices associated to patches
                for ver = 1:numel(obj.patch)
                    % U = chol(locmat); %stores only the upper triangular : locA = U'*U
                    obj.locAhigh{ver} = chol(A(obj.patch{ver},obj.patch{ver}));
                end
            end

            freeDofs = getFreeDofs(obj.P.fes);
            setupLinearSystem@MGSolver(obj, A(freeDofs,freeDofs), b, x0);
        end



        % Geometric MultiGrid
        function [Cx, algEta2] = Vcycle(obj, res)
            assert(isequal(getDofs(obj.feslow).nDofs, size(obj.lowfinest, 1)), ...
                'Data for multilevel iteration not given!')

            Ahigh = obj.A; %matrix on finest level, potentially high-order
            lev = numel(obj.P.fes.mesh.intergrid); %here numbering starts 1 for coarsest level
            poldeg = obj.P.fes.finiteElement.order;
            intergrid = obj.P.fes.mesh.intergrid;
            
            %built-in estimator of the algebraic error
            algEta2 = zeros(1, size(res, 2));

            %Vcycle only when more than one level;
            %otherwise coarse solve
            if lev > 1
                freeDofslow = getFreeDofs(obj.feslow);
                %p to 1 interpolation on the finest level of the residual
                if poldeg > 1
                    p1SmoothLev = lev-1; %until which level should the local smoothing be done
                    
                    rup = FeFunction(obj.P.fes);
                    rup.setFreeData(res(:,1));
                    res1temp = nodalInterpolation(rup, obj.feslow);
                    resid{lev} = res1temp(freeDofslow);
                    
                    A{lev} = obj.lowfinest(freeDofslow,freeDofslow);
                else
                    p1SmoothLev = lev;
                    A{lev} = Ahigh;
                    resid{lev} = res;
                end

                % descending cascade in P1: no smoothing
                for k = lev:-1:2
                    I = intergrid{k};
                    newfree = obj.P.fes.mesh.freeVert{k};
                    oldfree = obj.P.fes.mesh.freeVert{k-1};
                    Ifree = I(newfree,oldfree);
                    A{k-1} = Ifree'*A{k}*Ifree;
                    resid{k-1}  = Ifree'*resid{k};
                end

                % exact solve on coarsest level
                % accumulative lifting of the residual: sigma
                sigma = A{1}\resid{1};
                algEta2 = algEta2 + sigma'*A{1}*sigma;

                % ascending cascade in P1 AND local smoothing
                for k = 2:p1SmoothLev
                    I = intergrid{k};
                    newfree = obj.P.fes.mesh.freeVert{k};
                    oldfree = obj.P.fes.mesh.freeVert{k-1};
                    Ifree = I(newfree,oldfree);

                    sigma = Ifree*sigma;
                    uptres = resid{k} - A{k}*sigma; %updated residual 

                    localVer = obj.P.fes.mesh.locVert{k}; %numbering on all vertices
                    [~,vershift] = ismember(localVer,newfree);
                    vershift = vershift(vershift>0); %numbering on all inner vertices

                    obj.D = diag(A{k});
                    rho = zeros(size(sigma)); 
                    rho(vershift,:) = obj.D(vershift).^(-1).*uptres(vershift,:);

                    %error correction with optimal stepsize 
                    lambda = (uptres'*rho)/(rho'*A{k}*rho);

                    if (lambda > 3)
                       warning('MG step-sizes no longer bound by d+1!')
                    end

                    sigma = sigma + lambda*rho;
                    algEta2 = algEta2 + (rho'*A{k}*rho)*(lambda^2);
                end

                %finest level high-order patch-problems ONLY for p > 1
                if poldeg > 1
                    %back to finest level: still P1 at this point
                    I = intergrid{lev};
                    newfree = obj.P.fes.mesh.freeVert{lev};
                    oldfree = obj.P.fes.mesh.freeVert{lev-1};
                    Ifree = I(newfree,oldfree);
                    sigma = Ifree*sigma;

                    %interpolation of the error correction from 1 to p
                    freeDofshigh = getFreeDofs(obj.P.fes);
                    
                    sigma1func = FeFunction(obj.feslow);
                    sigma1func.setFreeData(sigma(:,1));
                    sigma1p = nodalInterpolation(sigma1func, obj.P.fes);
                    sigma = sigma1p(freeDofshigh);

                    %updating the high-order residual
                    rho = zeros(size(sigma)); 
                    uptres = res - Ahigh*sigma;
                    for ver = 1: obj.P.fes.mesh.nCoordinates
                            
                        cholA = obj.locAhigh{ver};
                        dofspatch = obj.patch{ver};
                        [~,dofshift] = ismember(dofspatch,freeDofshigh);
                        dofshift = dofshift(dofshift>0); %renumbering on freedofs vector

                        %solve local patch problems
                        temp = (cholA')\(uptres(dofshift,:));
                        rhoa = cholA\temp;
                        rho(dofshift,:) = rho(dofshift,:) + rhoa;%additive
                    end

                    %optimal stepsize
                    lambda = (uptres'*rho)/(rho'*Ahigh*rho);

                    if (lambda > 3)
                       warning('MG step-sizes no longer bound by d+1!')
                    end

                    %error correction
                     sigma = sigma + lambda*rho;
                     algEta2 = algEta2 + (rho'*Ahigh*rho)*(lambda^2);
                end

                Cx = sigma;

            else %one coarse level setting only: exact solve
                sigma = Ahigh \ res;
                Cx = sigma;
                algEta2 = sigma'*Ahigh*sigma;
            end
        end
    end
end