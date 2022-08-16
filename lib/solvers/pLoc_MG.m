% pLoc_MG (subclass of MGSolver) Solves linear equations
%   iteratively using high order p optimal multilevel
%   local (additive Schwarz/block Jacobi) smoothing
%   one iteration = Vcycle(0,1) : only 1 post-smoothing

classdef pLoc_MG < MGSolver
    %% properties
    properties (GetAccess=public,SetAccess=protected)
        nLevels
    end

    properties (Access=public)
        Acoarse
        feslow
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
        end

        function setupLinearSystem(obj, A, b, x0)
            mesh = obj.P.fes.mesh;
            polDeg = obj.P.fes.finiteElement.order;

            %TO BE IMPROVED : later 1->p on the finest level in matrix form
            obj.feslow = FeSpace(obj.P.fes.mesh, HigherOrderH1Fe(1), 'dirichlet', 'all');
            blflow = BilinearForm(obj.feslow);
            blflow.a = Constant(mesh, 1);
            blflow.qra = QuadratureRule.ofOrder(max(2*1-2, 1));
            obj.lowfinest = assemble(blflow);

            %if high order, assemble patch Dofs
            if (polDeg>1)

                %adding the patch dof indices
                dofs = assemblePatchDofs(obj.P.fes,polDeg);

                %storing cholesky decomposition of local matrices associated to patches
                for ver = 1:mesh.nCoordinates
                    % U = chol(locmat); %stores only the upper triangular : locA = U'*U
                    obj.locAhigh{ver} = chol(A(dofs.patch2Dofs{ver},dofs.patch2Dofs{ver}));
                end

                obj.patch = dofs.patch2Dofs;
            end

            freeDofs = getFreeDofs(obj.P.fes);
            setupLinearSystem@MGSolver(obj, A(freeDofs,freeDofs), b, x0);
        end



        % Geometric MultiGrid
        function [Cx, algEta2] = Vcycle(obj, res)

            assert(~isempty(obj.nLevels), 'Data for multilevel iteration not given!')


            Ahigh = obj.A; %matrix on finest level, potentially high-order
           % res = obj.b - Ahigh*obj.x;
            lev = numel(obj.P.fes.mesh.intergrid); %here numbering starts 1 for coarsest level
            poldeg = obj.P.fes.finiteElement.order;
            intergrid = obj.P.fes.mesh.intergrid;

            dual = (size(res(1,:)) == 2);
            
            %built-in estimator of the algebraic error
            algEta2 = zeros(size(res(1,:)));

            %Vcycle only when more than one level;
            %otherwise coarse solve
            if lev > 1

                freeDofslow = getFreeDofs(obj.feslow);
                %p to 1 interpolation on the finest level of the residual
                if poldeg > 1
                    p1SmoothLev = lev-1; %until which level should the local smoothing be done
                    if length(res(1,:))==2 %only primal or prima/dual solve
    
                        rup = FeFunction(obj.P.fes);
                        rup.setFreeData(res(:,1));
                        res1temp = nodalInterpolation(rup, obj.feslow);
    
                        rzp = FeFunction(obj.P.fes);
                        rzp.setFreeData(res(:,2));
                        res2temp = nodalInterpolation(rzp, obj.feslow);
                       
                        resid{lev} = [res1temp(freeDofslow), res2temp(freeDofslow)];
                    else 
                        rup = FeFunction(obj.P.fes);
                        rup.setFreeData(res(:,1));
                        res1temp = nodalInterpolation(rup, obj.feslow);
                        resid{lev} = res1temp(freeDofslow); 

                    end
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
                
                algEta2(1) = algEta2(1) + sigma(:,1)'*A{1}*sigma(:,1);
                if dual
                    algEta2(2) = algEta2(2) + sigma(:,2)'*A{1}*sigma(:,2);
                end

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

                    D = diag(A{k});
                    rho = zeros(size(sigma)); 
                    rho(vershift,:) = D(vershift).^(-1).*uptres(vershift,:);


                    %error correction with optimal stepsize 
                    lambda = (uptres(:,1)'*rho(:,1))/(rho(:,1)'*A{k}*rho(:,1));
                    if dual
                       lambda = [lambda, (uptres(:,2)'*rho(:,2))/(rho(:,2)'*A{k}*rho(:,2))];
                    end

                    if (lambda > 3)
                       warning('MG step-sizes no longer bound by d+1!')
                    end

                    sigma(:,1) = sigma(:,1) + lambda(1)*rho(:,1);
                    algEta2(1) = algEta2(1) + (rho(:,1)'*A{k}*rho(:,1))*(lambda(1)^2);
                    if dual
                       sigma(:,2) = sigma(:,2) + lambda(2)*rho(:,2);
                       algEta2(2) = algEta2(2) + (rho(:,2)'*A{k}*rho(:,2))*(lambda(2)^2);
                    end
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
                     if length(res(1,:))==2 %only primal or prima/dual solve
                         sigma1func = FeFunction(obj.feslow);
                         sigma1func.setFreeData(sigma(:,1));
                         sigma1p = nodalInterpolation(sigma1func, obj.P.fes);

                         sigma2func = FeFunction(obj.feslow);
                         sigma2func.setFreeData(sigma(:,2));
                         sigma2p = nodalInterpolation(sigma2func, obj.P.fes);
                         sigma = [sigma1p(freeDofshigh), sigma2p(freeDofshigh)];
                     else
                         sigma1func = FeFunction(obj.feslow);
                         sigma1func.setFreeData(sigma(:,1));
                         sigma1p = nodalInterpolation(sigma1func, obj.P.fes);
                         sigma = sigma1p(freeDofshigh);
                     end

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
                    lambda = (uptres(:,1)'*rho(:,1))/(rho(:,1)'*Ahigh*rho(:,1));
                    if dual
                       lambda = [lambda, (uptres(:,2)'*rho(:,2))/(rho(:,2)'*Ahigh*rho(:,2))];
                    end

                    if (lambda > 3)
                       warning('MG step-sizes no longer bound by d+1!')
                    end

                    %error correction
                     sigma(:,1) = sigma(:,1) + lambda(1)*rho(:,1);
                     algEta2(1) = algEta2(1) + (rho(:,1)'*Ahigh*rho(:,1))*(lambda(1)^2);
                    if dual
                       sigma(:,2) = sigma(:,2) + lambda(2)*rho(:,2);
                       algEta2(2) = algEta2(2) + (rho(:,2)'*Ahigh*rho(:,2))*(lambda(2)^2);
                    end
                end



                Cx = sigma;

            else %one coarse level setting only: exact solve
                

                sigma = Ahigh \ res;
                Cx = sigma;
                algEta2(1) = sigma(:,1)'*Ahigh*sigma(:,1);
                if dual
                    algEta2(2) = sigma(:,2)'*Ahigh*sigma(:,2);
                end

            end


        end
    end
end