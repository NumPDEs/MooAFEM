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
        blf
        hoFes
        loFes
        lowfinest
        intergridMatrix
        D
        patchwiseChol
        patch
        changedPatches
        freeVertices
        freeVerticesOld
        P
    end
    
    %% event data
    properties (Access=private)
        listenerHandle
    end

    %% methods
    methods (Access=public)
        function obj = pLoc_MG(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end
            obj = obj@MGSolver();
            
            obj.hoFes = fes;
            obj.nLevels = 0;
            obj.D = {};
            obj.Acoarse = [];
            obj.patchwiseChol = {};
            obj.lowfinest = [];
            obj.intergridMatrix = cell(0);
            obj.changedPatches = cell(1);
            
            mesh = obj.hoFes.mesh;
            obj.loFes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', ':');
            obj.P = LoFeProlongation(obj.loFes);
            obj.blf = blf;
            % TODO: check coefficients of blf for compatibility with theoretical results!
            
            obj.listenerHandle = mesh.listener('IsAboutToRefine', @obj.getChangedPatches);
        end

        function setupLinearSystem(obj, A, b, x0)
            obj.nLevels = obj.nLevels + 1;            
            obj.lowfinest = assemble(obj.blf, obj.loFes);
            obj.freeVerticesOld = obj.freeVertices;
            obj.freeVertices = getFreeDofs(obj.loFes);
            obj.lowfinest = obj.lowfinest(obj.freeVertices, obj.freeVertices);
            if obj.nLevels >= 2
                obj.intergridMatrix{obj.nLevels} = ...
                    obj.P.matrix(obj.freeVertices, obj.freeVerticesOld);
                
                localVer = obj.changedPatches{obj.nLevels}; %numbering on all vertices
                [~,localVer] = ismember(localVer,obj.freeVertices);
                obj.changedPatches{obj.nLevels} = localVer(localVer>0); %numbering on all inner vertices
            end

            freeDofs = getFreeDofs(obj.hoFes);
            global2freeDofs = zeros(getDofs(obj.hoFes).nDofs, 1);
            global2freeDofs(freeDofs) = 1:numel(freeDofs);
            %if high order, assemble patch Dofs
            if obj.hoFes.finiteElement.order > 1
                obj.patch = assemblePatchDofs(obj.hoFes);

                % store cholesky decomposition of local matrices associated
                % to patches (only upper triangle)
                for ver = 1:numel(obj.patch)
                    idx = global2freeDofs(obj.patch{ver});
                    obj.patch{ver} = idx;
                    obj.patchwiseChol{ver} = chol(A(idx,idx));
                end
            end

            setupLinearSystem@MGSolver(obj, A, b, x0);
        end



        % Geometric MultiGrid
        function [Cx, algEta2] = Vcycle(obj, res)
            assert(isequal(size(res, 1), size(obj.A, 1)), ...
                'Setup for multilevel iteration not complete!')

            Ahigh = obj.A; %matrix on finest level, potentially high-order
            lev = obj.nLevels; %here numbering starts 1 for coarsest level
            poldeg = obj.hoFes.finiteElement.order;
            
            %built-in estimator of the algebraic error
            algEta2 = zeros(1, size(res, 2));

            %Vcycle only when more than one level;
            %otherwise coarse solve
            if lev > 1
                %p to 1 interpolation on the finest level of the residual
                if poldeg > 1
                    p1SmoothLev = lev-1; %until which level should the local smoothing be done
                    
                    resid{lev} = interpolateFreeData(res, obj.hoFes, obj.loFes);
                    
                    A{lev} = obj.lowfinest;
                else
                    p1SmoothLev = lev;
                    A{lev} = Ahigh;
                    resid{lev} = res;
                end

                % descending cascade in P1: no smoothing
                for k = lev:-1:2
                    A{k-1} = obj.intergridMatrix{k}'*A{k}*obj.intergridMatrix{k};
                    resid{k-1}  = obj.intergridMatrix{k}'*resid{k};
                end

                % exact solve on coarsest level
                % accumulative lifting of the residual: sigma
                sigma = A{1}\resid{1};
                algEta2 = algEta2 + sigma'*A{1}*sigma;

                % ascending cascade in P1 AND local smoothing
                for k = 2:p1SmoothLev
                    sigma = obj.intergridMatrix{k}*sigma;
                    uptres = resid{k} - A{k}*sigma; %updated residual 
                    
                    vershift = obj.changedPatches{k};
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
                    sigma = obj.intergridMatrix{lev}*sigma;
                    
                    sigma = interpolateFreeData(sigma, obj.loFes, obj.hoFes);

                    %updating the high-order residual
                    rho = zeros(size(sigma)); 
                    uptres = res - Ahigh*sigma;
                    for ver = 1: obj.hoFes.mesh.nCoordinates
                            
                        cholA = obj.patchwiseChol{ver};
                        dofshift = obj.patch{ver};

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
    
    methods (Access=protected)
        % callback function to be executed before mesh refinement:
        % -> patches change for all vertices adjacent to refined edges and
        %    all new vertices
        function getChangedPatches(obj, mesh, bisecData)
            nCOld = mesh.nCoordinates;
            nCNew = mesh.nCoordinates + nnz(bisecData.bisectedEdges) + ...
                bisecData.nRefinedElements'*bisecData.nInnerNodes;
            bisectedEdgeNodes = unique(mesh.edges(:,bisecData.bisectedEdges));
            obj.changedPatches{obj.nLevels+1} = ...
                [unique(bisectedEdgeNodes); ((nCOld+1):nCNew)'];
        end
    end
end

%% auxiliary functions
function interpolatedData = interpolateFreeData(data, fromFes, toFes)
    freeDofs = getFreeDofs(toFes);
    nComponents = size(data, 2);
    interpolatedData = zeros(numel(freeDofs), nComponents);
    feFunctionWrapper = FeFunction(fromFes);
    for k = 1:nComponents
        feFunctionWrapper.setFreeData(data(:,k));
        wholeData = nodalInterpolation(feFunctionWrapper, toFes);
        interpolatedData(:,k) = wholeData(freeDofs);
    end
end