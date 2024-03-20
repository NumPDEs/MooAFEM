% choose Return suitable solver (instance of subclass of IterativeSolver) and
%   suitable Prolongation P.
%
%   solver = choose(fes, blf, class) returns a solver of class 'class' for given
%       finite element space fes and bilinear form blf. For each class, sensible
%       default variants are chosen.
%
%   solver = choose(fes, blf, class, variant) returns a solver of class 'class'
%       of specific variant 'variant'.
%
%   [solver, P] = choose(__) additionally to the solver, returns a suitable
%       Prolongation P, that can be used to prolongate the solution of the
%       solver to the next finer mesh.

function [solver, P] = choose(fes, blf, method, variant)
arguments
    fes FeSpace
    blf BilinearForm
    method (1,1) string {mustBeMember(method, ["cg", "pcg", "multigrid", "direct"])}
    variant (1,1) string {mustBeMember(variant, [ "", ...
        "iChol", "jacobi", ...
        "additiveSchwarzLowOrder", "additiveSchwarzHighOrder" ...
        "lowOrderVcycle", "highOrderVcycle", ...
        "multiplicativeCR", "vcycleCR", "additiveCR"])} = ""
end

order = fes.finiteElement.order;
P = Prolongation.chooseFor(fes);

switch method
    % non-preconditioned CG
    case "cg"
        solver = PcgSolver(NoPreconditioner());
        
    % preconditioned CG family
    case "pcg"
        switch variant
            case "iChol"
                preconditioner = IncompleteCholesky();
            case "jacobi"
                if order == 1
                    preconditioner = DiagonalJacobi();
                else
                    preconditioner = BlockJacobi(fes, blf);
                end
            case {"", "additiveSchwarzLowOrder"}
                if order == 1
                    preconditioner = OptimalMLAdditiveSchwarz(P1JacobiCascade(fes, blf, P));
                else
                    preconditioner = OptimalMLAdditiveSchwarz(JacobiLowOrderCascade(fes, blf));
                end
            case "additiveSchwarzHighOrder"
                if order == 1
                    preconditioner = OptimalMLAdditiveSchwarz(P1JacobiCascade(fes, blf, P));
                else
                    preconditioner = OptimalMLAdditiveSchwarz(JacobiHighOrderCascade(fes, blf, P));
                end
            case "additiveCR"
                Av = LoMeshAveraging(fes);  % TODO: use proper choose function
                smoother = CRJacobiCascade(fes, blf, Av, P);
                preconditioner = LocalAdditiveMultilevelCRPreconditioner(smoother);
            otherwise
                error('No PCG variant %s!', variant)
        end
        solver = PcgSolver(preconditioner);
        
    % geometric , Pmultigrid family
    case "multigrid"
        switch variant
            case {"", "lowOrderVcycle"}
                if order == 1
                    smoother = P1JacobiCascade(fes, blf, P);
                else
                    smoother = JacobiLowOrderCascade(fes, blf);
                end
                solver = OptimalVcycleMultigridSolver(smoother);
            case "highOrderVcycle"
                if order == 1
                    smoother = P1JacobiCascade(fes, blf, P);
                else
                    smoother = JacobiHighOrderCascade(fes, blf, P);
                end
                solver = OptimalVcycleMultigridSolver(smoother);
            case "multiplicativeCR"
                Av = LoMeshAveraging(fes);  % TODO: use proper choose function
                smoother = CRJacobiCascade(fes, blf, Av, P);
                solver = LocalMultiplicativeMultigridCRSolver(smoother);
            case "vcycleCR"
                Av = LoMeshAveraging(fes);  % TODO: use proper choose function
                smoother = CRJacobiCascade(fes, blf, Av, P);
                solver = LocalMultigridVcycleCRSolver(smoother);
            otherwise
                error('No multigrid variant %s!', variant)
        end
        
    % direct solve (for testing purposes)
    case "direct"
        solver = DirectSolver();
   
    % default
    otherwise
        error('No iterative solver of class %s!', method)
end
    
end
