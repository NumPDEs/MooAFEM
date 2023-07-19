% chooseIterativeSolver Return suitable solver (instance of subclass of
%   IterativeSolver) and suitable Prolongation P.

function [solver, P] = chooseIterativeSolver(fes, blf, class, variant)
arguments
    fes FeSpace
    blf BilinearForm
    class (1,1) string {mustBeMember(class, ["cg", "pcg", "multigrid", "direct"])}
    variant (1,1) string {mustBeMember(variant, [ "", ...
        "iChol", "jacobi", ...
        "additiveSchwarzLowOrder", "additiveSchwarzHighOrder" ...
        "lowOrderVcycle", "highOrderVcycle"])} = ""
end

order = fes.finiteElement.order;
P = Prolongation.chooseFor(fes);

switch class
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
            otherwise
                error('No PCG variant %s!', variant)
        end
        solver = PcgSolver(preconditioner);
        
    % geometric multigrid family
    case "multigrid"
        switch variant
            case {"", "lowOrderVcycle"}
                if order == 1
                    smoother = P1JacobiCascade(fes, blf, P);
                else
                    smoother = JacobiLowOrderCascade(fes, blf);
                end
            case "highOrderVcycle"
                if order == 1
                    smoother = P1JacobiCascade(fes, blf, P);
                else
                    smoother = JacobiHighOrderCascade(fes, blf, P);
                end
            otherwise
                error('No multigrid variant %s!', variant)
        end
        solver = OptimalVcycleMultigridSolver(smoother);
        
    % direct solve (for testing purposes)
    case "direct"
        solver = DirectSolver();
   
    % default
    otherwise
        error('No iterative solver of class %s!', class)
end
    
end


