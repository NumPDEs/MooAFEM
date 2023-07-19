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
                    preconditioner = OptimalMLAdditiveSchwarz(P1JacobiSmoother(fes, blf, P));
                else
                    preconditioner = OptimalMLAdditiveSchwarz(JacobiLowOrderCascadeSmoother(fes, blf));
                end
            case "additiveSchwarzHighOrder"
                if order == 1
                    preconditioner = OptimalMLAdditiveSchwarz(P1JacobiSmoother(fes, blf, P));
                else
                    preconditioner = OptimalMLAdditiveSchwarz(JacobiHighOrderCascadeSmoother(fes, blf, P));
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
                    smoother = P1JacobiSmoother(fes, blf, P);
                else
                    smoother = JacobiLowOrderCascadeSmoother(fes, blf);
                end
            case "highOrderVcycle"
                if order == 1
                    smoother = P1JacobiSmoother(fes, blf, P);
                else
                    smoother = JacobiHighOrderCascadeSmoother(fes, blf, P);
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


