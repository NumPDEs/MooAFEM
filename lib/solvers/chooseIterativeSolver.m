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
        solver = PcgSolverTEMP(IdentityPC());
        
    % preconditioned CG family
    case "pcg"
        switch variant
            case "iChol"
                solver = ICholPcgSolver(fes);
            case "jacobi"
                if order == 1, solver = JacobiPcgSolver(fes);
                else, solver = BlockJacobiPcgSolver(fes, blf); end
            case {"", "additiveSchwarzLowOrder"}
                if order == 1, solver = LowestOrderAdditiveSchwarzPcg(fes, blf, P);
                else
                    solver = AdditiveSchwarzLowOrderPcg(fes, blf);
                end
            case "additiveSchwarzHighOrder"
                solver = AdditiveSchwarzHighOrderPcg(fes, blf, P);
            otherwise
                error('No PCG variant %s!', variant)
        end
        
    % geometric multigrid family
    case "multigrid"
        switch variant
            case {"", "lowOrderVcycle"}
                if order == 1, solver = LowestOrderLocalMg(fes, blf, P);
                else, solver = LocalMgLowOrderVcycle(fes, blf); end
            case "highOrderVcycle"
                if order == 1, solver = LowestOrderLocalMg(fes, blf, P);
                else, solver = LocalMgHighOrderVcycle(fes, blf, P); end
            otherwise
                error('No multigrid variant %s!', variant)
        end
        
    % direct solve (for testing purposes)
    case "direct"
        solver = DirectSolver();
   
    % default
    otherwise
        error('No iterative solver of class %s!', class)
end
    
end


