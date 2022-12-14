% chooseIterativeSolver Return suitable solver (instance of subclass of
%   IterativeSolver) and suitable Prolongation P.

function [solver, P] = chooseIterativeSolver(fes, blf, class, variant)
arguments
    fes FeSpace
    blf BilinearForm
    class (1,1) string {mustBeMember(class, ["cg", "pcg", "multigrid"])}
    variant (1,1) string {mustBeMember(variant, [ "", ...
        "jacobi", "additiveSchwarz", ...
        "lowOrderVcycle", "highOrderVcycle"])} = ""
end

order = fes.finiteElement.order;
P = Prolongation.chooseFor(fes);

switch class
    % non-preconditioned CG
    case "cg"
        solver = CgSolver();
        
    % preconditioned CG family
    case "pcg"
        switch variant
            case "jacobi"
                if order == 1, solver = JacobiPcgSolver(fes);
                else, solver = BlockJacobiPcgSolver(fes, blf); end
            case {"", "additiveSchwarz"}
                solver = AdditiveSchwartzPcg(fes, blf);
            otherwise
                error('No PCG variant %s!', variant)
        end
        
    
    % geometric multigrid family
    case "multigrid"
        switch variant
            case {"", "lowOrderVcycle"}
                if order == 1, solver = LowestOrderLocalMg(fes, blf);
                else, solver = LocalMgLowOrderVcycle(fes, blf); end
            case "highOrderVcycle"
                if order == 1, solver = LowestOrderLocalMg(fes, blf);
                else, solver = LocalMgHighOrderVcycle(fes, blf, P); end
            otherwise
                error('No multigrid variant %s!', variant)
        end
   
    % default
    otherwise
        error('No iterative solver of class %s!', class)
end
    
end


