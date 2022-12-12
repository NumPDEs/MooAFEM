% chooseIterativeSolver Return suitable subclass of IterativeSolver.

function solver = chooseIterativeSolver(fes, blf, class, variant)
arguments
    fes FeSpace
    blf BilinearForm
    class (1,1) string {mustBeMember(class, ["cg", "pcg", "multigrid"])}
    variant (1,1) string {mustBeMember(variant, [ "", ...
        "jacobi", "additiveSchwarz", ...
        "lowOrderVcycle", "highOrderVcycle"])} = ""
end

switch class
    % non-preconditioned CG
    case "cg"
        solver = CgSolver();
        
    % preconditioned CG family
    case "pcg"
        switch variant
            case "jacobi"
                solver = JacobiPcgSolver(fes, blf);
            case {"", "additiveSchwarz"}
                solver = AdditiveSchwartzPcg(fes, blf);
            otherwise
                error('No PCG variant %s!', variant)
        end
        
    
    % geometric multigrid family
    case "multigrid"
        switch variant
            case {"", "lowOrderVcycle"}
                solver = OptimalLocalMGSolver(fes, blf);
            case "highOrderVcycle"
                solver = OptimalLocalMGSolver1pp(fes, blf);
            otherwise
                error('No multigrid variant %s!', variant)
        end
   
    % default
    otherwise
        error('No iterative solver of class %s!', class)
end
    
end


