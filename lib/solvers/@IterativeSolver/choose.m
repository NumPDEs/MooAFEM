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

function [solver, P] = choose(fes, blf, class, variant)
arguments
    fes FeSpace
    blf BilinearForm
    class (1,1) string {mustBeMember(class, ["cg", "pcg", "multigrid", "direct"])}
    variant (1,1) string {mustBeMember(variant, [ "", ...
        "iChol", "jacobi", ...
        "additiveSchwarzLowOrder", "additiveSchwarzHighOrder" ...
        "lowOrderVcycle", "highOrderVcycle", "CRjacobi"])} = ""
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
        
    % geometric , Pmultigrid family
    case "multigrid"
        switch variant
            case {"", "lowOrderVcycle"}
                if order == 1
                    switch class(fes)
                        case "LowestOrderH1Fe"
                            smoother = P1JacobiCascade(fes, blf, P);
                        case "LowestOrderCRFe"
                            smoother = CRJacobiCascade(fes, blf, P);
                    otherwise
                        error('No multigrid for finite element %s!', class(fes))
                    end
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
