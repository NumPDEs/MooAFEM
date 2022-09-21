classdef TestPatchwiseMatrix < matlab.unittest.TestCase
    
properties
    p
end

methods (TestMethodSetup)
    function initializeForms(testCase)
        % expected values below are computed by hand for unit triangle!
        testCase.p = 1:5;
    end
end

methods (Test)
    function lowestOrderIsDiagonalInverse(testCase)
        mesh = Mesh.loadFromGeometry('Lshape');
        fes = FeSpace(mesh, LowestOrderH1Fe(), 'dirichlet', [], 'neumann', ':');
        blf = BilinearForm();
        blf.a = Constant(mesh, 1);
        blf.c = Constant(mesh, 1);
        for ell = 1:5
            mesh.refineUniform();
            A = assemble(blf, fes);
            B = assemblePatchwise(blf, fes);
            res = full(diag(A)) .* (B \ ones(getDofs(fes).nDofs, 1));
            assert(all(abs(res - 1) < 10*eps))
        end
    end
end
    
end
