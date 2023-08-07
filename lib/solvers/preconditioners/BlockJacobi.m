% BlockJacobi (subclass of Preconditioner) block-diagonal preconditioner.
%
% See also: Preconditioner, PcgSolver

classdef BlockJacobi < Preconditioner
    %% properties
    properties (Access=private)
        C
        fes
        blf
    end
    
    %% methods
    methods (Access=public)
        % preconditioner action: patchwise inverse of diagonal
        function obj = BlockJacobi(fes, blf)
            arguments
                fes FeSpace
                blf BilinearForm
            end
            
            obj.fes = fes;
            obj.blf = blf;
        end

        function setup(obj, A)
            obj.C = assemblePatchwise(obj.blf, obj.fes, ':');
        end
           
        function Cx = apply(obj, x)
            Cx = obj.C \ x;
        end
    end
end