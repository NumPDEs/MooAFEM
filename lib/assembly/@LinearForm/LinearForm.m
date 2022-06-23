% LinearForm (handle class) Collects data to assemble the FEM vector
%   F_i = int_Omega f phi_i + fvec*Dphi_i dx + int_{Gamma_N} neumann*phi_i ds
%          + int_{Gamma_R} robin*phi_i ds.

classdef LinearForm < handle
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        fes
    end
    
    properties (Access=public)
        f {mustBeEvaluableOrEmpty} = []         % Volume load
        fvec {mustBeEvaluableOrEmpty} = []      % Generalized volume load
        robin {mustBeEvaluableOrEmpty} = []     % Robin boundary load
        neumann {mustBeEvaluableOrEmpty} = []   % Neumann boundary load
        qrf (1,1) QuadratureRule = QuadratureRule.ofOrder(1)
        qrfvec (1,1) QuadratureRule = QuadratureRule.ofOrder(1)
        qrRobin (1,1) QuadratureRule = QuadratureRule.ofOrder(1, '1D')
        qrNeumann (1,1) QuadratureRule = QuadratureRule.ofOrder(1, '1D')
        bndRobin {mustBeIndexVector} = []
        bndNeumann {mustBeIndexVector} = []
    end
    
    %% methods
    methods
        function obj = LinearForm(fes)
        % LinearForm Construct bilinear form from FeSpace.
            arguments
                fes (1,1) FeSpace
            end
            
            obj.fes = fes;
        end
    end
    
    methods (Access=public)
        vec = assemble(obj);
    end
    
    methods (Static, Access=protected)
        val = scalarPart(f, phi);
        val = vectorPart(fvec, Dphi);
        val = neumannPart(neumann, phi);
        val = robinPart(robin, phi);
    end
end
