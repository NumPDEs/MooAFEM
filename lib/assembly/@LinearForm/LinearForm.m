% LinearForm (handle class) Collects data to assemble the FEM vector
%   F_i = int_Omega f phi_i + fvec*Dphi_i dx + int_{Gamma_N} neumann*phi_i ds
%          + int_{Gamma_R} robin*phi_i ds.

classdef LinearForm < handle
    %% properties
    properties (Access=public)
        f {mustBeEvaluableOrEmpty} = []         % Volume load
        fvec {mustBeEvaluableOrEmpty} = []      % Generalized volume load
        robin {mustBeEvaluableOrEmpty} = []     % Robin boundary load
        neumann {mustBeEvaluableOrEmpty} = []   % Neumann boundary load
        qrf (1,1) QuadratureRule = UnspecifiedQR
        qrfvec (1,1) QuadratureRule = UnspecifiedQR
        qrRobin (1,1) QuadratureRule = UnspecifiedQR
        qrNeumann (1,1) QuadratureRule = UnspecifiedQR
    end
    
    %% methods
    methods (Access=public)
        vec = assemble(obj, fes);
    end
    
    methods (Static, Access=protected)
        val = scalarPart(f, phi);
        val = vectorPart(fvec, Dphi);
        val = neumannPart(neumann, phi);
        val = robinPart(robin, phi);
    end
end
