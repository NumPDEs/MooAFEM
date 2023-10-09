% BilinearForm (handle class) Collects data to assemble the FEM matrix
%   A_ij = int_Omega Dphi_i'*a*Dphi_j + b*Dphi_j*phi_i + c*phi_j*phi_i dx
%          + int_{Gamma_R} robin*phi_j*phi_i ds.

classdef BilinearForm < handle    
    %% properties
    properties (Access=public)
        a {mustBeEvaluableOrEmpty} = []     % Diffusion matrix (or scalar)
        b {mustBeEvaluableOrEmpty} = []     % Convection vector field
        c {mustBeEvaluableOrEmpty} = []     % Reaction coefficient
        robin {mustBeEvaluableOrEmpty} = [] % Robin coefficient
        qra (1,1) QuadratureRule = UnspecifiedQR
        qrb (1,1) QuadratureRule = UnspecifiedQR
        qrc (1,1) QuadratureRule = UnspecifiedQR
        qrRobin (1,1) QuadratureRule = UnspecifiedQR
    end
    
    %% methods
    methods (Access=public)
        mat = assemble(obj, fes);
        mat = assemblePatchwise(obj, fes);
    end
    
    methods (Access=protected)
        data = computeVolumeData(obj, fes);
    end
    
    methods (Static, Access=protected)
        val = diffusionPart(a, Dphi);
        val = convectionPart(b, Dphi, phi);
        val = reactionPart(c, phi);
        val = robinPart(robin, phi);
    end
end
