% computeVolumeData Computes raw data of bilinear form for volume terms.
%
%   data = computeData(obj, fes) computes the elementwise data of the given
%       BilinearForm for a FeSpace fes for every non-trivial dof-pairing.
%
%   See also: BilinearForm

function data = computeVolumeData(obj, fes)
arguments
    obj
    fes FeSpace
end

if all([isempty(obj.a), isempty(obj.b), isempty(obj.c), isempty(obj.robin)])
    error('BilinearForm:NoCoefficients', 'No coefficients are given.')
end

%% integrate element data
data = 0;
phi = TestFunction(fes);
Dphi = TestFunctionGradient(fes);

% diffusion
if ~isempty(obj.a)
    f = CompositeFunction(@BilinearForm.diffusionPart, obj.a, Dphi);
    data = data + integrateElement(f, obj.qra);
end

% convection
if ~isempty(obj.b)
    f = CompositeFunction(@BilinearForm.convectionPart, obj.b, Dphi, phi);
    data = data + integrateElement(f, obj.qrb);
end

% reaction
if ~isempty(obj.c)
    f = CompositeFunction(@BilinearForm.reactionPart, obj.c, phi);
    data = data + integrateElement(f, obj.qrc);
end

end
