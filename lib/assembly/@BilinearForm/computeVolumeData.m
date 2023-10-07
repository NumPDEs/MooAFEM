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
p = fes.finiteElement.order;

% diffusion
if ~isempty(obj.a)
    qr = QuadratureRule.ifEmptyOfOrder(obj.qra, 2*(p-1));
    f = CompositeFunction(@BilinearForm.diffusionPart, obj.a, Dphi);
    data = data + integrateElement(f, qr);
end

% convection
if ~isempty(obj.b)
    qr = QuadratureRule.ifEmptyOfOrder(obj.qrb, 2*p-1);
    f = CompositeFunction(@BilinearForm.convectionPart, obj.b, Dphi, phi);
    data = data + integrateElement(f, qr);
end

% reaction
if ~isempty(obj.c)
    qr = QuadratureRule.ifEmptyOfOrder(obj.qrc, 2*p);
    f = CompositeFunction(@BilinearForm.reactionPart, obj.c, phi);
    data = data + integrateElement(f, qr);
end

end
