% assemblePatchwise Assembles matrix of bilinear form patchwise.
%
%   A = assemblePatchwise(obj, fes) assembles the data of the given
%       BilinearForm for a FeSpace fes and returns a PatchwiseMatrix A.
%
%   See also: BilinearForm

function mat = assemblePatchwise(obj, fes)
arguments
    obj
    fes FeSpace
end

if ~isempty(obj.robin)
    error('Patchwise assembly not implemented for Robin boundary!')
end

% in the long run, this should work just with non-assembled data:
data = computeVolumeData(obj, fes);
mat = PatchwiseMatrix(fes, data);

end