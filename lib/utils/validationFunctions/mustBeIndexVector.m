% mustBeIndexVector Check if parameter is valid for indexing
%
%   mustBeIndexVector(a)

function mustBeIndexVector(a)

if ~(isempty(a) ...
     || (isvector(a) && (islogical(a) || isnumeric(a))) ...
     || (isscalar(a) && a==':'))
    eid = 'mustBeIndexVector:notAnIndexVector';
    msg = 'Index vector must be either vector of logicals, vector of numericals, or '':''.';
    throwAsCaller(MException(eid,msg))
end

end