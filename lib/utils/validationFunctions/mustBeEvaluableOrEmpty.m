% mustBeEvaluableOrEmpty Check if child of Evaluable or empty
%
%   mustBeEvaluableOrEmpty(a)

function mustBeEvaluableOrEmpty(a)

if ~(isempty(a) || isa(a, 'Evaluable'))
    eid = 'mustBeEvaluableOrEmpty:NotEvaluableOrEmpty';
    msg = 'Coefficients must be either empty or subclasses of Evaluable.';
    throwAsCaller(MException(eid,msg))
end

end
