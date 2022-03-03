% vectorProduct Product of two objects represented by row-wise vectors
%
%   The columns of the input arguments represent matrices of dimension [ma, ra]
%   and [rb, mb], respectively, which are stored in a column-wise fashion. This
%   function multiplies these matrices, resulting in a [ma, mb] matrix, which is
%   also stored column-wise. The dimensions of the matrices can be given as
%   optional row-vector parameter. Transposing the dimension to a column vector
%   parameter results in the matrix to be interpreted as transposed.
%   This function can simultaneously multiply ka and kb (columns of A and B) of
%   such matrices, where ka=1, kb=1, or ka=kb.
%   Additional dimensions (3 and onward) are preserved.
%
%   val = vectorProduct(A, B) Take column-wise dot product of A and B
%
%   val = vectorProduct(A, B, dimA, dimB) Take column-wise product of matrices A
%       and B with dimensions dimA and dimB (row-vectors), respectively. If dimA
%       is transposed, this realizes A'*B column-wise. Analogously for dimB.

function val = vectorProduct(A, B, dimA, dimB)

arguments
    A double
    B double
    dimA double = [size(A, Dim.Vector), 1]';
    dimB double = [size(B, Dim.Vector), 1];
end

[trA, nA] = checkForTransposition(dimA);
[trB, nB] = checkForTransposition(dimB);

assert(length(dimA) == 2 && length(dimB) == 2, 'vectorProduct:dimensionMismatch', ...
    'Dimension parameters must be two-vectors.')
assert(nA(2) == nB(1), 'vectorProduct:dimensionMismatch', ...
    'Dimensions must be suitable for matrix-matrix multiplication.')
assertCorrectFirstDimension(A, dimA)
assertCorrectFirstDimension(B, dimB)

n = max(ndims(A), ndims(B));
sizeAB = max(padSize(A, n), padSize(B, n));
val = reshape( ...
        pagemtimes(prepare(A, dimA), trA, prepare(B, dimB), trB), ...
      [nA(1)*nB(2), sizeAB(2:end)]);
end

%% local functions
function assertCorrectFirstDimension(mat, dim)
    if size(mat, Dim.Vector) ~= prod(dim)
        eidType = 'vectorProduct:dimensionMismatch';
        msgType = 'Number of rows must match given dimension.';
        throwAsCaller(MException(eidType,msgType))
    end
end

function [transposeString, n] = checkForTransposition(dim)
    n = dim(:)';
    isTransposed = (size(dim, Dim.Vector) == 2);
    if isTransposed
        n = flip(n);
        transposeString = 'transpose';
    else
        transposeString = 'none';
    end
end

function res = padSize(mat, n)
    res = [size(mat), zeros(1,n-ndims(mat))];
end

function res = prepare(mat, dim)
    actualSize = size(mat);
    res = reshape(mat, [dim(:)', actualSize(2:end)]);
end
