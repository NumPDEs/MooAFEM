classdef TestVectorProduct < matlab.unittest.TestCase
    
properties (TestParameter)
    A = struct('vector', [1;2], 'matrix', [1,2;3,4]);
    B = struct('vector', [1;2], 'matrix', [1,2;3,4]);
end

methods (Test, ParameterCombination='exhaustive')
    function repmatAndMultiplicationCommutate(testCase, A, B)
        [Aext, dimA] = testCase.extendAlongDim(A, 2, 10);
        [Bext, dimB] = testCase.extendAlongDim(B, 2, 10);
        actual = vectorProduct(Aext, Bext, dimA', dimB);
        expected = repmat(asVector(A' * B), 1, 10);
        testCase.verifyEqual(actual, expected, 'AbsTol', eps)
    end

    function sanityCheckDetectsDimensionMismatch(testCase)
        testCase.verifyError(@() vectorProduct([1;2], [1;2;3]), ...
            'vectorProduct:dimensionMismatch')
        testCase.verifyError(@() vectorProduct([1;2;3], [1;2;3], [3,1], [3,1]), ...
            'vectorProduct:dimensionMismatch')
    end

    function outerProductWorksCorrectly(testCase)
        u = [1;2];
        v = [1;2;3];
        [Aext, ~] = testCase.extendAlongDim(u, 2, 10);
        [Bext, ~] = testCase.extendAlongDim(v, 2, 10);
        actual = vectorProduct(Aext, Bext, [2,1], [1,3]);
        expected = repmat(asVector(u * v'), 1, 10);
        testCase.verifyEqual(actual, expected, 'AbsTol', eps)
    end

    function trailingDimensionsHandledCorrectly(testCase, A, B)
        d = (numel(A)+numel(B))/2;
        [Aext, dimA] = testCase.extendAlongDim(A, 2, 10);
        [Bext, dimB] = testCase.extendAlongDim(B, d, 10);
        repetitions = ones(1, d);
        repetitions([2,d]) = 10;
        actual = vectorProduct(Aext, Bext, dimA', dimB);
        expected = repmat(asVector(A' * B), repetitions);
        testCase.verifyEqual(actual, expected, 'AbsTol', eps)
    end
end

methods (Access='protected')
    function [Aext, dimA] = extendAlongDim(~, A, dim, n)
        dimA = size(A);
        repetitions = ones(1, dim);
        repetitions(end) = n;
        Aext = repmat(asVector(A), repetitions);
    end
end
    
end