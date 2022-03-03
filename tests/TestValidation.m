classdef TestValidation < matlab.unittest.TestCase

methods (Test)
    function constantIsEvaluable(testCase)
        testCase.verifyWarningFree(@() mustBeEvaluableOrEmpty(Constant([], 1)))
    end
    
    function emptyVectorIsEmpty(testCase)
        testCase.verifyWarningFree(@() mustBeEvaluableOrEmpty([]))
    end
    
    function vectorIsValidIndexVector(testCase)
        testCase.verifyWarningFree(@() mustBeIndexVector([1,2,3]))
    end
    
    function logicalsAreValidIndexVector(testCase)
        testCase.verifyWarningFree(@() mustBeIndexVector([true,false,true]))
    end
    
    function colonIsValidIndexVector(testCase)
        testCase.verifyWarningFree(@() mustBeIndexVector(':'))
    end
    
    function numberIsNotEvaluable(testCase)
        testCase.verifyError(@() mustBeEvaluableOrEmpty(1), 'mustBeEvaluableOrEmpty:NotEvaluableOrEmpty')
    end 
       
    function matrixIsNotIndexVector(testCase)
        testCase.verifyError(@() mustBeIndexVector([1,2;3,4]), 'mustBeIndexVector:notAnIndexVector')
    end
    
    function stringsAreNotIndexVectors(testCase)
        testCase.verifyError(@() mustBeIndexVector(['a', 'b']), 'mustBeIndexVector:notAnIndexVector')
    end
end
    
end
