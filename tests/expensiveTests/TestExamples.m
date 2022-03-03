classdef TestExamples < matlab.unittest.TestCase
    
properties (TestParameter)
    example = TestExamples.getExamplesInDirectory();
end

methods (Test)
    function exampleRunsWithoutErrorOrWarning(testCase, example)
        evalc(example);
        close all
        testCase.verifyTrue(true);
    end
end

methods (Static)
    function examples = getExamplesInDirectory()
        myDirectory = fileparts(mfilename('fullpath'));
        rootDirectory = fileparts(fileparts(myDirectory));
        files = {dir(fullfile(rootDirectory, 'examples', '*.m')).name};
        [~, examples, ~] = fileparts(files);
        examples = repelem(examples, 1, 2);
        examples = struct(examples{:});
    end
end
    
end