import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin

suite = TestSuite.fromFolder('test');

runner = TestRunner.withTextOutput;
runner.addPlugin(CodeCoveragePlugin.forFolder(string(pwd) + filesep + "code"))
result = runner.run(suite);

% Note where the results are stored, and run (replace temp name):
%copyfile('/tmp/tp28ccf01e_b
% d54_48d1_a992_5e69865dbf1a/index.html', pwd)
% for example to bring it here. This can be read then in a browser to check
% code coverage and also run all the tests in the 'test' folder.