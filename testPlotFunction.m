%% testPlotFunction
%
%  This function operates in a similar manner to student plot checker.  It takes
%  in a string or a cell array of two strings and compares name name_soln or the
%  two given strings.  Passed along are the rest of the arguments to both
%  functions.
%
%  Inputs:
%    1. (char or cell) Either the name of a function, or the name of two
%      functions.
%    2-n. (any) The arguments to be passed to the other functions
%
%  Outputs:
%    1. (logical) Whether or not the function's plot outputs match
%    2. (char) Detailed output passed from testFigure

function [out, outStr] = testPlotFunction(name, varargin)
  if ischar(name)
    testName = name;
    refName = [name '_soln'];
  elseif iscell(name) && length(name) == 2 && ischar(name{1}) && ischar(name{2})
    testName = name{1};
    refName = name{2};
  else
    error 'Function name is invalid.  Require either a string or two strings in a cell array'.
  end

  testFunc = str2func(testName);
  refFunc = str2func(refName);

  testHandle = figure;
  testFunc(varargin{:});

  % have to do this in case the function does something unexpected, like creates
  % its own handle... (I'm looking at you, plotShapes_soln)
  testHandleReal = gcf;

  refHandle = figure;
  refFunc(varargin{:});

  refHandleReal = gcf;

  % default to displaying 20 outputs
  [out, outStr] = testFigure(testHandleReal, refHandleReal, 20);

  % clean up a bit
  if (isvalid(testHandle))
    close(testHandle);
  end

  if (isvalid(testHandleReal))
    close(testHandleReal);
  end

  if (isvalid(refHandle))
    close(refHandle);
  end

  if (isvalid(refHandleReal))
    close(refHandleReal);
  end
end
