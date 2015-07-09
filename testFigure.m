%% testFigure
%
%  This function is intended to test two figures for identical axes and similar
%  points graphed, along with line style and color.
%
%  Inputs:
%    1. (Figure) The figure to test
%    2. (Figure) The reference figure
%    3. (int32) Stop after finding this many errors, default is 5.  -1 is
%       infinite
%    4. (uint32) An integer representing the number of epsilon to accept,
%       default is 1000.
%    5. (double) A double representing the largest difference to reject
%
%  Outputs:
%    1. (bool) true if the two figures match
%    2. (char) A detailed description of what matches and what does not.

function [out, outStr] = testFigure(testFig, refFig, keepGoing, epsilons, absoluteDiff)

  if ~isequal(class(testFig), 'matlab.ui.Figure') ||...
     ~isequal(class(refFig),  'matlab.ui.Figure')
    error 'Unexpected class type passed in as arguments.';
  end

  if nargin < 3
    keepGoing = 5;
  end

  if nargin < 4
    epsilons = 1000;
  end

  % this is required because sometimes when we get near 0 we're still a little
  % off, but the eps on these values is so low that our eps difference is
  % ridiculous (10^12 territory)
  if nargin < 5
    absoluteDiff = 10^-15;
  end

  % if user enters a 0, we want to fail after one
  if ~keepGoing
    keepGoing = 1;
  end

  if epsilons <= 1
    epsilons = 1;
  end

  epsilons = floor(epsilons);

  % grab the axes.  There may be several of each
  testAxes = get(testFig, 'Children');
  refAxes = get(refFig, 'Children');

  if ~isequal(class(testAxes), 'matlab.graphics.axis.Axes') ||...
     ~isequal(class(refAxes),  'matlab.graphics.axis.Axes')
    error 'Unexpected class type inside of figures.';
  end

  % innocent until proven guilty
  out = true;
  outStr = '';

  % next, if the number of axes are the different, set out to false and an outStr
  % explaining the issue
  if length(testAxes) ~= length(refAxes)
    out = false;
    outStr = [outStr sprintf(['The number of subplots are different.  The ',...
      'test figure has %d plots, but the reference figure has %d.\n'],...
      length(testAxes), length(refAxes))];
    keepGoing = keepGoing - 1;
    if ~keepGoing
      return;
    end
  end

  % the number of axes to test
  axesToTest = min(length(testAxes), length(refAxes));

  for ii = 1:length(axesToTest)
    [axesOut, axesStr, keepGoing] = testAxis(testAxes(ii), refAxes(ii),...
      epsilons, absoluteDiff, keepGoing);
    out = out && axesOut;
    if ~axesOut
      outStr = [outStr sprintf('Subplot %d:\n', ii) axesStr];
      if ~keepGoing
        return;
      end
    end
  end

  if out
    outStr = [outStr sprintf('Everything seems to match!  Congratulations!\n')];
  end
end

% testAxis tests each graph in a similar manner to testFigure
function [out, outStr, keepGoing] = testAxis(testAxes, refAxes, epsilons, absoluteDiff, keepGoing)

  out = true;
  outStr = '';

  % first, test the limits of the axes
  testXLim = get(testAxes, 'XLim');
  refXLim = get(refAxes, 'XLim');

  testYLim = get(testAxes, 'YLim');
  refYLim = get(refAxes, 'YLim');
  
  testZLim = get(testAxes, 'ZLim');
  refZLim = get(refAxes, 'ZLim');

  [xOut, xOutStr] = testLimits(testXLim, refXLim, 'XLim', epsilons, absoluteDiff);
  out = out && xOut;
  outStr = [outStr xOutStr];

  if ~xOut
    keepGoing = keepGoing - 1;
    if ~keepGoing
      return;
    end
  end

  [yOut, yOutStr] = testLimits(testYLim, refYLim, 'YLim', epsilons, absoluteDiff);

  out = out && yOut;
  outStr = [outStr yOutStr];
  
  if ~yOut
    keepGoing = keepGoing - 1;
    if ~keepGoing
      return;
    end
  end

  if length(refZLim)
    [zOut, zOutStr] = testLimits(testZLim, refZLim, 'ZLim', epsilons, absoluteDiff);
    
    out = out && zOut;
    outStr = [outStr zOutStr];

    if ~zOut
      keepGoing = keepGoing - 1;
      if ~keepGoing
        return;
      end
    end
  end

  % now we test the view

  testView = get(testAxes, 'View');
  refView = get(refAxes, 'View');

  if ~isequal(testView, refView)
    out = false;

    testViewText = sprintf('[%f, %f]', testView(1), testView(2));
    refViewText = sprintf('[%f, %f]', refView(1), refView(2));

    % view(2) special case
    if (isequal(testView, [0 90]))
      testViewText = 'view(2)';
    end

    if (isequal(refView, [0 90]))
      refViewText = 'view(2)';
    end

    % view(3) special case
    if (isequal(testView, [-37.5 30]))
      testViewText = 'view(3)';
    end

    if (isequal(refView, [-37.5 30]))
      refViewText = 'view(3)';
    end
    
    outStr = [outStr sprintf('View mismatch.  Got %s, expected %s.\n',...
      testViewText, refViewText)];
      
    keepGoing = keepGoing - 1;
    if ~keepGoing
      return;
    end
  end

  % now we test the aspect ratio, for axis square, axis equal etc.

  testDataAsp = get(testAxes, 'DataAspectRatio');
  refDataAsp = get(refAxes, 'DataAspectRatio');

  testPlotAsp = get(testAxes, 'PlotBoxAspectRatio');
  refPlotAsp = get(refAxes, 'PlotBoxAspectRatio');

  if ~isequal(testDataAsp, refDataAsp) || ~isequal(testPlotAsp, refPlotAsp)
    out = false;
    outStr = [outStr sprintf(['Aspect ratio mismatch.  Did you forget to '...
      'set axis square or axis equal or something similar?\n'])];
  end

  % next, get every line and test for sanity

  testLines = get(testAxes, 'Children');
  refLines = get(refAxes, 'Children');
  
  if ~isequal(class(testLines), 'matlab.graphics.chart.primitive.Line') ||...
     ~isequal(class(refLines),  'matlab.graphics.chart.primitive.Line')
    error 'Unexpected class type inside of axes.';
  end

  % next, if the number of axes are the different, set out to false and an outStr
  % explaining the issue
  if length(testLines) ~= length(refLines)
    out = false;
    outStr = [outStr sprintf(['The number of lines are different.  The ',...
      'test axis has %d lines, but the reference axis has %d.\n'],...
      length(testLines), length(refLines))];
    keepGoing = keepGoing - 1;
    if ~keepGoing
      return;
    end
  end

  linesToTest = min(length(testLines), length(refLines));

  for ii = 1:linesToTest
    [lineOut, lineStr, keepGoing] = testLine(testLines(ii), refLines(ii),...
      epsilons, absoluteDiff, keepGoing);
    out = out && lineOut;
    
    if ~lineOut
      outStr = [outStr sprintf('Line %d:\n', ii) lineStr];

      if ~keepGoing
        return;
      end
    end
  end
end

% test limit tests an axes limit pair
function [out, outStr] = testLimits(testLim, refLim, limName, epsilons, absoluteDiff)
  out = true;
  outStr = '';

  [lowOut, lowOutStr] = testLimit(testLim(1), refLim(1),...
    ['Low' limName], epsilons, absoluteDiff);
  out = out && lowOut;
  outStr = [outStr lowOutStr];

  [highOut, highOutStr] = testLimit(testLim(2), refLim(2),...
    ['High' limName], epsilons, absoluteDiff);

  out = out && highOut;
  outStr = [outStr highOutStr];

end

% test one limit
function [out, outStr] = testLimit(testLim, refLim, limName, epsilons, absoluteDiff)
  out = true;
  outStr = '';

  [xEq, xEps] = almostEqual(testLim, refLim, epsilons, absoluteDiff);
  if ~xEq
    out = false;
    if xEps < 100 * epsilons
      outStr = [outStr sprintf(['%s is near expected value, a difference of'...
        ' %d epsilons.\n'], limName, xEps)];
    else
      outStr = [outStr sprintf('%s is incorrect. Expected %0.6f got %0.6f.\n',...
        limName, testLim, refLim)];
    end
  end
end

% test one line
function [out, outStr, keepGoing] = testLine(testLine, refLine, epsilons, absoluteDiff, keepGoing)

  out = true;
  outStr = '';

  % first test the properties of the line to make sure they're the same.
  testColor = get(testLine, 'Color');
  refColor = get(refLine, 'Color');

  if ~isequal(testColor, refColor)
    out = false;
    outStr = [outStr sprintf('Colors of line differ.\n')];
    keepGoing = keepGoing - 1;
    if ~keepGoing
      return;
    end
  end

  % and the line type
  testType = get(testLine, 'LineStyle');
  refType = get(refLine, 'LineStyle');

  if ~isequal(testType, refType)
    out = false;
    outStr = [outStr sprintf(['Style of the lines differ.  Got ''%s'' '...
      'expected ''%s''.\n'], testType, refType)];
    keepGoing = keepGoing - 1;
    if ~keepGoing
      return;
    end
  end

  % finally, test the points
  testX = get(testLine, 'XData');
  testY = get(testLine, 'YData');
  testZ = get(testLine, 'ZData');
  testPoints = [];
  if ~length(testZ)
    testPoints = [testX; testY;];
  else
    testPoints = [testX; testY; testZ;];
  end

  refX = get(refLine, 'XData');
  refY = get(refLine, 'YData');
  refZ = get(refLine, 'ZData');
  refPoints = [];
  if ~length(refZ)
    refPoints = [refX; refY;];
  else
    refPoints = [refX; refY; refZ;];
  end
 
  if ~isequal(size(refPoints), size(testPoints))
    out = false;
    outStr = [outStr sprintf('The number of dimensions disagree.\n')];
    keepGoing = keepGoing - 1;
    return; % don't need to check, we're done with this line
  end

  zEnabled = length(testZ) ~= 0;
  if isempty(testX)
    if ~isempty(refX)
      out = false;
      outStr = [outStr sprintf(['Empty data in line where data was '...
        'expected.\n'])];
      keepGoing = keepGoing - 1;
      if ~keepGoing
        return;
      end
    end
  % dataY will always be the same length as dataX, so we only need to test X
  elseif length(testX) ~= length(refX)
    out = false;
    outStr = [outStr sprintf(['Data length mismatch.  Got length %d'...
      ' expected %d.\n'], length(testX), length(refX))];
    keepGoing = keepGoing - 1;
    if ~keepGoing
      return;
    end
  else
    % test polygon in this case
    if refPoints(:,1) == refPoints(:,end)
      if testPoints(:,1) ~= testPoints(:,end)
        out = false;
        outStr = [outStr 'Poly not closed.  Last point does not match first '...
          'point.\n'];
        keepGoing = keepGoing - 1;
        if ~keepGoing
          return;
        end
      else
        % strip out the end points
        testPoints = testPoints(:,1:end-1);
        refPoints = refPoints(:,1:end-1);
        testX = testX(1:end-1);
        testY = testY(1:end-1);
        refX = refX(1:end-1);
        refY = refY(1:end-1);

        if zEnabled
          testZ = testZ(1:end-1);
          refZ = refZ(1:end-1);
        end

        dataLen = length(testX);

        % now figure out where, if anywhere the points match to come up with an
        % offset

        testOffset = 0;

        for ii = 1:dataLen
          [~, xEps] = almostEqual(testX(ii), refX, epsilons, absoluteDiff);
          [~, yEps] = almostEqual(testY(ii), refY, epsilons, absoluteDiff);
          offsetsX = find(almostEqual(testX(ii), refX, epsilons, absoluteDiff));
          offsetsY = find(almostEqual(testY(ii), refY, epsilons, absoluteDiff));

          offsets = intersect(offsetsX, offsetsY);

          if zEnabled
            offsetsZ = find(almostEqual(testZ(ii), refZ, epsilons, absoluteDiff));
            offsets = intersect(offsets, offsetsZ);
          end

          if ~isempty(offsets)
            testOffset = ii - 1;
            break;
          end
        end

        if isempty(offsets)
          out = false;
          outStr = [outStr sprintf(['Could not find a match for any '...
            'data point to start the polygon in the reference data.\n'])];

          keepGoing = keepGoing - 1;
          if ~keepGoing
            return;
          end
        else
          best = 0;
          bestEps = 0;
          bestStr = '';
          matchFound = false;

          % for each match
          for offset = offsets

            % if we've found a match on a previous iteration we don't need to do
            % this anymore
            if matchFound
              break;
            end

            % and for each direction
            for direction = [-1 1]
              % does everything match?

              testRot = mod((1:dataLen) + testOffset - 1, dataLen) + 1;

              % this rotation is a bit complicated, I have to run it 0 indexed to
              % get the math right, which is where the -1 comes from.
              % Other than that it's reasonably simple.
              refRot = mod(...
                (0:direction:(direction*(dataLen-1))) + offset - 1,...
                dataLen) + 1;

              thisTestX = testX(testRot);
              thisTestY = testY(testRot);
              thisRefX = refX(refRot);
              thisRefY = refY(refRot);

              thisTestPoints = [thisTestX; thisTestY];
              thisRefPoints = [thisRefX; thisRefY];

              if zEnabled
                thisTestZ = testZ(testRot);
                thisRefZ = refZ(refRot);
                thisTestPoints = [thisTestPoints; thisTestZ];
                thisRefPoints = [thisRefPoints; thisRefZ];
              end

              [xEq, xEps] = almostEqual(thisTestPoints, thisRefPoints, epsilons, absoluteDiff);

              matches = all(xEq, 1);
              maxEps = max(max(xEps));
              count = sum(matches);
              if count > best || (count == best && maxEps < bestEps)
                best = count;
                bestEps = maxEps;
                bestStr = '';
                misses = find(~matches);
                for ii = misses
                  if all(xEps(:,ii) < 100 * epsilons)
                    zCoord = '';
                    if zEnabled
                      zCoord = sprintf(', %0.17f', thisRefZ(ii));
                    end
                    bestStr = [bestStr,...
                      sprintf(['Data mismatch, index %d is near expected '...
                        'value of [%0.17f, %0.17f%s], with a max difference '...
                        'of %d epsilons.\n'],...
                        mod(ii + testOffset, dataLen) + 1, thisRefX(ii),...
                        thisRefY(ii), zCoord, maxEps)];
                  else
                    zCoordTest = '';
                    zCoordRef = '';
                    if zEnabled
                      zCoordTest = sprintf(', %0.6f', thisTestZ(ii));
                      zCoordRef = sprintf(', %0.6f', thisRefZ(ii));
                    end
                    bestStr = [bestStr,...
                      sprintf(['Data mismatch, index %d is incorrect. '...
                        'Expected [%0.6f, %0.6f%s] got [%0.6f, %0.6f%s].\n'],...
                        mod(ii + testOffset, dataLen), thisRefX(ii),...
                        thisRefY(ii), zCoordRef, thisTestX(ii),...
                        thisTestY(ii), zCoordTest)];
                  end
                end
              end % if count > best

              % if we have a full match
              if best == dataLen
                matchFound = true;
                break;
              end

            end % for direction = [-1 1]
          end % for offset = offsets

          out = out && matchFound;
          if ~matchFound
            outStr = [outStr bestStr]
            keepGoing = keepGoing - 1;
            if ~keepGoing
              return;
            end
          end

        end % if isempty(offsets)
      end
    else % if refPoints(:,1) == refPoints(:,end)
      % if we're not testing direction and offset agnostic then the situation is
      % much easier

      for ii = 1:length(testX)
        [xEq, xEps] = almostEqual(testPoints(:, ii), refPoints(:, ii), epsilons, absoluteDiff);

        if ~all(xEq)
          out = false;

          maxEps = max(xEps(1:end));

          if maxEps < 100 * epsilons
            zCoord = '';
            if zEnabled
              zCoord = sprintf(', %0.17f', refZ(ii));
            end

            outStr = [outStr,...
              sprintf(['Data mismatch, index %d is near expected '...
                'value of [%0.17f, %0.17f%s], with a max difference of '...
                '%d epsilons.\n'],...
                ii, refX(ii), refY(ii), zCoord, max(xEps(ii), yEps(ii)))];
          else
            zCoordTest = '';
            zCoordRef = '';
            if zEnabled
              zCoordTest = sprintf(', %0.6', testZ(ii));
              zCoordRef = sprintf(', %0.6', refZ(ii));
            end

            outStr = [outStr,...
            sprintf(['Data mismatch, index %d is incorrect. '...
              'Expected [%0.6f, %0.6f%s] got [%0.6f, %0.6f%s].\n'],...
              ii, refX(ii), refY(ii), zCoordRef,...
              testX(ii), testY(ii), zCoordTest)];
          end

          keepGoing = keepGoing - 1;
          if ~keepGoing
            return;
          end
        end
      end
    end % if refPoints(:,1) == refPoints(:,end)
  end % if isempty(testX)
end

% almostEqual tests floats for near equality.
function [out, epsDiff] = almostEqual(a, b, epsilons, absoluteDiff)
  diff = abs(a - b);
  epsDiff = floor(diff ./ (epsilons * eps(b)));
  out = epsDiff <= epsilons;
  out = out | diff < absoluteDiff;
end
