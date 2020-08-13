function [xMin, yMin, pMin, endStr] = minimize2DSpline(P, durs, strs, durStartSub, strStartSub)
% This function will find the minimum point in a space defined by x and y,
% that is undersampled by matrix P.

tic;
size = 100;

durStart = durs(durStartSub);
strStart = strs(strStartSub);

durSmallInds = max((durStartSub - size),1):min((durStartSub + size),length(durs));
strSmallInds = max((strStartSub - size),1):min((strStartSub + size),length(strs));

PSmall = P(strSmallInds,durSmallInds);
durSmall = durs(durSmallInds);
strSmall = strs(strSmallInds);

options = optimoptions('fmincon', 'Display', 'off');
[X, pMin] = fmincon(@(pIn) interp2(durSmall,strSmall,PSmall,pIn(1),pIn(2), 'cubic'), [durStart strStart], [], [], [], [], [durSmall(1) strSmall(1)], [durSmall(end) strSmall(end)], [], options);

% fprintf('Took %0.2f seconds.\n', toc);
xMin = X(1);
yMin = X(2);
endStr = '';
% % Settings
% indsAround = 100; % Number of indicies around the center point of each location to take
% debug = false; % Output information on how close to tolerance each 
% 
% % Tolerances
% fTol = 1e-20;
% distTol = 1e-20;
% 
% % Find initial starting place
% maxX = length(x);
% maxY = length(y);
% minVal = min(P(:));
% minInd = find(P(:) == minVal, 1);
% [lastMinSub, last2MinSub] = ind2sub(size(P), minInd);
% lastMinLoc = y(lastMinSub); % The location along the other axis of the previous minimum
% last2MinLoc = x(last2MinSub); % The location along this axis of the minimum from two runs ago
% lastMinVal = inf;
% options = optimoptions('fmincon', 'Display', 'off');
% 
% curAxisX = true; % Is the current axis the x axis?
% thisAxisDone = false(2,1);
% count = 1;
% while true
%     if curAxisX
%         lastMinSubVec = max(lastMinSub-indsAround,1):min(lastMinSub+indsAround, maxY);
%         last2MinSubVec = max(last2MinSub-indsAround,1):min(last2MinSub+indsAround, maxX);
%         thisAxis = x;
%         otherAxis = y;
%         thisValsMat = P(lastMinSubVec, last2MinSubVec)';
%     else
%         lastMinSubVec = max(lastMinSub-indsAround,1):min(lastMinSub+indsAround, maxX);
%         last2MinSubVec = max(last2MinSub-indsAround,1):min(last2MinSub+indsAround, maxY);
%         thisAxis = y;
%         otherAxis = x;
%         thisValsMat = P(last2MinSubVec, lastMinSubVec);
%     end
%     thisSubAxis = thisAxis(last2MinSubVec); % The axis that the minimum will be found along
%     otherSubAxis = otherAxis(lastMinSubVec); % The opposite axis
%     
%     splineVals = zeros(size(thisValsMat,1),1);
%     for i = 1:size(thisValsMat,1)
%         % For each row along the axis that will be optimized
%         pp = spline(otherSubAxis, thisValsMat(i,:));
%         splineVals(i) = ppval(pp, lastMinLoc);
%     end
%     %% Check this axis
%     pp = spline(thisSubAxis, splineVals);
%     
%     startPos = max(min(last2MinLoc, thisSubAxis(end)), thisSubAxis(1));
%     thisMinLoc = fmincon(@(x) ppval(pp, x), startPos, [], [], [], [], thisSubAxis(1), thisSubAxis(end), [], options);
%     thisMinVal = ppval(pp, thisMinLoc);
%     
%     %% Check if tolerances have been met
%     if debug
%         fprintf('Iter number %d:\nDistance is %0.14f\n',count, abs(thisMinLoc - last2MinLoc));
%     end
%     
%     if abs(thisMinLoc - last2MinLoc) < distTol && ~thisAxisDone(curAxisX + 1)
%         % Correct this, consider optimizing along border case
%         thisAxisDone(curAxisX + 1) = true;
%         
%         if all(thisAxisDone)
%             endStr = 'Distance change between minimums';
%             break;
%         end
%     end
%     
%     if debug
%         fprintf('Difference is %0.14f\n\n',abs(thisMinVal - lastMinVal));
%         count = count + 1;
%     end
%     
%     if abs(thisMinVal - lastMinVal) < fTol
%         endStr = 'Change in P';
%         break;
%     end
%     
%     %% Prepare for next loop
%     lastMinVal = thisMinVal;
%     last2MinLoc = lastMinLoc;
%     lastMinLoc = thisMinLoc;
%     
%     last2MinSub = lastMinSub;
%     lastMinSub = find(min(abs(lastMinLoc - thisAxis)) == abs(lastMinLoc - thisAxis), 1);
%     
%     curAxisX = ~curAxisX; % Toggle which axis is being minimized
% end
% 
% if curAxisX
%     yMin = lastMinLoc;
%     xMin = thisMinLoc;
% else
%     yMin = thisMinLoc;
%     xMin = lastMinLoc;
% end
% 
% pMin = thisMinVal;