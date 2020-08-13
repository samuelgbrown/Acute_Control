function sdStrs = calculateSDStrs(firingProbFcn, sdCurveProb, durs, str0, numFcns)
% This function will calculate the SD curve for a given neuron
% firingProbFcn = @(GT, thetaInd)
% sdCurveProb = .5, usually
% durs = linspace(durMin, durMax, numDurs)
% str0 = 30, usually
% numFcns = 1, usually (upper limit of indices entered into firingProbFcn

probDiffFromSD = .01; % The probability calculated for each corrupted SD curve must be at most this far away from the target SD curve value
% diffCheckDist = 1e-2;

% numFcns = size(firingProbFcn, 1);
numSDDurs = length(durs);

if length(str0) == 1 && numFcns > 1
    str0 = str0*ones(1, numFcns);
end

str0Default = str0; % Save the initial strength, to return to later if neededs
ms = MultiStart;

sdStrs = zeros(numFcns, numSDDurs);
for durInd = 1:numSDDurs
    if mod(durInd, round(durs/10)) == 0
        fprintf('Starting duration %d out of %d...\n', durInd, numSDDurs);
    end
    
    % Get this duration
    %         thisDur = allDursMesh(sdIndsToUse(durIndInd));
    thisDur = durs(durInd);
    
    thisSs = nan(numFcns, 1);
    for fcnNum = 1:numFcns
        %         for i = 1:10
        % %             testProb1 = approxFun([(str0(fcnNum)) thisDur], theta(fcnNum, :));
        % %             testProb2 = approxFun([(str0(fcnNum) + diffCheckDist) thisDur], theta(fcnNum, :));
        %
        %             testProb1 = firingProbFcn([(str0(fcnNum)) thisDur], fcnNum);
        %             testProb2 = firingProbFcn([(str0(fcnNum) + diffCheckDist) thisDur], fcnNum);
        %
        %             if testProb1 ~= testProb2
        %                 % The two values are different enough that a
        %                 % gradient is defined
        %                 break;
        %             else
        %                 % Not gradient is defined at this point, try a
        %                 % smaller strength
        %                 str0(fcnNum) = str0(fcnNum)/2;
        %             end
        %         end
        
        % Calculate the strength for this duration and this theta
        %         thisSs(fcnNum) = fminunc(@(s) (approxFun([s thisDur], theta(fcnNum, :)) - sdCurveProb)^2, str0(fcnNum));
        %         thisP = approxFun([thisSs(fcnNum) thisDur], theta(fcnNum, :)); % Get the probability of firing
        
        prob = createOptimProblem('fminunc', 'objective', @(s) (firingProbFcn([s thisDur], fcnNum) - sdCurveProb)^2, 'x0', str0(fcnNum));
        allStartStr = CustomStartPointSet([linspace(.01, 2, 10)';str0(fcnNum)]);
        thisSs(fcnNum) = run(ms, prob, allStartStr);
        
        %         thisSs(fcnNum) = fminunc(@(s) (firingProbFcn([s thisDur], fcnNum) - sdCurveProb)^2, );
        thisP = firingProbFcn([thisSs(fcnNum) thisDur], fcnNum); % Get the probability of firing
        
        if abs(thisP - sdCurveProb) > probDiffFromSD
            % If the strength calculated doesn't REALLY
            % give the SD curve probability level, then
            % throw it out (we have a minimum, not a zero)
            thisSs(fcnNum) = nan;
        end
    end
    
    % Save this duration in the SD curve for all functions
    sdStrs(:, durInd) = thisSs;
    
    % If any of the strengths could not be found for this duration, return
    % them to the default initial strength, for the next loop
    if any(isnan(str0))
        str0(isnan(str0)) = str0Default;
    end
end
end