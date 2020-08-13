% This script will model a neuron's SD curve (neuron A) being characterized
% data being corrupted by a second neuron (neuron B).

clc;
close all;

% global PTableDurationLookup PTableLookup ProbOrCalc diffFun

ProbUnion = @(P1, P2) P1 + P2 - (P1.*P2);
diffFun = @(avgProb1, avgProb2) avgProb2 - avgProb1 + (1 - avgProb1).^2;

%% Parameters for simulation
% Time constants: simulation t is about 1/3 of real t
% Inverse alpha should be on order of 10s
% Membrane capacitance is pretty constant from cell to cell - usually
% F/cm^2, diameter is about 10-20um (can assume sphere).
% Choose g to be non-dimensional, beta has units of 1/s. What beta*g should
% be rheobase?
% Calibrate noise in FP:
% FP: Set initial condition to steady state of OU process, cut off at
% boundaries, then set mass to 1.  Distribution will change according to
% alpha.  Fix sigma according to smallest alpha (which will have the most
% spontaneous spiking)
% Check if there is low error in running FP for only the length of the
% stimulus (try running FP on a number of G's and T's for the full 15s or
% so, then plot mass and Phi [value of last bin on boundary] over time).


isProduction = false; % (NOT USED) Are we in production mode (i.e. incredibly precise calculations) or testing mode (faster calculations)
useFPInterpolation = true; % Should a FP table be used to find the simulation values? (alternative: the stimulation values being used will be simulated in FP directly, and saved to disk)
useFPCostOptimization = true; % Should the newer version of cost optimization (to find the best GT for each pair of neurons) be used, instead of the original version?
useAutoBetaOpt = true; % Should the beta values for neuron C be calculated (alternative: they are set manually)
betaPlotTest = false; % Should a beta test plot be made?
doPlotSecondCurve = false; % Should the second SD curve be plotted (in the case of a doubled curve, which sometimes happens in the deleted spikes case)
useMultistart = true; % Should multistart global optimization be used?

doPlot = true;
doSaveFile = true; % Should the analysis results be saved?
savePlot = true;
plotPath = 'D:\samuelgb\Documents\MATLAB\Acute Control\Figures\Corrupted_SD_Curves\';
saveFilePath = 'D:\samuelgb\Documents\MATLAB\Acute Control\Scripts\Analysis\Test\Control with Corruption\corruptedSpikeSorting_updatedResults.mat';
numCorruptions = 10; %1000

FPTableFile = 'logisticFPTable.mat'; % The name of the file with the FP lookup table in it

GTOpt  = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'none'); % Optimization options;
GTLB = [0 0];

% Neuron A
alphaA = .7; %2
% betaA = .35; %1
betaA = .035; %1
sigmaA = .05;
thetaA = [alphaA betaA sigmaA];
GTStart1 = [2 .5]; % Strength duration pair to start with

% Neuron C (will test multiple values of each parameter
alphaC = [alphaA alphaA alphaA];
betaC = [betaA .045 .025]; % Previous (Changed 1/20/20): [betaA .5 .1]
% betaC = [.07 .045 .025]; % TODO EMERGENCY POSTER STUFF
sigmaC = [.05 .05 .05];
thetaC = [alphaC' betaC' sigmaC']; % numCs x 3
numCs = length(alphaC);
% numCs = 1; % TODO EMERGENCY POSTER STUFF

% Old
% alphaC = .15;
% betaC = .1;
% sigmaC = .05;
% thetaC = [alphaC betaC sigmaC];

% Neuron B Prime
alphaB = .05; %.3
% betaB = .08; %.4
betaB = .02; %.4
sigmaB = .05;
thetaB = [alphaB betaB sigmaB];
GTStart2 = [.5 5]; % Strength duration pair to start with

% Unified model parameters
unifiedPAExc = .1;
unifiedCInd = 1;
unifiedPIncC = .2;
unifiedHashInd = 2;
unifiedPIncHashA = .3;
unifiedPIncHashC = .3;

% FP Controls
doNeuronA = false;
doNeuronC = true;
doNeuronB = false;

% SD curve firing rate
sdCurveProb = .5; % Would like to find the SD curve corresponding to this firing rate

% SD curve fitting durations
% allTestedDurs = [1 2 5 10 15]; % The durations that are used in the experimen
durLow = 0;
durHigh = 15; % The maximum duration tested
strLow = 0;
strHigh = 10;

% Resolution of the simulation
numV = 301;%1001; % Default is 301
numT = 10001;%20001; % Default is 5001

% Range/resolution of the table
numStrs = 5001;
allStrsMesh = linspace(strLow, strHigh, numStrs);
%     numDurs = 100;
numDurs = numT;
allDursMesh = linspace(durLow, durHigh, numDurs);

% Parameters of SD curve to calculate
numSDCurveDurs = 100;
% allMainSDDurs = [1 2 5 10 15];
allMainSDDurs = linspace(0, 15, numSDCurveDurs);
numMainSDDurs = length(allMainSDDurs);
diffCheckDist = 1e-2;
probDiffFromSD = .01; % The probability calculated for each corrupted SD curve must be at most this far away from the target SD curve value
str0Default = 30;

% Parameters for using automatic beta calculation
% desiredRFPerBeta = [.5 .7 .3]; % The desired response fractions for neuron C at a collection of SD points (will use the pure SD curve points)
desiredRFPerBeta = [.5 .85 .15]; % The desired response fractions for neuron C at a collection of SD points (will use the pure SD curve points)
numBetasToCalc = length(desiredRFPerBeta);

% Define a logistic function that we can control by simply defining where
% "zero" and "one" are (to within some delta difference).
logDelta = 1e-4;
logisticDefByHandles = @(x, h1, h2, delta) logistic(x, log((1/delta) - 1)/((h2 - h1*(log((1/(1 - delta)) - 1)/log((1/delta) - 1)))/(1 - (log((1/(1 - delta)) - 1)/log((1/delta) - 1))) - h1), (h2 - h1*(log((1/(1 - delta)) - 1)/log((1/delta) - 1)))/(1 - (log((1/(1 - delta)) - 1)/log((1/delta) - 1))));
% pOccurrenceNew = @(G, T, kOccur, eHalf, ePMin, ePMax, kStr) ePMin + (abs(ePMax - ePMin)./(1 + exp(-kOccur*((G).*(T) - eHalf) - kStr*G)));
pOccurrence = @(G, T, kOccur, eHalf, ePMin, ePMax) ePMin + (abs(ePMax - ePMin)./(1 + exp(-kOccur*((G).*(T) - eHalf))));
pOccurrenceLogHandles = @(G, T, handle0, handle1) logisticDefByHandles(G.*T, handle0, handle1, logDelta);

% numTestedDurs = length(allTestedDurs); % The number of durations

vThreshLims = [.25 .5]; % Common voltage threshold (scaled such that a voltage of the lowest value will have probability 0 of a spike, and the largest value will have probability 1 of a spike

% Calculation parameters
numParams = 2; % G and T
numNeurons = 2; % Neuron A (in combination with Neuron C), and Neuron B
numPlottedCorruptions = 10;

if useMultistart
   ms = MultiStart; 
end


if useFPInterpolation
    %     sdIndsToUse = 1:100:numDurs;
    %     numSDInds = length(sdIndsToUse);
    
    % Get the true SD curves for each pure neuron
    %     sdCurveAStrengthsPure = zeros(1, numSDInds);
    %     sdCurveCStrengthsPure = zeros(numCs, numSDInds);
    %     sdCurveBStrengthsPure = zeros(1, numSDInds);
    
    %     sdCurveAStrengthsPure = zeros(1, numMainSDDurs);
    %     sdCurveCStrengthsPure = zeros(numCs, numMainSDDurs);
    %     sdCurveBStrengthsPure = zeros(1, numMainSDDurs);
    
    %% Prepare the FokkerPlanck dataset
    % Load the Fokker-Planck data
    %     load(FPTableFile);
    prepareFPWorkspace; % Prepare the FP interpolation equations using irradiance values
    GTUB = [max(allStrs) max(allDurs)];
    
    % Prepare the approximation function
%     approxFun = @(GT, t) lininterpn(allStrs, allDurs, allAlphas, allSigmas, PAll, GT(:,1)*(t(2)/beta), GT(:,2), t(1), t(3));
    fErr = @FPErrorCalc;
    
    gPts = linspace(0, GTUB(1), 5); %[5 20 50];
    tPts = linspace(0, GTUB(2), 5); %[2 7 15];
    %     gPts = gPts;
    tPts = tPts(2:end);
    [gX, tX] = ndgrid(gPts, tPts);
    gtX = [gX(:) tX(:)];
    
    % Cost function
    costFunLambda = 1e-5; %.01;
    costFun = @(P1, P2, G) -(P1 .* (1 - P2)) + costFunLambda*(G(:, 1).^2); % Cost function to optimize using control
    
    %% Find the pure SD curves for each neuron
    % Create a function to return all of the strengths form a probability table
    % for a given duration (finding the duration that in the table that is the
    % closest match).
    
    sdCurveAStrengthsPure = calculateSDStrs(@(GT, thetaInd) approxFun(GT, thetaA), sdCurveProb, allMainSDDurs, str0Default, 1);
    sdCurveBStrengthsPure = calculateSDStrs(@(GT, thetaInd) approxFun(GT, thetaB), sdCurveProb, allMainSDDurs, str0Default, 1);
    %     sdCurveCStrengthsPure = zeros(numBetasToCalc, numMainSDDurs);
    
    % Next, calculate the best beta for neuron C
    if useAutoBetaOpt
        allBetas = zeros(1, numBetasToCalc);
        allFVals = zeros(1, numBetasToCalc); % Just for testing
        
        for thetaInd = 1:numBetasToCalc
            thisIdealRF = desiredRFPerBeta(thetaInd);
            
            if thisIdealRF == sdCurveProb && all([(thetaA(1) == thetaC(thetaInd, 1)) (thetaA(3) == thetaC(thetaInd, 3))])
                % If we want a neuron that will behave the same as neuron
                % A, then just use the same beta value
                thisBeta = thetaA(2);
                fVal = 0;
                
            else
                % Using the SD curve for neuron A, find the beta value that
                % generates a neuron C response closest to the desired RF on
                % each point.
                [thisBeta, fVal] = fminunc(@(beta) sum(approxFun([sdCurveAStrengthsPure' allMainSDDurs'], [thetaC(thetaInd, 1) beta thetaC(thetaInd, 3)]) - thisIdealRF).^2, thetaC(thetaInd, 2));
            end
            
            allBetas(thetaInd) = thisBeta;
            allFVals(thetaInd) = fVal;
            
            if betaPlotTest
                % For each SD pair in the original neuron A SD curve, calculate the
                % RF of each
                rfsAtThisBeta = approxFun([sdCurveAStrengthsPure' allMainSDDurs'], [thetaC(thetaInd, 1) thisBeta thetaC(thetaInd, 3)])' - thisIdealRF;
                figure;
                patch([allMainSDDurs nan], [sdCurveAStrengthsPure nan], [rfsAtThisBeta nan], [rfsAtThisBeta nan], 'Edgecolor', 'interp', 'Facecolor', 'none', 'linewidth', 4);
                xlabel('Duration (ms)');
                ylabel('Strength (mW/mm^s)');
                zlabel(sprintf('Response Fraction Error (from %0.1f)', thisIdealRF));
                title(sprintf('Response Fraction Error with Fitted Beta %d for Neuron C', thetaInd));
                caxis([-.05 .05]);
                colorbar;
            end
        end
        
        % Replace the theta values with these betas
        thetaC(:, 2) = allBetas;
    end
    
    % Calculate the SD curves from the new betas
    sdCurveCStrengthsPure = calculateSDStrs(@(GT, thetaInd) approxFun(GT, thetaC(thetaInd, :)), sdCurveProb, allMainSDDurs, str0Default, numCs); % First calculate the SD curves
    
    % Find the powers of stimulations on the SD curve
    allSDCurveAPows = sdCurveAStrengthsPure.*allMainSDDurs;
    allSDCurveBPows = sdCurveBStrengthsPure.*allMainSDDurs;
    allSDCurveCDurs = repmat(allMainSDDurs, numBetasToCalc, 1);
    allSDCurveCPows = sdCurveCStrengthsPure.*allSDCurveCDurs;
    
    allSDCurveAMeanPows = nanmean(allSDCurveAPows);
    allSDCurveBMeanPows = nanmean(allSDCurveBPows);
    allSDCurveCMeanPows = nanmean(allSDCurveCPows, 2);
    
    % Using the SD curve for neuron C, find the power that should produce
    % enough hash that an A spike will be blocked half of the time (when
    % not considering duration effects).
    % ** Will try using only stims below 5ms **
    betaDurInds = allMainSDDurs <= 5;
    betaDurs = allMainSDDurs(betaDurInds);
    hashHalfProbPows = nanmean(allSDCurveCPows(:, betaDurInds), 2);
    
    % Find the k_occurrence parameter for each neuron C, as well as the min
    % and max probabilities (just so there are more DOF's for the optimizer
    % to play with)
    betaStrs = linspace(0, 2*max(sdCurveAStrengthsPure), length(betaDurs));
    [betaStrMesh, betaDurMesh] = meshgrid(betaStrs, betaDurs);
    betaStrMeshCol = betaStrMesh(:);betaDurMeshCol = betaDurMesh(:);
    
    
    kOccurrence = zeros(numBetasToCalc, 1);
    fValKOccur = zeros(numBetasToCalc, 1);
    for betaInd = 1:numBetasToCalc
        neuronCProbs = approxFun([betaStrMeshCol betaDurMeshCol], thetaC(betaInd, :));
        useInds = neuronCProbs < 1; % Only use strength duration pairs that produce probabilities that are NOT exactly 1 (to capture more useful dynamics)

        % Was having an issue with the below line initially, because the
        % squared error calculation was malformed.  Now, with proper
        % squared error, we are getting a good fit
        %                 [pOccurrenceParams(betaInd, :), fValKOccur(betaInd)] = fminunc(@(k) sum((neuronCProbs(useInds) - pOccurrenceNew(betaStrMeshCol(useInds), betaDurMeshCol(useInds), k(1), hashHalfProbPows(betaInd), k(2), k(3), 0)).^2), [1 0 1]);
        
        % Technically, we may want to be using a log-likelihood (as in
        % FPErrorCalc.m) to find the maximum likelihood, however SQE is
        % giving good results, and MLE will be troublesome to implement
        % (will likely have to remove 0's and 1' from data).
        [kOccurrence(betaInd, :), fValKOccur(betaInd)] = fminunc(@(k) sum((neuronCProbs(useInds) - pOccurrence(betaStrMeshCol(useInds), betaDurMeshCol(useInds), k(1), hashHalfProbPows(betaInd), 0, 1)).^2), 1);
        
        %         [handle1, fValKOccur(betaInd)] = fminunc(@(handle1) sum(neuronCProbs(useInds) - pOccurrenceLogHandles(betaStrMeshCol(useInds), betaDurMeshCol(useInds), 0, handle1).^2), 2);
        %         pOccurrenceParams(betaInd, :) = [];
    end
    
    %%
    if betaPlotTest
        % Plot the current settings for the hash compared to its
        % corresponding neuron C
        testPlotStrs = linspace(0, 2*max(sdCurveAStrengthsPure), length(allMainSDDurs) + 1);
        [testPlotDurMesh, testPlotStrMesh] = meshgrid(allMainSDDurs, testPlotStrs);
        testPlotStrMeshCol = testPlotStrMesh(:);testPlotDurMeshCol = testPlotDurMesh(:);
        neuronCZ = zeros(size(testPlotStrMesh));
        hashZ = zeros(size(testPlotStrMesh));
        for betaInd = 1:numBetasToCalc
            neuronCZCol = approxFun([testPlotStrMeshCol testPlotDurMeshCol], thetaC(betaInd, :));
            %            hashZCol = pOccurrenceNew(testPlotStrMeshCol, testPlotDurMeshCol, pOccurrenceParams(betaInd, 1), hashHalfProbPows(betaInd), pOccurrenceParams(betaInd, 2), pOccurrenceParams(betaInd, 3), 0);
            hashZCol = pOccurrence(testPlotStrMeshCol, testPlotDurMeshCol, kOccurrence(betaInd, 1), hashHalfProbPows(betaInd), 0, 1);
            
            neuronCZ(:) = neuronCZCol;
            hashZ(:) = hashZCol;
            
            figure;
            a = axes;
            hold(a, 'on');
            surf(a, allMainSDDurs, testPlotStrs, neuronCZ, 'FaceColor', [0 0 1], 'Edgealpha', .2);
            surf(a, allMainSDDurs, testPlotStrs, hashZ, 'FaceColor', [1 0 0], 'Edgealpha', .2);
            xlabel('Duration (ms)');
            ylabel('Strength (mw/mm^2)');
            zlabel('Probability');
            title(sprintf('Hash (Red) Approximation of Fitted Neuron C%d (Blue)', betaInd));
        end
        
        
        % For each SD pair in the original neuron A SD curve, calculate the
        % RF of each
        maxPow = max([allSDCurveAPows(:); allSDCurveBPows(:); allSDCurveCPows(:)]);
        maxStr = max([sdCurveAStrengthsPure(:); sdCurveBStrengthsPure(:); sdCurveCStrengthsPure(:)]);
        
        % Neuron A
        figure;
        patch([allMainSDDurs nan], [sdCurveAStrengthsPure nan], [allSDCurveAPows nan], [allSDCurveAPows nan], 'Edgecolor', 'interp', 'Facecolor', 'none', 'linewidth', 4);
        xlabel('Duration (ms)');
        ylabel('Strength (mW/mm^s)');
        zlabel('Stimulation Energy (mWms/mm^2)');
        title('Stimulation Energy of the Neuron A SD Curve');
        ylim([0 maxStr]);
        caxis([0 maxPow]);
        colorbar;
        
        % Neuron B
        figure;
        patch([allMainSDDurs nan], [sdCurveBStrengthsPure nan], [allSDCurveBPows nan], [allSDCurveBPows nan], 'Edgecolor', 'interp', 'Facecolor', 'none', 'linewidth', 4);
        xlabel('Duration (ms)');
        ylabel('Strength (mW/mm^s)');
        zlabel('Stimulation Energy (mWms/mm^2)');
        title('Stimulation Energy of the Neuron B SD Curve');
        ylim([0 maxStr]);
        caxis([0 maxPow]);
        colorbar;
        
        % Neuron C
        fillNan = nan(numBetasToCalc, 1);
        figure;
        patch([allSDCurveCDurs fillNan]', [sdCurveCStrengthsPure fillNan]', [allSDCurveCPows fillNan]', [allSDCurveCPows fillNan]', 'Edgecolor', 'interp', 'Facecolor', 'none', 'linewidth', 4);
        xlabel('Duration (ms)');
        ylabel('Strength (mW/mm^s)');
        zlabel('Stimulation Energy (mWms/mm^2)');
        title('Stimulation Energy of Each Neuron C SD Curve');
        ylim([0 maxStr]);
        caxis([0 maxPow]);
        colorbar;
    end
    
    %%
    % Corruption options
    corruptionTypes = {'Excluded Spikes', 'Added Spikes', 'Deleted Spikes'}; %, 'Unified Corruption'};
%     corruptionTypes = {'Added Spikes'}; % TODO EMERGENCY POSTER STUFF
    numCorruptionTypes = length(corruptionTypes);
    hashFcn = @(GT, corrInstanceInd) pOccurrence(GT(:, 1), GT(:, 2), kOccurrence(corrInstanceInd), hashHalfProbPows(corrInstanceInd), 0, 1);
    
    % Define the functions used to calculate the final observed PA firing
    % probability, qhich receives as input only the strength/duration that
    % is being tested, and the value of corruption (between 0 and 1)
    exclusionCorruptionFcn = @(GT, corrVal, corrInstanceInd) approxFun(GT, thetaA)*(1 - corrVal);
    additionCorruptionFcn = @(GT, corrVal, corrInstanceInd) ProbUnion(approxFun(GT, thetaA), corrVal*approxFun(GT, thetaC(corrInstanceInd, :)));
    deletionCorruptionFcn = @(GT, corrVal, corrInstanceInd) approxFun(GT, thetaA).*(1 - corrVal*hashFcn(GT, corrInstanceInd));
    unifiedCorruptionFcn = @(GT, corrVal, corrInstanceInd) ProbUnion(approxFun(GT, thetaA)*(1 - (unifiedPAExc*(corrInstanceInd ~= 1) + corrVal*(corrInstanceInd == 1)))*(1 - (unifiedPIncHashA*(corrInstanceInd ~= 3) + corrVal*(corrInstanceInd == 3))*hashFcn(GT, unifiedHashInd)), approxFun(GT, thetaC(unifiedCInd, :))*(unifiedPIncC*(corrInstanceInd ~= 2) + corrVal*(corrInstanceInd == 2))*(1 - (unifiedPIncHashC*(corrInstanceInd ~= 3) + corrVal*(corrInstanceInd == 3))*hashFcn(GT, unifiedHashInd)));
    
    %     exclusionCorruptionFcn = @(PA, PC, C) PA*(1 - C); % Neuron A spikes AND NOT lambda [A spike excluded from cluster]
    %     additionCorruptionFcn = @(PA, PC, C) ProbUnion(PA, C*PC); % Neuron A spikes OR (A' spikes AND lambda [A' spike is mischaracterized as A])
    %     deletionCorruptionFcn = @(PA, PC, C) PA*(1 - C*PC); % Neuron A spikes AND NOT (A' spikes AND lambda [A' spike interferes with A spike])
    %     pPow = @(G, T, kPow, powHalf, powPMin, powPMax) powPMin + (abs(powPMax - powPMin)./(1 + exp(-kPow*(G.*T - powHalf))));
    %     pDur = @(T, kDur, durHalf, durPMin, durPMax) durPMin + (abs(durPMax - durPMin)./(1 + exp(-kDur*(-T + durHalf))));
    %     hashDeletionCorruptionFcn = @(PA, G, T, pDur, pPow) PA .* (1 - pDur(T) .* pPow(G, T));
    
    allCorruptionFcns = {...
        exclusionCorruptionFcn;...
        additionCorruptionFcn;...
        deletionCorruptionFcn;...
        unifiedCorruptionFcn};
    
%     allCorruptionFcns = {additionCorruptionFcn}; % TODO EMERGENCY POSTER STUFF
    
    allCorruptions = repmat(linspace(0, 1, numCorruptions), numCorruptionTypes, 1);
%     allCorruptions(1, :) = linspace(0, .49, numCorruptions);
    tmpCorrVals = .5 - logspace(-4, log10(.5), numCorruptions);
    allCorruptions(1, :) = tmpCorrVals(end:-1:1);
    
%     allCalcDoubleSD = [false false true]; % Should two values of the SD curve be calculated (needed when the SD curve is non-monotonic, as in the deleted spikes case)
    allCalcDoubleSD = [false false true true]; % Should two values of the SD curve be calculated (needed when the SD curve is non-monotonic, as in the deleted spikes case)
    allCalcUseCorruptingSignal = [false true true true]; % Is a second signal required for the corruption (i.e. can we use only a single "C" value, or should we go through all numCs worth?)
    
    quietOpt = optimoptions('fmincon', 'Display', 'none', 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-16, 'FunctionTolerance', 1e-12);
    
    fprintf('Done calculating raw SD curve.\n\n');
else
    % Corruption options
    corruptionTypes = {'Linear Interpolation', 'Added Spikes', 'Deleted Spikes'};
    numCorruptionTypes = length(corruptionTypes);
    
    linCorruptions = [linspace(1,0,numCorruptions)' linspace(0,1,numCorruptions)'];
    addedSpikesCorruptions = [ones(numCorruptions,1) linspace(0,1,numCorruptions)'];
    deletedSpikesCorruptions = [linspace(1,0,numCorruptions)' zeros(numCorruptions,1)]; % Removed 10/17/19, to update deleted spikes such that it is affected by the activity of the corrupting neuron
    % deletedSpikesCorruptions = [ones(numCorruptions,1) linspace(0,-1,numCorruptions)'];
    allCorruptions = cat(3,cat(3,linCorruptions,addedSpikesCorruptions),deletedSpikesCorruptions); % numCorruptions x neuronNum (alpha or alpha-prime) x corruptionType
    %% Calculate the Fokker-Planck spiking probabilities for the two neurons
    if ~(exist('PAllA', 'var') || exist('ProbabilityTablesCorruptionTest.mat', 'file'))
        % If this data has not yet been computed (either actively in the
        % workspace or in a .mat file), then recalculate it (takes a LONG time,
        % 10-30 minutes).
        
        tic;
        
        PAllA = zeros(numStrs, numDurs);
        PAllC = zeros(numStrs, numDurs);
        PAllB = zeros(numStrs, numDurs);
        fprintf('Starting table calculation...\n');
        for strInd = 1:numStrs
            G = allStrsMesh(strInd);
            % Inform user
            if mod(strInd, round(numStrs/10)) == 0
                fprintf('Start strength %d out of %d...\n', strInd, numStrs);
            end
            
            if doNeuronA
                [~, PAllA(strInd,:)] = FokkerPlanckProbSim(G, durHigh, alphaA, betaA, vThreshLims(2), false, [numV numT]);
            end
            
            if doNeuronC
                [~, PAllC(strInd,:)] = FokkerPlanckProbSim(G, durHigh, alphaC(2), betaC(2), vThreshLims(2), false, [numV numT]);
            end
            
            if doNeuronB
                [~, PAllB(strInd,:)] = FokkerPlanckProbSim(G, durHigh, alphaB, betaB, vThreshLims(2), false, [numV numT]);
            end
        end
        fprintf('\nTable calculation took %0.2f seconds to calculate %d strengths (numV = %d, numT = %d).\n', toc, numStrs, numV, numT);
        
        if doPlot
            if doNeuronA
                h1 = figure;
                mesh(allDursMesh, allStrsMesh, PAllA);
                xlabel('Duration');
                ylabel('Strength');
                zlabel('Probability');
                title('Probability of Spiking (Neuron A)');
            end
            
            if doNeuronC
                h2 = figure;
                mesh(allDursMesh, allStrsMesh, PAllC);
                xlabel('Duration');
                ylabel('Strength');
                zlabel('Probability');
                title('Probability of Spiking (Neuron A'')');
            end
            
            if doNeuronB
                h3 = figure;
                mesh(allDursMesh, allStrsMesh, PAllB);
                xlabel('Duration');
                ylabel('Strength');
                zlabel('Probability');
                title('Probability of Spiking (Neuron B)');
            end
        end
        
        if all([doNeuronA doNeuronC doNeuronB])
            save('ProbabilityTablesCorruptionTest.mat', 'PAllA', 'PAllC', 'PAllB', 'allStrsMesh', 'allDursMesh');
        end
        
        if doNeuronA
            save('ProbabilityTablesCorruptionTestA.mat', 'PAllA', 'allStrsMesh', 'allDursMesh');
        end
        
        if doNeuronC
            save('ProbabilityTablesCorruptionTestC.mat', 'PAllC', 'allStrsMesh', 'allDursMesh');
        end
        
        if doNeuronB
            save('ProbabilityTablesCorruptionTestB.mat', 'PAllB', 'allStrsMesh', 'allDursMesh');
        end
    elseif ~exist('PDiffA', 'var') && exist('ProbabilityTablesCorruptionTest.mat', 'file')
        load('ProbabilityTablesCorruptionTest.mat');
    end
    
    %     PTableDurationLookup = @(PTable, dur) PTable(:, abs(allDursMesh - dur) == min(abs(allDursMesh - dur)));
    %     PTableLookup = @(PTable, str, dur) PTable(find(abs(allStrsMesh - str) == min(abs(allStrsMesh - str)), 1), find(abs(allDursMesh - dur) == min(abs(allDursMesh - dur)), 1));
    
    %% Find the pure SD curves for each neuron
    % Create a function to return all of the strengths form a probability table
    % for a given duration (finding the duration that in the table that is the
    % closest match).
    
    sdCurveAStrengthLogicals = min(abs(PAllA - sdCurveProb), [], 1) == abs(PAllA - sdCurveProb);
    sdCurveCStrengthLogicals = min(abs(PAllC - sdCurveProb), [], 1) == abs(PAllC - sdCurveProb);
    sdCurveBStrengthLogicals = min(abs(PAllB - sdCurveProb), [], 1) == abs(PAllB - sdCurveProb);
    sdCurveAStrengthInds = zeros(size(sdCurveAStrengthLogicals,2),1);
    sdCurveCStrengthInds = zeros(size(sdCurveCStrengthLogicals,2),1);
    sdCurveBStrengthInds = zeros(size(sdCurveBStrengthLogicals,2),1);
    
    for i = 1:size(sdCurveAStrengthLogicals,2)
        sdCurveAStrengthInds(i) = find(sdCurveAStrengthLogicals(:,i), 1, 'last');
    end
    for i = 1:size(sdCurveCStrengthLogicals,2)
        sdCurveCStrengthInds(i) = find(sdCurveCStrengthLogicals(:,i), 1, 'last');
    end
    for i = 1:size(sdCurveBStrengthLogicals,2)
        sdCurveBStrengthInds(i) = find(sdCurveBStrengthLogicals(:,i), 1, 'last');
    end
    
    sdCurveAStrengthsPure = allStrsMesh(sdCurveAStrengthInds);
    sdCurveCStrengthsPure = allStrsMesh(sdCurveCStrengthInds);
    sdCurveBStrengthsPure = allStrsMesh(sdCurveBStrengthInds);
    
    % Fix errors in calculating the max strengths
    sdCurveAStrengthsPure(:,~any(diff(PAllA),1)) = max(allStrsMesh);
    sdCurveCStrengthsPure(:,~any(diff(PAllC),1)) = max(allStrsMesh);
    sdCurveBStrengthsPure(:,~any(diff(PAllB),1)) = max(allStrsMesh);
end

%% Simulate the SD curve calculations at different corruption levels
% First, find the SD curve of the unified model
% unifiedCorruptedSDCurve = calculateSDStrs(@(GT, thetaInd) unifiedCorruptionFcn(GT), sdCurveProb, allMainSDDurs, str0Default, 1);

% Find the sd curves for A (observed) and B
% [~, sdCurveBStrengthInds] = min(abs(PAllB - sdCurveProb), [], 1);
% sdCurveBStrengths = allStrsMesh(sdCurveBStrengthInds);
sdCurveAPerceivedStrengths = nan(numCorruptions, numCorruptionTypes, numCs, numMainSDDurs, 1 + useFPInterpolation);
optimalGTs = zeros(numNeurons, numCorruptions, numCorruptionTypes, numCs, numParams);
allStartStrs = linspace(.05, 1.5, 5)'; % Create an array of strengths to initialize optimization with

for corrInstanceNum = 1:numCs
    thisThetaC = thetaC(corrInstanceNum, :);
    
    for corruptionNum = 1:numCorruptions
        % Go through each amount of corruption
        if mod(corruptionNum, round(numCorruptions/10)) == 0
            fprintf('Starting corruption number %d out of %d...\n\n', corruptionNum, numCorruptions);
        end
        
        if useFPInterpolation
            corruption = allCorruptions(corrInstanceNum, corruptionNum);
        end
        
        % Hold onto the previously optimal GT, so that it's one possible start
        % point for the next round (smooths out the data)
        thisStartGT1 = ones(numCorruptionTypes, numParams);
        thisStartGT2 = ones(numCorruptionTypes, numParams);
        
        for corruptionTypeNum = 1:numCorruptionTypes
            % Get the corruption level for this series of measurements
            
            if useFPInterpolation
                if ~allCalcUseCorruptingSignal(corruptionTypeNum) && corrInstanceNum ~= 1
                   % If we only need to run this corruption type once
                   % (because it doesn't require all numCs external
                   % corruption signals) and we have already completed them
                   % (corruptionNum is greater than 1), then skip this
                   % round
                   continue;
                end
                thisCorFcn = allCorruptionFcns{corruptionTypeNum};
                thisDoCalcDoubleSD = allCalcDoubleSD(corruptionTypeNum);
            else
                corruption = allCorruptions(corruptionNum, :, corruptionTypeNum);
                
                if ~(useFPInterpolation && useFPCostOptimization)
                    % Find the probability of firing for both neurons A and A-prime
                    % individually, given this level of corruption
                    PA = PAllA*corruption(1);
                    PC = PAllC*corruption(2);
                    % Adjust corruption to real-world clustering errors. Over/under
                    % estimates (from adding sin-waves)
                    % 1 over alpha should be around order 10 (5ms - 20ms)
                    
                    % Calculate the probability that EITHER neuron will fire at a given
                    % strength
                    PAllAEither = ProbUnion(PA, PC);
                end
            end
            
            %         % Test each duration to find the 50% strength level
            %         allTestedStrengths = zeros(size(allTestedDurs));
            %         for durNum = 1:numTestedDurs
            %             % Get the duration for this measurement
            %             duration = allTestedDurs(durNum);
            %
            %             % Find the probability of either neurons A or A-Prime firing
            %             % for all strength levels given this duration
            %             PAEither = PTableDurationLookup(PAllAEither, duration);
            %
            %             % Find the index in PEither that is closest to ideal SD curve line
            %             % for the chosen firing probability
            %             [~, strengthInd] = min(abs(PAEither - sdCurveProb));
            %
            %             % Use the index to get the strength
            %             allTestedStrengths(durNum) = allStrsMesh(strengthInd);
            %         end
            
            %                 sdCurveAStrengthLogicals = min(abs(PAllAEither - sdCurveProb), [], 1) == abs(PAllAEither - sdCurveProb);
            %                 sdCurveAStrengthInds = zeros(size(sdCurveBStrengthLogicals,2),1);
            %                 for i = 1:size(sdCurveAStrengthLogicals,2)
            %                     sdCurveAStrengthInds(i) = find(sdCurveAStrengthLogicals(:,i), 1, 'last');
            %                 end
            
            % Find the strengths that define the strength duration curve for
            % the observed (corrupted) neuron A
            if useFPInterpolation
                str0 = str0Default;
                %                 for durIndInd = 1:numSDInds
                for durIndInd = 1:numMainSDDurs
                    % Get this duration
                    %                     thisDur = allDursMesh(sdIndsToUse(durIndInd));
                    thisDur = allMainSDDurs(durIndInd);
                    
                    % Find the initial str0 (try to make it close to
                    % previous SD's that were already calculated)
                    lastStr = sdCurveAPerceivedStrengths(max(1, corruptionNum - 1), corruptionTypeNum, corrInstanceNum, durIndInd, 1);
%                     if ~isnan(lastStr)
%                         % If the strength from the last calculated SD is
%                         % not nan, and a gradient is defined, then use it.
%                         % Otherwise, continue using the last str0
%                         testProb1 = thisCorFcn([lastStr thisDur], corruption, corrInstanceNum);
%                         testProb2 = thisCorFcn([(lastStr + diffCheckDist) thisDur], corruption, corrInstanceNum);
%                         
%                         if testProb1 ~= testProb2
%                             str0 = lastStr;
%                         end
%                     end
                    
                    % First, check if our strength is on the correct order
                    for i = 1:10
%                         testProb1 = thisCorFcn(approxFun([str0 thisDur], thetaA), approxFun([str0 thisDur], thisThetaC), corruption);
%                         testProb2 = thisCorFcn(approxFun([(str0 + diffCheckDist) thisDur], thetaA), approxFun([(str0 + diffCheckDist) thisDur], thisThetaC), corruption);
                        
                        testProb1 = thisCorFcn([str0 thisDur], corruption, corrInstanceNum);
                        testProb2 = thisCorFcn([(str0 + diffCheckDist) thisDur], corruption, corrInstanceNum);
                        
                        if testProb1 ~= testProb2
                            % The two values are different enough that a
                            % gradient is defined
                            break;
                        else
                            % Not gradient is defined at this point, try a
                            % smaller strength
                            str0 = str0/2;
                        end
                        
                        %                        if i == 10
                        %                           warning(fprintf('\n*****\n**********\n***************\n\nCould not get a gradient defined at: CInd = %d, corrNum = %d, corrType = %d, durInd = %d\n\n***************\n**********\n*****\n\n',  corrInstanceNum, corruptionNum, corruptionTypeNum, durIndInd));
                        %                        end
                    end
                    
%                     foundDouble = false;
                    if thisDoCalcDoubleSD
                        % If we are expecting a double SD curve, then first
                        % find the strength that gives the maximum firing
                        % probability
                        %                         thisStrMax = fmincon(@(s) -thisCorFcn(approxFun([s thisDur], thetaA), approxFun([s thisDur], thisThetaC), corruption), str0, [], [], [], [], 0, GTUB(1), [], quietOpt);
                        if useMultistart
                            probMax = createOptimProblem('fmincon', 'objective', @(s)-thisCorFcn([s thisDur], corruption, corrInstanceNum), 'x0', str0, 'lb', 0, 'ub', GTUB(1));
                            thisStartStrs = CustomStartPointSet([allStartStrs;str0;lastStr]);
                            thisStrMax = run(ms, probMax, thisStartStrs);
                        else
                            thisStrMax = fmincon(@(s) -thisCorFcn([s thisDur], corruption, corrInstanceNum), str0, [], [], [], [], 0, GTUB(1), [], quietOpt);
                        end
                        
                        if isnan(thisStrMax)
                            thisStrMax = GTUB(1);
                        end
                        
                        % Find the sdCurveLevel BELOW this maximum firing
                        % probability
                        %                         thisSBLow = fmincon(@(s)(thisCorFcn(approxFun([s thisDur], thetaA), approxFun([s thisDur], thisThetaC), corruption) - sdCurveProb)^2, thisStrMax/2, [], [], [], [], 0, thisStrMax, [], quietOpt);
                        %                         thisPLow = thisCorFcn(approxFun([thisSBLow thisDur], thetaA), approxFun([thisSBLow thisDur], thisThetaC), corruption);
                        if useMultistart
                            probLow = createOptimProblem('fmincon', 'objective', @(s) (thisCorFcn([s thisDur], corruption, corrInstanceNum) - sdCurveProb)^2, 'x0', thisStrMax/2, 'lb', 0, 'ub', thisStrMax);
                            thisStartStrs = CustomStartPointSet([linspace(.01, thisStrMax, 5)';str0;lastStr]);
                            thisSBLow = run(ms, probLow, thisStartStrs);
                        else
                            thisSBLow = fmincon(@(s)(thisCorFcn([s thisDur], corruption, corrInstanceNum) - sdCurveProb)^2, thisStrMax/2, [], [], [], [], 0, thisStrMax, [], quietOpt);
                        end
                        
                        thisPLow = thisCorFcn([thisSBLow thisDur], corruption, corrInstanceNum);
                        if abs(thisPLow - sdCurveProb) > probDiffFromSD
                            % If the strength calculated doesn't REALLY
                            % give the SD curve probability level, then
                            % throw it out (we have a minimum, not a zero)
                            thisSBLow = nan;
                        end
                        
                        % Find the sdCurveLevel ABOVE this maximum firing
                        % probability
                        if thisStrMax < GTUB(1)*.95
                            % Assuming
                            %                             thisSBHigh = fmincon(@(s)(thisCorFcn(approxFun([s thisDur], thetaA), approxFun([s thisDur], thisThetaC), corruption) - sdCurveProb)^2, thisStrMax*2, [], [], [], [], thisStrMax, GTUB(1), [], quietOpt);
                            %                             thisPHigh = thisCorFcn(approxFun([thisSBHigh thisDur], thetaA), approxFun([thisSBHigh thisDur], thisThetaC), corruption);
                            
                            if useMultistart
                                probHigh = createOptimProblem('fmincon', 'objective', @(s) (thisCorFcn([s thisDur], corruption, corrInstanceNum) - sdCurveProb)^2, 'x0', thisStrMax*2, 'lb', thisStrMax, 'ub', GTUB(1));
                                thisStartStrs = CustomStartPointSet([linspace(thisStrMax, GTUB(1), 5)';str0;lastStr]);
                                thisSBHigh = run(ms, probHigh, thisStartStrs);
                            else
                                thisSBHigh = fmincon(@(s)(thisCorFcn([s thisDur], corruption, corrInstanceNum) - sdCurveProb)^2, thisStrMax*2, [], [], [], [], thisStrMax, GTUB(1), [], quietOpt);
                            end
                            
                            thisPHigh = thisCorFcn([thisSBHigh thisDur], corruption, corrInstanceNum);
                            if abs(thisPHigh - sdCurveProb) > probDiffFromSD
                                % If the strength calculated doesn't REALLY
                                % give the SD curve probability level, then
                                % throw it out (we have a minimum, not a zero)
                                thisSBHigh = nan;
                            end
                        else
                            thisSBHigh = nan;
                        end
                        
                        thisSB = [thisSBLow thisSBHigh];
                        foundDouble = any(~isnan(thisSB));
                        %                     end
                    else
                        
                        %                     if ~thisDoCalcDoubleSD || (thisDoCalcDoubleSD && ~foundDouble)
                        % If we are only looking for a single SD curve (or
                        % we're looking for a double, but it failed for
                        % some reason)
                        try
                            % Find the strength that correponds closest to the SD curve
                            %                 thisSB = fzero(@(s) sum(corruption.*[approxFun([s thisDur], thetaA) approxFun([s thisDur], thetaC)]) - sdCurveProb, str0);
                            %                             thisSB = fmincon(@(s) (thisCorFcn(approxFun([s thisDur], thetaA), approxFun([s thisDur], thisThetaC), corruption) - sdCurveProb)^2, str0, [], [], [], [], 0, GTUB(1), [], quietOpt);
                            %                             thisP = thisCorFcn(approxFun([thisSB thisDur], thetaA), approxFun([thisSB thisDur], thisThetaC), corruption);
                            
                            thisSB = fmincon(@(s) (thisCorFcn([s thisDur], corruption, corrInstanceNum) - sdCurveProb)^2, str0, [], [], [], [], 0, GTUB(1), [], quietOpt);
                            thisP = thisCorFcn([thisSB thisDur], corruption, corrInstanceNum);
                            
                            if abs(thisP - sdCurveProb) > probDiffFromSD
                                % If the strength calculated doesn't REALLY
                                % give the SD curve probability level, then
                                % throw it out (we have a minimum, not a zero)
                                thisSB = nan;
                            end
                        catch
                            thisSB = NaN;
                        end
                    end
                    
                    % Save the value, and keep track of it to quicken the
                    % optimization later
                    if allCalcUseCorruptingSignal(corruptionTypeNum)
                        sdCurveAPerceivedStrengths(corruptionNum, corruptionTypeNum, corrInstanceNum, durIndInd, :) = thisSB;
                    else
                        % If we only are using one "C" for this corruption
                        % type, then save this value to "all Cs"
                        sdCurveAPerceivedStrengths(corruptionNum, corruptionTypeNum, :, durIndInd, :) = thisSB;
                    end
                    if isnan(any(thisSB))
                        str0 = str0Default;
                    else
                        str0 = mean(thisSB);
                    end
                    
                    if isnan(str0)
                       str0 = str0Default; 
                    end
                end
            else
                [allMins, sdCurveAStrengthInds] = min(abs(PAllAEither - sdCurveProb), [], 1);
                
                if min(allMins) < 1e-3
                    sdCurveAPerceivedStrengths(corruptionNum, corruptionTypeNum, corrInstanceNum, :) = allStrsMesh(sdCurveAStrengthInds);
                else
                    sdCurveAPerceivedStrengths(corruptionNum, corruptionTypeNum, corrInstanceNum, :) = nan(length(sdCurveAStrengthInds),1);
                end
                
                % Fix errors in calculating the max strengths
                sdCurveAPerceivedStrengths(:, corruptionTypeNum, ~any(diff(PAllA), corrInstanceNum ,1)) = max(allStrsMesh);
            end
            
            
            % With a complete SD curve, find the ideal GT pair
            % NEED A METHOD TO FIND THE GT PAIRS FROM THE SD CURVE
            % Use an optimization technique to find the optimum GT pair for the
            % pair of neurons A/A-prime and B
            
            if useFPCostOptimization
                % Perform a cost optimization with a corrupted neuron A
                %             fObsA = @(GT) sum(corruption.*[approxFun(GT, thetaA) approxFun(GT, thetaC)]);
%                 fObsA = @(GT) thisCorFcn(approxFun(GT, thetaA), approxFun(GT, thisThetaC), corruption);
                fObsA = @(GT) thisCorFcn(GT, corruption, corrInstanceNum);
                
                fThisCost1 = @(GT) costFun(fObsA(GT), approxFun(GT, thetaB), GT);
                fThisCost2 = @(GT) costFun(approxFun(GT, thetaB), fObsA(GT), GT);
                
                thisGT1 = multipleInitialPoint(@(x) fmincon(@(T) fThisCost1(T), x, [], [], [], [], GTLB, GTUB, [], GTOpt), [gtX; thisStartGT1(corruptionTypeNum, :)], false);
                thisGT2 = multipleInitialPoint(@(x) fmincon(@(T) fThisCost2(T), x, [], [], [], [], GTLB, GTUB, [], GTOpt), [gtX; thisStartGT2(corruptionTypeNum, :)], false);
                
                optimalGTs(1, corruptionNum, corruptionTypeNum, corrInstanceNum, :) = thisGT1;
                optimalGTs(2, corruptionNum, corruptionTypeNum, corrInstanceNum, :) = thisGT2;
                
                thisStartGT1(corruptionTypeNum, :) = thisGT1;
                thisStartGT2(corruptionTypeNum, :) = thisGT2;
            else
                % Perform the optimization
                diff1 = diffFun(PAllAEither, PAllB);
                diff2 = diffFun(PAllB, PAllAEither);
                %     if any(corruptionNum == [4 5])
                %         figure;
                %         mesh(diff1);
                %         title(num2str(corruptionNum));
                %     end
                thisMin1 = min(diff1(:));
                thisMin2 = min(diff2(:));
                diffMinInd1 = find(diff1(:) == thisMin1);
                diffMinInd2 = find(diff2(:) == thisMin2);
                
                if length(diffMinInd1) == 1
                    [GSub1, TSub1] = ind2sub(size(diff1), diffMinInd1);
                    [optT1, optG1] = minimize2DSpline(diff1, allDursMesh, allStrsMesh, TSub1, GSub1);
                    optimalGTs(1, corruptionNum, corruptionTypeNum, corrInstanceNum, :) = permute([optG1 optT1], [5 4 3 2 1]);
                    
                    %             optimalGTs(1, corruptionNum, corruptionTypeNum, :) = permute([allStrsMesh(GSub1); allDursMesh(TSub1)], [4 3 2 1]);
                else
                    optimalGTs(1, corruptionNum, corruptionTypeNum, corrInstanceNum, :) = nan(1,1,1,1,2);
                end
                
                if length(diffMinInd2) == 1
                    [GSub2, TSub2] = ind2sub(size(diff1), diffMinInd2);
                    [optT2, optG2] = minimize2DSpline(diff2, allDursMesh, allStrsMesh, TSub2, GSub2);
                    optimalGTs(2, corruptionNum, corruptionTypeNum, corrInstanceNum, :) = permute([optG2 optT2], [5 4 3 2 1]);
                    
                    %             optimalGTs(2, corruptionNum, corruptionTypeNum, :) = permute([allStrsMesh(GSub2); allDursMesh(TSub2)], [4 3 2 1]);
                else
                    optimalGTs(2, corruptionNum, corruptionTypeNum, corrInstanceNum, :) = nan(1,1,1,1,2);
                end
            end
            
            %     fprintf('corruptionNum = %d\nmin: %0.3f\nloc: %d\n\n', corruptionNum, min(diff1(:)), GSub1);
            %         optimalGTs(1, :, corruptionNum) = fminunc(@(x) minimizationFunction(x, PAllAEither, PAllB),GTStart1,options);
        end
    end
end

% Save the data, if desired
if doSaveFile
    doSaveFile = false;
    save(saveFilePath);
end

%% Display the results
if doPlot
    close all;
    fontSize = 40;
    thinLineWidth = 2;
    normLineWidth = 4;
    boldLineWidth = 6;
    %     roundedGTs = round(optimalGTs, 2);
    tol = 1e-2;
    nearMaxInds = abs(repmat(permute([strHigh;durHigh], [5 4 3 2 1]), size(optimalGTs,1) , size(optimalGTs, 2), size(optimalGTs, 3), size(optimalGTs, 4), 1) - optimalGTs) < tol;
    roundedGTs = optimalGTs;
    roundedGTs(nearMaxInds) = round(optimalGTs(nearMaxInds),2);
    
    if doPlotSecondCurve
        sdDoubleCurveInds = [1 2];
        dblStr = 'Dbl';
    else
        sdDoubleCurveInds = 1;
        dblStr = '';
    end
    
    plotLabels = {'Exc', 'Add', 'Del'};
    
    h1 = zeros(numCorruptionTypes, 1);
    h2 = zeros(numCorruptionTypes, 1);
    h3 = zeros(numCorruptionTypes, 1);
    corruptionsAxis = linspace(0,1,numCorruptions);
    for corruptionTypeNum = 1:numCorruptionTypes
        if allCalcDoubleSD(corruptionTypeNum)
           thisDblStr = dblStr;
        else
            thisDblStr = '';
        end
        
        for corrInstanceNum = 1:numCs
            
            %             % GT results
            %             h1(corruptionTypeNum) = figure;
            %
            %             a1 = subplot(2,1,1);
            %             hold(a1, 'on');
            %
            %             a2 = subplot(2,1,2);
            %             hold(a2, 'on');
            %
            %             h2(corruptionTypeNum) = figure;
            %             % h3 = figure;
            %             a3 = subplot(2,1,1);
            %             hold(a3, 'on');
            %
            %             % h4 = figure;
            %             a4 = subplot(2,1,2);
            %             hold(a4, 'on');
            %
            %             %             plot(a1, corruptionsAxis*100, squeeze(roundedGTs(1,:, corruptionTypeNum, corrInstanceNum,1)), 'b', 'linewidth', normLineWidth);
            % %             plot(a1, allCorruptions*100, squeeze(optimalGTs(2,1,:)));
            %
            %             title(a1, sprintf('Ideal G1 with Corruption Type: %s', corruptionTypes{corruptionTypeNum}))
            %             xlabel(a1, 'Corruption (%)');
            %             ylabel(a1, 'G');
            %             % legend(a1, 'G1', 'G2');
            %
            %             plot(a2, corruptionsAxis*100, squeeze(roundedGTs(1,:, corruptionTypeNum,corrInstanceNum,2)), 'b', 'linewidth', normLineWidth);
            %             % plot(a2, allCorruptions*100, squeeze(optimalGTs(2,2,:)));
            %
            %             title(a2, sprintf('Ideal T1 with Corruption Type: %s', corruptionTypes{corruptionTypeNum}))
            %             xlabel(a2, 'Corruption (%)');
            %             ylabel(a2, 'T');
            %             % legend(a2, 'T1', 'T2');
            %
            %             plot(a3, corruptionsAxis*100, squeeze(roundedGTs(2,:, corruptionTypeNum,corrInstanceNum,1)), 'g', 'linewidth', normLineWidth);
            %             % plot(a3, allCorruptions*100, squeeze(optimalGTs(2,1,:)));
            %
            %             title(a3, sprintf('Ideal G2 with Corruption Type: %s', corruptionTypes{corruptionTypeNum}))
            %             xlabel(a3, 'Corruption (%)');
            %             ylabel(a3, 'G');
            %             % legend(a3, 'G1', 'G2');
            %
            %             plot(a4, corruptionsAxis*100, squeeze(roundedGTs(2,:, corruptionTypeNum,corrInstanceNum,2)), 'g', 'linewidth', normLineWidth);
            %             % plot(a4, allCorruptions*100, squeeze(optimalGTs(2,2,:)));
            %
            %             title(a4, sprintf('Ideal T2 with Corruption Type: %s', corruptionTypes{corruptionTypeNum}))
            %             xlabel(a4, 'Corruption (%)');
            %             ylabel(a4, 'T');
            %             % legend(a4, 'T1', 'T2');
            %
            %             set(findall(h1(corruptionTypeNum),'-property','FontSize'),'FontSize',fontSize)
            %             set(findall(h2(corruptionTypeNum),'-property','FontSize'),'FontSize',fontSize)
            if allCalcUseCorruptingSignal(corruptionTypeNum) || corrInstanceNum == 1
                % SD curves
                h3(corruptionTypeNum) = figure;
                a = axes;
                hold(a, 'on');
                %             c = [zeros(numPlottedCorruptions,2) linspace(0,.8, numPlottedCorruptions)'];
                c = jet(1000);
                
                %             line2 = plot(a, allDursMesh(sdIndsToUse), sdCurveAStrengthsPure, 'k:', 'linewidth', boldLineWidth);
                %             line3 = plot(a, allDursMesh(sdIndsToUse), sdCurveCStrengthsPure(corrInstanceNum, :), 'b', 'linewidth', boldLineWidth);
                %             line4 = plot(a, allDursMesh(sdIndsToUse), sdCurveBStrengthsPure, 'g:', 'linewidth', boldLineWidth);
                
                %             line2 = plot(a, allMainSDDurs, sdCurveAStrengthsPure, 'k:', 'linewidth', boldLineWidth);
                %             line3 = plot(a, allMainSDDurs, sdCurveCStrengthsPure(corrInstanceNum, :), 'b', 'linewidth', boldLineWidth);
                %             line4 = plot(a, allMainSDDurs, sdCurveBStrengthsPure, 'g:', 'linewidth', boldLineWidth);
                
                line1 = zeros(numPlottedCorruptions, 1);
                for corruptionNum = 1:numPlottedCorruptions
                    corruptionInd = floor((corruptionNum/numPlottedCorruptions)*numCorruptions);
                    if any(corruptionNum == [1 numPlottedCorruptions])
                        lineWidth = boldLineWidth;
                        %                 color = [0 0 1];
                    else
                        lineWidth = normLineWidth;
                        %                 color = c(corruptionNum,:);
                    end
                    colorInd = round(allCorruptions(corruptionTypeNum, corruptionNum)*(length(c) - 1) + 1);
                    color = c(colorInd,:);
                    %                 thisPlotLines = plot(a, allDursMesh(sdIndsToUse), squeeze(sdCurveAPerceivedStrengths(corruptionInd, corruptionTypeNum, corrInstanceNum,:, :)), 'color', color, 'linewidth', lineWidth);
                    thisPlotLines = plot(a, allMainSDDurs, squeeze(sdCurveAPerceivedStrengths(corruptionInd, corruptionTypeNum, corrInstanceNum,:, sdDoubleCurveInds)), 'color', color, 'linewidth', lineWidth);
                    line1(corruptionNum) = thisPlotLines(1);
                end
                allStrData = sdCurveAPerceivedStrengths(:, corruptionTypeNum, :, :, sdDoubleCurveInds);
                
                if corruptionTypeNum == 2
                    thisStrData = sdCurveCStrengthsPure(corrInstanceNum, :);
                    plot(a, allMainSDDurs, thisStrData, 'k--', 'linewidth', 2);
                    allStrData = [allStrData(:); thisStrData(:)];
                end
                
                maxStr = max(allStrData(:));
                colormap(a, c);
                colorbar('peer', a);
                xlabel(a, 'Duration ($ms$)', 'Interpreter', 'Latex');
                ylabel(a, 'Strength ($\frac{mW}{mm^2}$)', 'Interpreter', 'Latex');
                title(a, sprintf('SD curves with Corruption Type: %s', corruptionTypes{corruptionTypeNum}));
                ylim(a, [0 1.1*maxStr]);
                xlim(a, [0 15]);
                % legend([line1(end) line2], {'Unit A', 'Unit B'});
                
                set(findall(h3(corruptionTypeNum),'-property','FontSize'),'FontSize',fontSize)
                
                if savePlot
                    fullPlotPath = [plotPath corruptionTypes{corruptionTypeNum} filesep];
                    mkdir(fullPlotPath);
                    %                 saveFigure(h1(corruptionTypeNum), [fullPlotPath sprintf('IdealGT1OverCorruption_Sim%d.png', corrInstanceNum)]);
                    %                 savefig(h1(corruptionTypeNum), [fullPlotPath sprintf('IdealGT1OverCorruption_Sim%d.fig', corrInstanceNum)]);
                    %                 saveFigure(h2(corruptionTypeNum), [fullPlotPath sprintf('IdealGT2OverCorruption_Sim%d.png', corrInstanceNum)]);
                    %                 savefig(h2(corruptionTypeNum), [fullPlotPath sprintf('IdealGT2OverCorruption_Sim%d.fig', corrInstanceNum)]);
                    saveFigure(h3(corruptionTypeNum), [fullPlotPath sprintf('SDCurveCorruption_%s_Sim%d%s.png', plotLabels{corruptionTypeNum}, corrInstanceNum, thisDblStr)]);
                    savefig(h3(corruptionTypeNum), [fullPlotPath sprintf('SDCurveCorruption_%s_Sim%d%s.fig', plotLabels{corruptionTypeNum}, corrInstanceNum, thisDblStr)]);
                end
            end
        end
    end
end

% function sdStrs = calculateSDStrs(theta, sdCurveProb, durs, str0, approxFun)



% function cost = minimizationFunction(thisGT, probTable1, probTable2)
% % This function will input a GT pair, and output a cost function analysis
% % based on the two probability tables that are also input.
% global PTableLookup diffFun
%
% P1 = PTableLookup(probTable1, thisGT(1), thisGT(2));
% P2 = PTableLookup(probTable2, thisGT(1), thisGT(2));
%
% if (P1 == 1) && (P2 == 1)
%     cost = inf; % Not ideal, but make sure that optimization doesn't get off course
% else
%     cost = diffFun(P1, P2);
% end
% end