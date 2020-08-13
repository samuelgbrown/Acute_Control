% This script will be used to explore the parameter space of various types
% of corruption and how they affect control.  Specifically, this will be
% done by rescaling the IAF equation's time axis and sensitivity (and
% possibly noise level) so that all calculations will be performed relative
% to a "standard bearer" neuron N_S

%% Set Parameters
close all;
% TODO: Make script autosave figures
doPlot = true; % Should the data be plotted
doPlotSave = true; % Should the plots be saved?
doSave = true; % Should the calculated data be saved?

figureDir = 'D:\samuelgb\Documents\MATLAB\Acute Control\Scripts\Analysis\Test\Control with Corruption\Parameter Exploration\';
matFileName = 'D:\samuelgb\Documents\MATLAB\Acute Control\Scripts\Analysis\Test\Control with Corruption\parameterExploreRescale.mat';

% Prepare the workspace to use the FP interpolation
if ~exist('approxFun', 'var')
    prepareFPWorkspace;
end

% Standard bearer true parameters
alphaS = .3;
betaS = 2*alphaS*beta; % Scale beta to 2
sigmaS = .05;
thetaS = [alphaS betaS sigmaS];

trueAlphaVal = max(allAlphas)/2; % The value of alpha that will actually be fed into the approxFun function

% Parameters for the target neurons, generated procedurally
alphaFactor = 1.5;
alphaBetaSlope = .5;

rAlpha = [alphaFactor (1/alphaFactor)];
rBeta = (rAlpha - 1)*alphaBetaSlope + 1;
rSigma = [1 1];
rTheta = [rAlpha' rBeta' rSigma'];

numThetaRatios = size(rTheta, 1);

numCorruptionPlots = 3; % Number of "ghost" curves to plot on the demonstration figure

gtSpikeRemovalInds = [8 95]; % Define the GTs that will be used for when the test neurons have spikes removed, and therefore interfere with their stimulation (short durations for fast neurons, and vice versa)
gtSpikeAdditionInds = [95 8]; % Define the GTs that will be used for when the test neurons have spikes added, and therefore interfere with the other neuron's stimulation (long durations for fast neurons, and vice versa)

% Hash settings
maxK = 1;

% Parameters for SD curve
sdCurveProb = .5;
durMin = 0;
durMax = 15;
numDurs = 100;

% Plot options
fontSize = 25;
doTitle = true;

% Create a mesh of IAF neuron parameters (alpha and beta), to be used as
% either target neurons or as corruption neurons
paramGridSize = 100;
testAlphas = linspace(.5, 2, paramGridSize + 1); % Plus one to make sure that our dimensions are always consistent (will get an error if we mix something up)
testBetas = linspace(.5, 2, paramGridSize);
[alphaMesh, betaMesh] = meshgrid(testAlphas, testBetas);
alphasCol = alphaMesh(:);
betasCol = betaMesh(:);
numIAFConds = length(alphasCol);

%% Prepare values
% Prepare the rescaled FP approximation (nondimensionalized in alpha, and
% scaled to trueAlphaVal, to stay in range of allAlphas)
gamma = (trueAlphaVal*betaS)/alphaS;
epsilon = sigmaS/sqrt(alphaS/trueAlphaVal);
approxFunRS = @(GT, thetaRatios) approxFun([GT(:, 1) ((GT(:, 2).*alphaS)./trueAlphaVal)], thetaRatios.*[trueAlphaVal gamma epsilon]);

%% Calculate the SD curve of the standard bearer neuron
durs = linspace(durMin, durMax, numDurs);

% sdCurveS = calculateSDStrs(@(GT, thetaInd) approxFun(GT, thetaS), sdCurveProb, durs, 30, 1);
sdCurveS = calculateSDStrs(@(GT, thetaInd) approxFunRS(GT, [1 1 1]), sdCurveProb, durs, 30, 1); % Should give the same result at the line above it

%% Generate the data needed for a demonstrative corruption figure, to show how corruption affects control in SD space
thisTargetThetaRatio = [2 1.5 1];
thisTargetSDCurve = calculateSDStrs(@(GT, thetaInd) approxFunRS(GT, thisTargetThetaRatio), sdCurveProb, durs, 15, 1);
thisNonTargetSDCurve = sdCurveS;

thisGTInd = gtSpikeRemovalInds(1);
otherGTInd = gtSpikeRemovalInds(2);
thisDur = durs(thisGTInd);
otherDur = durs(otherGTInd);

nonTargetStimStr = (thisNonTargetSDCurve(otherGTInd) + thisTargetSDCurve(otherGTInd))/2;
targetStimStr = (thisNonTargetSDCurve(thisGTInd) + thisTargetSDCurve(thisGTInd))/2;

allPExclusions = linspace(0, .47, numCorruptionPlots);
ghostTargetCurves = zeros(numDurs, numCorruptionPlots);
for curveInd = 1:numCorruptionPlots
    ghostTargetCurves(:, curveInd) = calculateSDStrs(@(GT, thetaInd) approxFunRS(GT, thisTargetThetaRatio)*(1 - allPExclusions(curveInd)), sdCurveProb, durs, 15, 1);
end
nonTargetCorruptedStimStr = (thisNonTargetSDCurve(otherGTInd) + ghostTargetCurves(otherGTInd, :))/2; % The stimulation strengths calculated for each corrupted target curve, assuming a simple cost function
targetCorruptedStimStr = (thisNonTargetSDCurve(thisGTInd) + ghostTargetCurves(thisGTInd, :))/2; % The stimulation strengths calculated for each corrupted target curve, assuming a simple cost function


h = figure;
a = axes;
hold(a, 'on');
colormap(jet);
colorbar;
allC = jet(2*numCorruptionPlots); % Multiplied by two, because we are only using half of the color space, because nothing above 50% corruption is plotted
plot(a, durs, thisTargetSDCurve, 'b-', durs, thisNonTargetSDCurve, 'r-', otherDur, nonTargetStimStr, 'dk', thisDur, targetStimStr, 'dk', 'linewidth', 3, 'markersize', 20);
for curveInd = 1:numCorruptionPlots
    plot(a, durs, ghostTargetCurves(:, curveInd), 'k-', 'color', allC(curveInd, :));
    plot(a, otherDur, nonTargetCorruptedStimStr(curveInd), 'd', 'color', allC(curveInd, :), 'linewidth', 2, 'markersize', 15);
    plot(a, thisDur, targetCorruptedStimStr(curveInd), 'd', 'color', allC(curveInd, :), 'linewidth', 2, 'markersize', 15);
end
caxis([0 1]);
xlabel(a, 'Duration (ms)');
ylabel(a, 'Strength (mW/mm^2)');
title(a, 'Demonstration of Exclusion Effect on Control');
legend(a, 'Target Neuron', 'Non-Target Neuron', 'Optimal Stimulations');

hStr = sprintf('%sExc_Demo.', figureDir);
saveFigure(h, [hStr 'png']);
savefig(h, [hStr 'fig']);

%% Create functions that will be used to calculate the max tolerable corruption
% Generate functions to find the maximum corruption probability that the
% target neuron could tolerate of each type, and still be controllable with
% the other neuron
fExcludeSpike = @(alpha, beta, gtInd) 1 - .5./(approxFunRS(repmat([sdCurveS(gtInd) durs(gtInd)], length(alpha), 1), [alpha, beta, ones(size(alpha))]));
fDeleteSpike = @(targetThetaRatio, gtInd) 1 - .5./(approxFunRS(repmat([sdCurveS(gtInd) durs(gtInd)], size(targetThetaRatio, 1), 1), targetThetaRatio));
fAddSpike = @(targetThetaRatio, gtInd) (.5 - (approxFunRS(repmat([sdCurveS(gtInd) durs(gtInd)], size(targetThetaRatio, 1), 1), targetThetaRatio)))./(1 - (approxFunRS(repmat([sdCurveS(gtInd) durs(gtInd)], size(targetThetaRatio, 1), 1), targetThetaRatio)));

% Create a function to simulate the hash "firing rate"
% % hashFun = @(k, eHalf, gtInd) 1./(1 + exp(-k*(sdCurveS(gtInd)*durs(gtInd)) + eHalf)); % Old (pre 7/15/20)
hashFun = @(k, eHalf, gtInd) 1./(1 + exp(-k*(sdCurveS(gtInd)*durs(gtInd) - eHalf))); % New (post 7/15/20)

%% Exclusion
% TODO: Expect the edges of the allowable alpha/beta region to have firing
% probabilities that approach .5, right?

pExcMax = zeros(size(alphaMesh, 1), size(alphaMesh, 2), numThetaRatios);
pExcT = zeros(size(alphaMesh, 1), size(alphaMesh, 2), numThetaRatios);

for thetaNum = 1:numThetaRatios
    thisGTInd = gtSpikeRemovalInds(thetaNum);
    
    if durs(thisGTInd) < (max(durs)/2)
        % If this is a small duration, then we should only be looking at
        % neurons that are faster than the standard bearer
        %         toCalc = alphasCol > 1 & betasCol > 1 & alphasCol./betasCol > 1; % For every alpha/beta combination where the ratio is greater than 1, calculate the maximum corruption probability
        toCalc = true(size(alphasCol));
        standType = 'Slow';
        durType = 'Short';
    else
        % If this is a small duration, then we should only be looking at
        % neurons that are slower than the standard bearer
        %         toCalc = alphasCol < 1 & betasCol < 1 & alphasCol./betasCol < 1; % For every alpha/beta combination where the ratio is less than 1, calculate the maximum corruption probability
        toCalc = true(size(alphasCol));
        standType = 'Fast';
        durType = 'Long';
    end
    
    pExcMaxCol = nan(numIAFConds, 1);
    pTCol = nan(numIAFConds, 1);
    for i = 1:numIAFConds
        if toCalc(i)
            pExcMaxCol(i) = fExcluSpike(alphasCol(i), betasCol(i),  thisGTInd);
            pTCol(i) = approxFunRS([sdCurveS(thisGTInd) durs(thisGTInd)], [alphasCol(i) betasCol(i) 1]);
        end
    end
    
    % Remove any P < 0 or 1 < P values
    pExcMaxCol(pExcMaxCol < 0 | 1 < pExcMaxCol) = nan;
    
    % Rearrange the values into a matrix
    thisPExcMax = zeros(size(alphaMesh));
    thisPT = zeros(size(alphaMesh));
    thisPExcMax(:) = pExcMaxCol;
    thisPT(:) = pTCol;
    
    pExcMax(:, :, thetaNum) = thisPExcMax;
    pExcT(:, :, thetaNum) = thisPT;
    
    %% Plot the results
    if doPlot
        h1 = figure;
        surf(testBetas, testAlphas, thisPExcMax');
        xlabel('$r_{\beta_T}$', 'Interpreter', 'Latex');
        ylabel('$r_{\alpha_T}$', 'Interpreter', 'Latex');
        zlabel('Max P_E_X_C');
        if doTitle
            title(sprintf('Maximum exclusion corruption tolerated for control at D=%0.2f on %s S SD curve (S=%0.2f)', durs(thisGTInd), standType, sdCurveS(thisGTInd)));
        end
        xlim([min(testBetas) max(testBetas)]);
        ylim([min(testAlphas) max(testAlphas)]);
        zlim([0 1]);
        caxis([0 1]);
        colorbar;
        view(2);
        set(findall(h1,'-property','FontSize'),'FontSize',fontSize);
        
        h2 = figure;
        surf(testBetas, testAlphas, thisPT');
        xlabel('$r_{\beta_T}$', 'Interpreter', 'Latex');
        ylabel('$r_{\alpha_T}$', 'Interpreter', 'Latex');
        zlabel('P_T');
        if doTitle
            title(sprintf('Probability of N_T firing at D=%0.2f on %s SD curve (S=%0.2f)', durs(thisGTInd), standType, sdCurveS(thisGTInd)));
        end
        xlim([min(testBetas) max(testBetas)]);
        ylim([min(testAlphas) max(testAlphas)]);
        zlim([0 1]);
        caxis([0 1]);
        colorbar;
        view(2);
        set(findall(h2,'-property','FontSize'),'FontSize',fontSize);
        
        if doPlotSave
            h1Str = sprintf('%sExc_%s_MaxExc.', figureDir, durType);
            h2Str = sprintf('%sExc_%s_Firing.', figureDir, durType);
%             h3Str = sprintf('%sExc_%s_Demo.', figureDir, durType);
            saveFigure(h1, [h1Str 'png']);
            savefig(h1, [h1Str 'fig']);
            saveFigure(h2, [h2Str 'png']);
            savefig(h2, [h2Str 'fig']);
%             saveFigure(h3, [h3Str 'png']);
%             savefig(h3, [h3Str 'fig']);
        end
    end
    
end

%% Addition
% Go through each target neuron that we are considering, and explore the
% three parameters that define the corruption
additionSDCurves = zeros(numDurs, numThetaRatios);
pIncMax = zeros(size(alphaMesh, 1), size(alphaMesh, 2), numThetaRatios);
pC = zeros(size(alphaMesh, 1), size(alphaMesh, 2), numThetaRatios);
pCPIncMax = zeros(numThetaRatios, 1);
for thetaNum = 1:numThetaRatios
    thisThetaRatio = rTheta(thetaNum, :);
    thisGTInd = gtSpikeAdditionInds(thetaNum);
    
    if thisThetaRatio(1) > 1
        % If this neuron is fast, then the standard bearer is comparatively
        % slow
        standType = 'Slow';
        durType = 'Long';
    else
        % If this neuron is slow, then the standard bearer is comparatively
        % fast
        standType = 'Fast';
        durType = 'Short';
    end
    
    % Define the SD curve of the test neuron
    additionSDCurves(:, thetaNum) = calculateSDStrs(@(GT, thetaInd) approxFunRS(GT, thisThetaRatio), sdCurveProb, durs, 15, 1);
    
    %% Explore the corrupting neuron's alpha/beta space
    thisPCPIncMax = fAddSpike(thisThetaRatio, thisGTInd); % Maximum tolerated corruption probability (P_Inc * P_C) pcpincMaxCol
    
    pIncMaxCol = nan(numIAFConds, 1); % Maximum tolerated P_Inc for each corruption neuron's alpha/beta
    pCCol = nan(numIAFConds, 1); % The firing probability of each corruptor neuron at this SD
    for i = 1:numIAFConds
        %         if toCalc(i)
        pCCol(i) = approxFunRS([sdCurveS(thisGTInd) durs(thisGTInd)], [alphasCol(i) betasCol(i) 1]);
        pIncMaxCol(i) = thisPCPIncMax/pCCol(i); % Find the maximum P_Inc, which must be equal to thisPCPIncMax when multiplied by P_C
        %         end
    end
    
    % Remove any P < 0 or 1 < P values
    pIncMaxCol(pIncMaxCol < 0 | 1 < pIncMaxCol) = nan;
    
    % Rearrange the values into a matrix
    thisPIncMax = zeros(size(alphaMesh));
    thisPC = zeros(size(alphaMesh));
    thisPIncMax(:) = pIncMaxCol;
    thisPC(:) = pCCol;
    
    pIncMax(:, :, thetaNum) = thisPIncMax;
    pC(:, :, thetaNum) = thisPC;
    pCPIncMax(thetaNum) = thisPCPIncMax;
    
    %% Plot the results
    if doPlot
        h1 = figure;
        surf(testBetas, testAlphas, thisPIncMax');
        xlabel('$r_{\beta_C}$', 'Interpreter', 'Latex');
        ylabel('$r_{\alpha_C}$', 'Interpreter', 'Latex');
        zlabel('Max P_I_N_C');
        if doTitle
            title(sprintf('Maximum inclusion corruption tolerated for control at D=%0.2f on %s S SD curve (S=%0.2f), with non-target neuron [%0.1f %0.1f %0.1f]', durs(thisGTInd), standType, sdCurveS(thisGTInd), thisThetaRatio(1), thisThetaRatio(2), thisThetaRatio(3)));
        end
        xlim([min(testBetas) max(testBetas)]);
        ylim([min(testAlphas) max(testAlphas)]);
        zlim([0 1]);
        caxis([0 1]);
        colorbar;
        view(2);
        set(findall(h1,'-property','FontSize'),'FontSize',fontSize);
        
        h2 = figure;
        surf(testBetas, testAlphas, thisPC');
        xlabel('$r_{\beta_C}$', 'Interpreter', 'Latex');
        ylabel('$r_{\alpha_C}$', 'Interpreter', 'Latex');
        zlabel('P_C');
        if doTitle
            title(sprintf('Probability of N_C firing at D=%0.2f on S SD curve (S=%0.2f)', durs(thisGTInd), sdCurveS(thisGTInd)));
        end
        xlim([min(testBetas) max(testBetas)]);
        ylim([min(testAlphas) max(testAlphas)]);
        zlim([0 1]);
        caxis([0 1]);
        colorbar;
        view(2);
        set(findall(h2,'-property','FontSize'),'FontSize',fontSize);
        
        if doPlotSave
            h1Str = sprintf('%sAdd_%s_MaxInc.', figureDir, durType);
            h2Str = sprintf('%sAdd_%s_Firing.', figureDir, durType);
            saveFigure(h1, [h1Str 'png']);
            savefig(h1, [h1Str 'fig']);
            saveFigure(h2, [h2Str 'png']);
            savefig(h2, [h2Str 'fig']);
        end
    end
end

%% Hash Deletion
% Create a mesh of logistic parameters (slope and offset), to be used as
% descriptors of the deletion hash
paramGridSize = 1000;
eHalfMaxRatio = 2; % The ratio between the highest eHalf to be tested and the maximum energy stimulation to be considered
[strMax, strMaxInd] = max(sdCurveS);
eMax = durs(strMaxInd)*strMax; % The maximum energy stimulation to be considered
% TODO: Figure out some principled way of selecting limits for k and EHalf
testKs = linspace(0, maxK, paramGridSize);
testEHalfs = linspace(0, eHalfMaxRatio*eMax, paramGridSize);
[kMesh, eHalfMesh] = meshgrid(testKs, testEHalfs);
ksCol = kMesh(:);
eHalfCol = eHalfMesh(:);
numHashConds = length(ksCol);

pDelMax = zeros(size(kMesh, 1), size(kMesh, 2), numThetaRatios);
pHash = zeros(size(kMesh, 1), size(kMesh, 2), numThetaRatios);
pHashPDelMax = zeros(numThetaRatios, 1);

deletionSDCurves = zeros(numDurs, numThetaRatios);
for thetaNum = 1:numThetaRatios
    thisThetaRatio = rTheta(thetaNum, :);
    thisGTInd = gtSpikeRemovalInds(thetaNum);
    
    if thisThetaRatio(1) > 1
        % If this neuron is fast, then the standard bearer is comparatively
        % slow
        standType = 'Slow';
        durType = 'Short';
    else
        % If this neuron is slow, then the standard bearer is comparatively
        % fast
        standType = 'Fast';
        durType = 'Long';
    end
    
    % Define the SD curve of the test neuron
    deletionSDCurves(:, thetaNum) = calculateSDStrs(@(GT, thetaInd) approxFunRS(GT, thisThetaRatio), sdCurveProb, durs, 15, 1);
    
    %% Explore the corrupting neuron's alpha/beta space
    thisPHashPDelMax = fDeleteSpike(thisThetaRatio, thisGTInd); % Maximum tolerated corruption probability (P_Del * P_Hash)
    
    pDelMaxCol = nan(numHashConds, 1); % Maximum tolerated P_Del for each hash instance
    pHashCol = nan(numHashConds, 1); % Firing probability of each hash corruption at this SD
    for i = 1:numHashConds
       pHashCol(i) = hashFun(ksCol(i), eHalfCol(i), thisGTInd);
       pDelMaxCol(i) = thisPHashPDelMax/pHashCol(i); % Find the maximum probability of P_Del, which must multiply with the "firing rate" of the hash to produce thisPDelPHashMax.
    end
    
    % Remove any P < 0 or 1 < P values
    pDelMaxCol(pDelMaxCol < 0) = nan;
    pDelMaxCol(1 < pDelMaxCol) = 1;
    
    % Rearrange the values into a matrix
    thisPDelMax = zeros(size(kMesh));
    thisPHash = zeros(size(kMesh));
    thisPDelMax(:) = pDelMaxCol;
    thisPHash(:) = pHashCol;
    
    pDelMax(:, :, thetaNum) = thisPDelMax;
    pHash(:, :, thetaNum) = thisPHash;
    pHashPDelMax(thetaNum) = thisPHashPDelMax;
    
    %% Plot the results
    if doPlot
        h1 = figure;
        surf(testEHalfs, testKs, thisPDelMax', 'edgealpha', 0);
        xlabel('E_H_a_l_f');
        ylabel('k');
        zlabel('Max P_D_E_L');
        if doTitle
            title(sprintf('Maximum deletion corruption tolerated for control at D=%0.2f on %s S SD curve (S=%0.2f), with target neuron [%0.1f %0.1f %0.1f]', durs(thisGTInd), standType, sdCurveS(thisGTInd), thisThetaRatio(1), thisThetaRatio(2), thisThetaRatio(3)));
        end
        xlim([min(testEHalfs) max(testEHalfs)]);
        ylim([min(testKs) max(testKs)]);
        zlim([0 1]);
        caxis([0 1]);
        colorbar;
        view(2);
        set(findall(h1,'-property','FontSize'),'FontSize',fontSize);
        
        h2 = figure;
        surf(testEHalfs, testKs, thisPHash', 'edgealpha', 0);
        xlabel('E_H_a_l_f');
        ylabel('k');
        zlabel('P_H_a_s_h');
        if doTitle
            title(sprintf('Probability of hash ''firing'' at D=%0.2f on %s S SD curve (S=%0.2f)', durs(thisGTInd), standType, sdCurveS(thisGTInd)));
        end
        xlim([min(testEHalfs) max(testEHalfs)]);
        ylim([min(testKs) max(testKs)]);
        zlim([0 1]);
        caxis([0 1]);
        colorbar;
        view(2);
        set(findall(h2,'-property','FontSize'),'FontSize',fontSize);
        
        if doPlotSave
            h1Str = sprintf('%sHash_%s_MaxDel.', figureDir, durType);
            h2Str = sprintf('%sHash_%s_Firing.', figureDir, durType);
            saveFigure(h1, [h1Str 'png']);
            savefig(h1, [h1Str 'fig']);
            saveFigure(h2, [h2Str 'png']);
            savefig(h2, [h2Str 'fig']);
        end
    end
end

%%
if doSave
   save(matFileName);
end