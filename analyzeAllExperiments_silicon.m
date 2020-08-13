% 2/16/18
% This script will analyze all data from across multiple acute control
% experiments
%
% BEFORE RUNNING THIS SCRIPT:
% Run "analyzeExperiment.m" on all data that is to be included in the
% resulting data analysis.  Set "saveAllData" to "true" before running the
% script to save the data to "analysis" .mat files.  This script assumes
% that all experimental days being analyzed have the same general format,
% including number of stims per sequence, sequences per run, runs per set,
% and sets per pair.  It can handle the experimental days having different
% adaptation information (i.e. can handle some having adaptive while others
% do not).

close all;

%% Set parameters
useOfflineSort = false; % Should the offline sorted data be displayed (alternative: using online sorted data)

doStatisticalTests = true; % Should the statistical tests be run to determine if control was effective, if stimulations produced significantly different results, and if a significant number of blocks were controllable?

reloadData = true; % Should the data from all of the experiments be reloaded?  Set this to true when analyzeExperiment must be re-run for some reason
% onlyISI100 = false; % Only analyze sets with an ISI of 100
useDiff = true; % Use the diff cost function instead of the sqe cost function (for performance plot)?
useWhiteList = false; % Should the analysis be restricted to only trials that are explicity while-listed in the dataToAnalyze.txt files within each data directory?
useClustParamLimit = false; % Should the data included in the analysis be limited by the clustering parameter?
doConsistentUnitLabeling = true; % Should the program ensure that unit A is always labeled such that it has a shorter stimulation, and vice versa for B?
doAcceptableControlabilityLabeling = true; % Should the figure show which units demonstrate "acceptable" control?

clustParam = 'Lratio'; % The parameter from the clustering statistics that should be used
lRatioLimit = .263; % The highest L_ratio that will be accepted in the performance plot

numExamples = 3; % Number of data sets to show examples of
exampleQuartilesFromControllable = true; % If doAcceptableControlabilityLabeling, should the examples be taken from the quartiles of the entire data set, or only of the controllable pairs (even the worst example will be controllable)
numTracesMax = 500; % (Maximum) number of traces to show in the plot quartiles charts
numSamples = 32; % Samples per waveform
numSamplesPrePeak = 5; % Samples before peak
numShuffles = 20000; % Number of dataset shuffles to do to generate a distribution to test hypotheses with

% Plotting
doStatisticalTestingPlot = true; % Should the statistical test results be plotted?
doPerformancePlot = false; % Should the overall statistics of the experiment be calculated?
doQuartilesPlot = false; % (If doPerformancePlot == true) Should examples in the CDF be plotted?
doQualityVsPerformancePlot = false; % (If doPerformancePlot == true) Should the quality vs performance plot be made?
doStringPlot = false; % Plot the string plot, a comparison between the costs from online to offline sorting
doProbSpacePlot = false; % Plot the probability space plot, showing where each stim-type (SA and SB) appears in probability space.
doCostFunctionPlot = false; % Should the results from different cost functions be compared?
doSummarySpreadsheet = false; % Should the clustering statistics and performances of each experiment be saved to a spreadsheet?

dataSaveFileBase = 'D:\samuelgb\Documents\MATLAB\Acute Control\Data\Silicon Probe Data\allExperimentsAnalysis';
doSaveData = true;
doSaveFigs = true; % Should the figures be saved?

plotFontSize = 24; % The font size in the plot
plotLineWidth = 2;
plotMarkerSize = 15;
smoothingSize = 100; % Number of stimuli to smooth over while plotting the POT

% Parameters from the experiment
sizeAdaptationBlock = 10; % Number of sequences per adaptation block
costFunLambda = 1e-5; %.01;
adaptiveCostFunction = @(P1, P2, G) -(P1 .* (1 - P2)) + costFunLambda*(G.^2); % Cost function to optimize using control
hdrMassFrac = .95; % The fraction of mass in a probability distribution that constitutes the high density region

% Quartiles plot options
doWaveformPlots = false; % Should the waveforms be plotted along with the POT in the quartile examples

% Probability plot parameters
doConnectPoints = true; % Should the SA and SB from a single dataset be connected?
doVectorField = true; % Should the vectorfield of the cost sensitivity be plotted?
vectorScale = .4;
vectorFieldMeshSize = 5;

% Analysis parameters
squaredMeanPower = 10; % Used in the computation of the squaredMeanCostCombine function (higher values seem to lead to better ordering of the data)
partialTimeLim = 60; % The time limit (in seconds) to include data for a "partial" sucess calculation
% performanceTiers = [-Inf 0 .1 .25 .4]; % The tiers of performance that will be checked (in terms of TP_desired - FA_undesired)
% performanceTiers = [-Inf 0 .05 .1 .15]; % The tiers of performance that will be checked (in terms of TP_desired - FA_undesired)
% numPerformanceTiers = length(performanceTiers); % Number of actual tiers to calculate
numPerformanceTiers = 5; % Number of actual tiers to calculate
numCostFunctionPerformanceTiers = 5; % The number of tiers to group data into according to cost (actual tier values are calculated on the fly)
numCostTimePeriods = 2; % Number of time periods to compare costs over

% Data to exclude
% Excluded for recording issues
% excludeDateStrs = {'20170620', '20170623', '20170711', '20170718', '20170817', '20170906', '20170912', '20170912', '20180102'};
% excludeTimeStrs = {'1249', '1620', '1403', '1312', '1236', '1220', '1237', '1359', '1402'};
excludeDateStrs = {'20190830', '20190830', '20191030', '20191030', '20191030', '20191030', '20191030'};
excludeTimeStrs = {'1715', '1829', '1249', '1341', '1421', '1520', '1605'};
% Excluded for spike sorting issues (TODO: REMOVE WHEN DATA HAS BEEN
% RESORTED
% excludeDateStrs = [excludeDateStrs, '20170406', '20170630', '20170711', '20170720', '20170727', '20171103'];
% excludeTimeStrs = [excludeTimeStrs, '1456', '1220', '1403', '1659', '1308', '1557'];

% Directories
mainDataDir = 'D:\samuelgb\Documents\MATLAB\Acute Control\Data\Silicon Probe Data';
parentFigureDir = 'D:\samuelgb\Documents\MATLAB\Acute Control\Figures\Silicon Probe Data\Total\';
metaDataFileName = 'dataToAnalyze.txt'; % The name of the text file that shows which .mat files to load
analysisDirName = 'AnalysisData';
analysisFileBase = 'analysis';
cumulativeDataDir = [mainDataDir filesep 'Cumulative Data'];
xlsFileName = [mainDataDir filesep 'clusterStatistics.xls'];
fileNameBase = 'cumulativeData';
if useOfflineSort
    sortSuffix = 'off';
else
    sortSuffix = 'on';
end

dataSaveFile = [dataSaveFileBase '_' sortSuffix '.mat'];

fileNameBase = [fileNameBase '_' sortSuffix];

if doWaveformPlots
    cumulativeDataFileName = [cumulativeDataDir filesep fileNameBase '.mat'];
else
    cumulativeDataFileName = [cumulativeDataDir filesep fileNameBase '.mat'];
end

% Conversion factor to turn the trace units to microvolts (they are in
% 1e-4V, or 1e2µV, now)
convToMicroV = 100;
difCostFn = @(difProb1, difProb2) -(difProb1 .* (1 - difProb2));
sqeCostFn  = @(difProb1, difProb2) (1/2)*(((1 - difProb1).^2) + difProb2.^2);

% Define the functions used to combine the costs from each electrode
meanCostCombine = @(costA, costB) (costA + costB)/2;
squaredMeanCostCombine = @(costA, costB) (costA.^squaredMeanPower + costB.^squaredMeanPower)/2; % Increasing the power in this expression seems to lead to better ordering of the data
diffCostCombine = @(costA, costB) abs(costA - costB);
minDiffAdjustedCostCombine = @(costA, costB) (2/3)*(diffCostCombine(costA, costB) + meanCostCombine(costA, costB)); % Raw expression has a range of [0 1.5], so scale it to fit on [0 1]

dirsToRead = 0; % Which directories to read(SET TO "0" TO ANALYZE ALL)

%% Check input
if length(excludeDateStrs) ~= length(excludeTimeStrs)
    error('Excluded data strings and excluded data times do not match (length(excludeDataStrs) ~= length(excludeTimeStrs))');
end

%% Load the data from all of the experiments, if needed
if ~exist(cumulativeDataFileName, 'file') || reloadData
    originalDir = pwd;
    isDataDir = @(a) ~isempty(regexp(a, '^\d{8}$', 'start')); % Is "a" formatted as a data directory (e.g. "20170723")
    isAnalysisFile = @(a) length(a) > length(analysisFileBase) && strcmp(a(1:length(analysisFileBase)),analysisFileBase);
    shouldAnalyzeThisDir = @(a) any(cell2mat(cellfun(@(c) strcmp(a, c), dirsToRead, 'uniformoutput', false))); % Is directory name "a" part of the string cell array "dirsToRead"
    readyAnyData = false;
    
    % Create variables to store the metadata
    includedDateStrs = cell(1); % All of the date strings
    includedTimeStrs = cell(1); % All of the time strings
    includedTrialNums = zeros(1); % All of the trial numbers
    includedInd = 1; % The index in the above "included" variables
    excluded = struct('dateStrs', cell(1), 'timeStrs', cell(1), 'trialNum', zeros(1), 'num', 0); % Keep track of all excluded pairs
    
    % Note that the data variables have not been created yet (only do so
    % when we know how large everything needs to be)
    preparedDataVars = false;
    
    readAllDirs = ~iscell(dirsToRead) && all(dirsToRead == 0); % Should all directories be read? (Is dirsToRead equal to 0?)
    cd(mainDataDir);
    mainDirContents = dir;
    numContents = length(mainDirContents);
    
    % Create online and offline structs to store both data sets
    sTemp = struct;
    
    for dirInd = 1:numContents
        thisDir = mainDirContents(dirInd);
        if isDataDir(thisDir.name)
            % This is a data directory
            cd(thisDir.name);
            pathParts = strsplit(pwd, filesep);
            dateStr = pathParts{end};
            if readAllDirs || any(strcmp(dateStr, dirsToRead))
                % If this directory should be analyzed (dirsToRead is 0, this
                % directory is on the whitelist)
                thisDataDir = [thisDir.folder filesep thisDir.name];
                thisAnalysisDir = [thisDataDir filesep analysisDirName];
                
                if exist(thisAnalysisDir, 'file')
                    % The analysis directory exists
                    
                    % First, check the dataToAnalyze file to find which
                    % units should be analyzed (if needed)
                    dataToAnalyzeFileName = [thisDataDir filesep metaDataFileName]; % The name of the data to analyze file (if it exists)
                    %             if ((length(fileName) >= length(metaDataFileName)) && strcmp(fileName(1:length(metaDataFileName)), metaDataFileName))
                    if useWhiteList && exist(dataToAnalyzeFileName, 'file')
                        % Load this meta data
                        fid = fopen(dataToAnalyzeFileName);
                        thisFilesToLoad = cell2mat(cellfun(@str2double, strsplit(fgetl(fid), ','), 'uniformoutput', false));
                        fclose(fid);
                        
                        if (length(thisFilesToLoad) == 1 && (thisFilesToLoad == -1))% && (length(thisMCFilesToLoad) == 1 && (thisMCFilesToLoad == -1))
                            % No files in this directory should be read, so skip
                            % it.
                            break;
                        end
                    else
                        thisFilesToLoad = 0; % If there is a metaDataFile in this directory, then use the input filesToLoad
                    end
                    
                    % Next, arrange the data in this directory to
                    % understand which experiments are being referenced by
                    % each file number
                    allNeuronStructFiles = dir('neuronStruct*.mat');
                    dataToAnalyzeTimeStrs = cell(1);
                    fileNum = 1;
                    % Get all the time strings, in order, for data on this
                    % date
                    for fileInd = 1:length(allNeuronStructFiles)
                        % Extract the time string
                        thisAnnoyingCell = regexp(allNeuronStructFiles(fileInd).name, '_(\d{4}).mat', 'tokens');
                        if isempty(thisAnnoyingCell)
                            continue; % This data set doesn't have a valid time string...for some reason...
                        end
                        thisTimeStr = thisAnnoyingCell{1};
                        dataToAnalyzeTimeStrs(fileNum) = thisTimeStr;
                        fileNum = fileNum + 1;
                    end
                    
                    % Finally, remove all time strings that are not
                    % white-listed
                    if thisFilesToLoad ~= 0
                        dataToAnalyzeTimeStrs = dataToAnalyzeTimeStrs(thisFilesToLoad);
                    end
                    
                    cd(thisAnalysisDir);
                    analysisDirContents = dir;
                    numPotentialAnalysisFiles = length(analysisDirContents);
                    fprintf('Starting %s...\n', dateStr);
                    
                    for analysisDirInd = 1:numPotentialAnalysisFiles
                        thisAnalysisFileName = analysisDirContents(analysisDirInd).name;
                        if isAnalysisFile(thisAnalysisFileName)
                            % Extract information from the file name
                            nameParts = strsplit(thisAnalysisFileName, '_');
                            thisDateStr = nameParts{3}; % should be the same as dateStr defined above
                            thisTimeStr = nameParts{4};
                            thisTrial = str2double(nameParts{5}(2)); % Get the second character (the first is 't')
                            sortType = nameParts{2};
                            
                            if strcmp(sortType, 'On') && ~useOfflineSort
                                % If this is an online dataset, check that
                                % the offline dataset exists, and then
                                % continue
                            elseif strcmp(sortType, 'Off') && useOfflineSort
                                
                            else
                                % If this is a sort type (online or
                                % offline) of the type that we are NOT
                                % using, then skip this file
                                continue;
                            end
                            
                            % Keep track of the file number (used to
                            % determine which experiments should be
                            % white-listed according to the dataToAnalyze
                            % file
                            if thisTrial == 1
                                % If this is the first trial of this
                                % experment, iterate the file number
                                fileNum = fileNum + 1;
                            end
                            
                            % Check if this data set is excluded or not
                            if ~any(strcmp(dataToAnalyzeTimeStrs, thisTimeStr)) || any(strcmp(thisDateStr, excludeDateStrs) & strcmp(thisTimeStr, excludeTimeStrs))
                                excluded.num = excluded.num + 1;
                                excluded.dateStrs{excluded.num} = thisDateStr;
                                excluded.timeStrs{excluded.num} = thisTimeStr;
                                excluded.trialNum(excluded.num) = thisTrial;
                                fprintf('Excluded %s_%s...\n', thisDateStr, thisTimeStr);
                                continue;
                            end
                            
                            % Load this data into the work space
                            load(thisAnalysisFileName);
                            
                            % Label the neurons A and B according to their
                            % stim types
                            if doConsistentUnitLabeling && round(mean(allT1 > allT2))
                                % Want neuron 1/stim A to be the SHORTER
                                % stimulation.  If this is NOT the case,
                                % and we want consistent labeling, then
                                % switch the data
                                [adaptiveInfo1, adaptiveInfo2] = deal(adaptiveInfo2, adaptiveInfo1);
                                [allT1, allT2] = deal(allT2, allT1);
                                if doWaveformPlots
                                    [vT1, vT2] = deal(vT2, vT1);
                                end
                                eRAll1 = ~eRAll1;
                                spikeHitInds = flipud(spikeHitInds);
                                [thisTrialResultsData1, thisTrialResultsData2] = deal(thisTrialResultsData2, thisTrialResultsData1);
                                
                                % Load the sorting statistics
                                clustStatAFileName = sprintf('%s%s%s_%s_B_01-CluQual.mat', thisDataDir, filesep, thisDateStr, thisTimeStr);
                                clustStatBFileName = sprintf('%s%s%s_%s_A_01-CluQual.mat', thisDataDir, filesep, thisDateStr, thisTimeStr);
                                fprintf('\nSwapped N_A N_B in data set %s_%s.\n', thisDateStr, thisTimeStr);
                            else
                                % Load the sorting statistics
                                clustStatAFileName = sprintf('%s%s%s_%s_A_01-CluQual.mat', thisDataDir, filesep, thisDateStr, thisTimeStr);
                                clustStatBFileName = sprintf('%s%s%s_%s_B_01-CluQual.mat', thisDataDir, filesep, thisDateStr, thisTimeStr);
                            end
                            
                            if useOfflineSort
                                if exist(clustStatAFileName, 'file') && exist(clustStatBFileName, 'file')
                                    load(clustStatAFileName);
                                    if ~isfield(CluSep, 'spikeInds')
                                        CluSep.spikeInds = nan;
                                    end
                                    CluSepA = CluSep;
                                    
                                    load(clustStatBFileName);
                                    if ~isfield(CluSep, 'spikeInds')
                                        CluSep.spikeInds = nan;
                                    end
                                    CluSepB = CluSep;
                                else
                                    fprintf('Missing cluster file from %s_%s\n', thisDateStr, thisTimeStr);
                                    CluSepA = [];
                                    CluSepB = [];
                                end
                            end
                            
                            % Store the loaded data
                            if ~preparedDataVars
                                % Create variables to store the data
                                sTemp.cAdaptiveInfo1{1} = adaptiveInfo2;
                                sTemp.cAdaptiveInfo2{1} = adaptiveInfo1;
                                sTemp.cAllT1 = allT1';
                                sTemp.cAllT2 = allT2';
                                sTemp.cDoAdapt = doAdapt;
                                sTemp.cSpikeHitInds{1} = spikeHitInds;
                                sTemp.cEpochStartInds{1} = epochStartInds;
                                sTemp.cERAll = eRAll1;
                                sTemp.cNumSequences = numSequences;
                                sTemp.cNumUnits = numUnits;
                                sTemp.cRunOrders = runOrdersT;
                                sTemp.cStimStrengthsThisTrial = stimStrengthsThisTrial;
                                sTemp.cStimDurationsThisTrial = stimDurationsThisTrial;
                                sTemp.cStimTimes = stimTimes1;
                                sTemp.cThisIS = iS;
                                sTemp.cThisTrialResultsData1 = thisTrialResultsData1;
                                sTemp.cThisTrialResultsData2 = thisTrialResultsData2;
                                
                                if useOfflineSort
                                    sTemp.cClustStat{1} = [CluSepA, CluSepB];
                                end
                                
                                if doWaveformPlots
                                    sTemp.cVT1{1} = vT1;
                                    sTemp.cVT2{1} = vT2;
                                end
                                
                                preparedDataVars = true;
                            else
                                % Store the data in the previously created
                                % variables
                                sTemp.cAdaptiveInfo1{includedInd} = adaptiveInfo1;
                                sTemp.cAdaptiveInfo2{includedInd} = adaptiveInfo2;
                                sTemp.cAllT1(includedInd, :) = allT1';
                                sTemp.cAllT2(includedInd, :) = allT2';
                                sTemp.cDoAdapt(includedInd) = doAdapt;
                                sTemp.cSpikeHitInds{includedInd} = spikeHitInds;
                                sTemp.cEpochStartInds{includedInd} = epochStartInds;
                                sTemp.cERAll(:,:,includedInd) = eRAll1;
                                sTemp.cNumSequences(includedInd) = numSequences;
                                sTemp.cNumUnits(includedInd) = numUnits;
                                sTemp.cRunOrders(:,:,includedInd) = runOrdersT;
                                sTemp.cStimStrengthsThisTrial(:,:,includedInd) = stimStrengthsThisTrial;
                                sTemp.cStimDurationsThisTrial(:,:,includedInd) = stimDurationsThisTrial;
                                sTemp.cStimTimes(:,:,includedInd) = stimTimes1;
                                sTemp.cThisIS(includedInd) = iS;
                                sTemp.cThisTrialResultsData1(:,:,includedInd) = thisTrialResultsData1;
                                sTemp.cThisTrialResultsData2(:,:,includedInd) = thisTrialResultsData2;
                                
                                if useOfflineSort
                                    sTemp.cClustStat{includedInd} = [CluSepA, CluSepB];
                                end
                                
                                if doWaveformPlots
                                    sTemp.cVT1{includedInd} = vT1;
                                    sTemp.cVT2{includedInd} = vT2;
                                end
                            end
                            
                            % If the data was successfully loaded, store
                            % the metadata
                            includedDateStrs{includedInd} = thisDateStr;
                            includedTimeStrs{includedInd} = thisTimeStr;
                            includedTrialNums(includedInd) = thisTrial;
                            
                            includedInd = includedInd + 1;
                            
                            % Mark that data has been read
                            readyAnyData = true;
                            
                            fprintf('Loaded %s_%s...\n', thisDateStr, thisTimeStr);
                        end
                    end
                    fprintf('\n');
                end
            end
        end
        cd(mainDataDir);
    end
    
    if readyAnyData
        % If data was loaded in from the analysis of the individual experiments, save the file so that this process does not need to occur every time
        % Assign all data to a struct
        
        % Record data
        s.adaptiveInds1 = sTemp.cAdaptiveInfo1;
        s.adaptiveInds2 = sTemp.cAdaptiveInfo2;
        s.allT1 = sTemp.cAllT1;
        s.allT2 = sTemp.cAllT2;
        s.doAdapt = sTemp.cDoAdapt;
        s.spikeHitInds = sTemp.cSpikeHitInds;
        s.epochStartInds = sTemp.cEpochStartInds;
        s.eRAll = sTemp.cERAll;
        s.numSequences = sTemp.cNumSequences;
        s.numUnits = sTemp.cNumUnits;
        s.runOrders = sTemp.cRunOrders;
        s.stimStrengths = sTemp.cStimStrengthsThisTrial;
        s.stimDurations = sTemp.cStimDurationsThisTrial;
        s.stimTimes = sTemp.cStimTimes;
        s.thisIS = sTemp.cThisIS;
        s.trialResultsData1 = sTemp.cThisTrialResultsData1;
        s.trialResultsData2 = sTemp.cThisTrialResultsData2;
        
        if useOfflineSort
            s.clustStat = sTemp.cClustStat;
        end
        
        if doWaveformPlots
            s.vT1 = sTemp.cVT1;
            s.vT2 = sTemp.cVT2;
        end
        
        s.numData = includedInd - 1;
        s.stimsPerSeq = size(eRAll1,2);
        s.seqsPerSet = size(eRAll1,1);
        s.dateStrs = includedDateStrs;
        s.timeStrs = includedTimeStrs;
        s.trialNums = includedTrialNums;
        s.excluded = excluded;
        
        if ~exist(cumulativeDataDir, 'dir')
            mkdir(cumulativeDataDir); % Make sure that the directory exists before saving the cumulative data file
        end
        save(cumulativeDataFileName, '-struct', 's', '-v7.3');
        
        fprintf('\n\n*\n**\n***\nSaved file %s\n***\n**\n*\n\n\n', cumulativeDataFileName);
        
    else
        fprintf('Something went wrong: No data ready.\n\n');
        return;
    end
    
    cd(originalDir); % Take the user back to the directory they were originally in
end

load(cumulativeDataFileName); % Load the cumulative data

%% Analyze data
numStimsPerSequence = size(eRAll,2);
numAllStims = numStimsPerSequence*numSequences(1);
numPreAdapt = sum(~doAdapt);
numPostAdapt = sum(doAdapt);
uniqueIS = 100;
isiStr = '';
includedISI = uniqueIS;
includedISIInds = ismember(thisIS, includedISI);
numPostAdaptISI = sum(includedISIInds & doAdapt);

numISIs = length(uniqueIS); % Number of unique inter-spike intervals
thisISPost = thisIS(doAdapt);

% allDataInds = 1:numData;
allDataInds = find(doAdapt);

[~, uniqueISLocs] = ismember(thisIS, uniqueIS);
[~, uniqueISLocsPost] = ismember(thisISPost, uniqueIS);

uniqueISInds = false(numISIs, numData);
uniqueISIndsPost = false(numISIs, numPostAdapt);
for i = 1:numISIs
    uniqueISInds(i,:) = uniqueISLocs == i;
    uniqueISIndsPost(i,:) = uniqueISLocsPost == i;
end

% Group all data into performance categories

% Organize the stimulation results into a stimsPerSet x stimsPerSequence x
% numData matrix
stimResults1 = permute(reshape(squeeze(trialResultsData1(:,3,:)), [stimsPerSeq, seqsPerSet, numData]), [2 1 3]);
stimResults2 = permute(reshape(squeeze(trialResultsData2(:,3,:)), [stimsPerSeq, seqsPerSet, numData]), [2 1 3]);
eR1 = eRAll;
eR2 = ~eRAll;

% Get data from offline sorting
stimResultsOnline1 = permute(reshape(squeeze(trialResultsData1(:,3,:)), [stimsPerSeq, seqsPerSet, numData]), [2 1 3]);
stimResultsOnline2 = permute(reshape(squeeze(trialResultsData2(:,3,:)), [stimsPerSeq, seqsPerSet, numData]), [2 1 3]);

% Manipulate eRAll in wierd ways because Matlab isn't really all that great
% about it
tmpERAll = permute(eRAll, [2 1 3]); % numStimsPSeq x numSeqs x numPostAdapt
spikeER1 = nan(numAllStims, numData);
spikeER1(:) = tmpERAll;
spikeER1 = logical(spikeER1)';
spikeER2 = ~logical(spikeER1);

tmpSR1 = permute(stimResults1, [2 1 3]);
tmpSR2 = permute(stimResults2, [2 1 3]);
spikeSR1 = nan(numAllStims, numData); % Temporary 
spikeSR2 = nan(numAllStims, numData); % Temporary 
spikeSR1(:) = tmpSR1;
spikeSR2(:) = tmpSR2;
spikeSR1 = logical(spikeSR1'); % numData x numAllStims
spikeSR2 = logical(spikeSR2'); % numData x numAllStims

tmpStimTimesAll = permute(stimTimes, [2 1 3]); % numStimsPSeq x numSeqs x numPostAdapt
stimTimesAll = nan(numAllStims, numData);
stimTimesAll(:) = tmpStimTimesAll;
stimTimesAll = stimTimesAll';

spikeSS = nan(numAllStims, numData);
spikeSS(:) = stimStrengths;
spikeSS = spikeSS';

spikeSD = nan(numAllStims, numData);
spikeSD(:) = stimDurations;
spikeSD = spikeSD';

if doPerformancePlot || doSummarySpreadsheet
    %% Prepare the time-series performance data
    % Get the logical matrices for true positive and false alarm for both units
    
    n = 0:(smoothingSize - 1);
    w = cos((pi*n)/(2*(smoothingSize - 1)) - pi/2); % Causal
    w = w/sum(w); % Normalize the window so it does not scale the data
    nTS = numAllStims - smoothingSize + 1; % Get the number of time steps
    
    % Get the response fraction values (both TP and FA) for each sequence
    tPFrac1 = nan(nTS, numData);
    fAFrac1 = nan(nTS, numData);
    tPFrac2 = nan(nTS, numData);
    fAFrac2 = nan(nTS, numData);
    
    TP1Avg = nan(numData,1);
    FA1Avg = nan(numData,1);
    TP2Avg = nan(numData,1);
    FA2Avg = nan(numData,1);
    
    for dataNum = 1:numData
        for timeStep = 1:nTS
            thisInds1 = find(spikeER1(dataNum, 1:(smoothingSize + timeStep - 1)), smoothingSize, 'last');
            thisInds2 = find(spikeER2(dataNum, 1:(smoothingSize + timeStep - 1)), smoothingSize, 'last');
            
            if length(thisInds1) >= smoothingSize
                thisSR1 = spikeSR1(dataNum, thisInds1);
                thisSR2 = spikeSR2(dataNum, thisInds1);
                
                tPFrac1(timeStep, dataNum) = sum(w.*thisSR1, 2);
                fAFrac2(timeStep, dataNum) = sum(w.*thisSR2, 2);
            end
            
            if length(thisInds2) >= smoothingSize
                thisSR1 = spikeSR1(dataNum, thisInds2);
                thisSR2 = spikeSR2(dataNum, thisInds2);
                
                fAFrac1(timeStep, dataNum) = sum(w.*thisSR1, 2);
                tPFrac2(timeStep, dataNum) = sum(w.*thisSR2, 2);
            end
        end
        
        % Calculate the average performance for each case
        TP1Avg(dataNum) = sum(spikeSR1(dataNum, logical(spikeER1(dataNum, :))))/sum(spikeER1(dataNum, :));
        FA1Avg(dataNum) = sum(spikeSR1(dataNum, logical(spikeER2(dataNum, :))))/sum(spikeER2(dataNum, :));
        TP2Avg(dataNum) = sum(spikeSR2(dataNum, logical(spikeER2(dataNum, :))))/sum(spikeER2(dataNum, :));
        FA2Avg(dataNum) = sum(spikeSR2(dataNum, logical(spikeER1(dataNum, :))))/sum(spikeER1(dataNum, :));
    end
    
    stimTimeInds = smoothingSize:(smoothingSize + nTS - 1); % Causal
    stimTimesAll = stimTimesAll(:, stimTimeInds);
    stimTimesAll = stimTimesAll - stimTimesAll(:,1); % Make all stim times start with 0
    
    %% Prepare the per-block performance data
    allBlockBorders = 1:sizeAdaptationBlock*numStimsPerSequence:(numAllStims + 1);
    numAdaptationBlocks = length(allBlockBorders) - 1;
    TP1_block = ones(numAdaptationBlocks, numData);
    TP2_block = ones(numAdaptationBlocks, numData);
    FA1_block = ones(numAdaptationBlocks, numData);
    FA2_block = ones(numAdaptationBlocks, numData);
    costs1_block = ones(numAdaptationBlocks, numData);
    costs2_block = ones(numAdaptationBlocks, numData);
    strs1_block = ones(numAdaptationBlocks, numData);
    durs1_block = ones(numAdaptationBlocks, numData);
    strs2_block = ones(numAdaptationBlocks, numData);
    durs2_block = ones(numAdaptationBlocks, numData);
%     timeAxis = ones(numAdaptationBlocks, 1);
    
    for dataNum = 1:numData
        for blockInd = 1:numAdaptationBlocks
            thisAllInds = allBlockBorders(blockInd):(allBlockBorders(blockInd + 1) - 1);
            eRInds1 = logical(spikeER1(dataNum, thisAllInds));
            eRInds2 = logical(spikeER2(dataNum, thisAllInds));
            eRInds1 = thisAllInds(eRInds1);
            eRInds2 = thisAllInds(eRInds2);
            
            % For neuron 1 stimulations
            thisSR11 = spikeSR1(dataNum, eRInds1);
            thisSR21 = spikeSR2(dataNum, eRInds1);
            thisStrs1 = spikeSS(dataNum, eRInds1);
            thisDurs1 = spikeSD(dataNum, eRInds1);
            TP1_block(blockInd, dataNum) = mean(thisSR11);
            FA2_block(blockInd, dataNum) = mean(thisSR21);
            costs1_block(blockInd, dataNum) = adaptiveCostFunction(TP1_block(blockInd, dataNum), FA2_block(blockInd, dataNum), thisStrs1(1));
            
            % For neuron 2 stimulations
            thisSR12 = spikeSR1(dataNum, eRInds2);
            thisSR22 = spikeSR2(dataNum, eRInds2);
            thisStrs2 = spikeSS(dataNum, eRInds2);
            thisDurs2 = spikeSD(dataNum, eRInds2);
            TP2_block(blockInd, dataNum) = mean(thisSR22);
            FA1_block(blockInd, dataNum) = mean(thisSR12);
            costs2_block(blockInd, dataNum) = adaptiveCostFunction(TP2_block(blockInd, dataNum), FA1_block(blockInd, dataNum), thisStrs2(1));
            
            % Get general block data
%             timeAxis(blockInd) = adjustedStimTimesCol(round(mean(thisAllInds)));
            strs1_block(blockInd, dataNum) = thisStrs1(1);
            durs1_block(blockInd, dataNum) = thisDurs1(1); % TODO
            strs2_block(blockInd, dataNum) = thisStrs2(1);
            durs2_block(blockInd, dataNum) = thisDurs2(1);
        end
    end
    
    %% Perform the cost source analysis
    % Calculate the costs for each stimulus
    SACosts = difCostFn([TP1Avg ones(dataNum,1) TP1Avg], [FA2Avg FA2Avg zeros(dataNum,1)]);
    SBCosts = difCostFn([TP2Avg ones(dataNum,1) TP2Avg], [FA1Avg FA1Avg zeros(dataNum,1)]);
    
    SACostsOnline = difCostFn([TP1Avg ones(dataNum,1) TP1Avg], [FA2Avg FA2Avg zeros(dataNum,1)]);
    SBCostsOnline = difCostFn([TP2Avg ones(dataNum,1) TP2Avg], [FA1Avg FA1Avg zeros(dataNum,1)]);
    
    % Get the adjusted costs for each stimulation (SA/SB) in each sorting
    % (online/offline)
    costSourceCompareA = [(SACosts(:,1) - SACosts(:,2)) (SACosts(:,1) - SACosts(:,3))]; % A numData x 2 matrix that shows the amount by which the cost would decrease if the target neuron was activated perfectly (column 1) or the non-target neuron was silent perfectly (column 2)
    costSourceCompareB = [(SBCosts(:,1) - SBCosts(:,2)) (SBCosts(:,1) - SBCosts(:,3))];
    
    costSourceCompareOnlineA = [(SACostsOnline(:,1) - SACostsOnline(:,2)) (SACostsOnline(:,1) - SACostsOnline(:,3))]; % A numData x 2 matrix that shows the amount by which the cost would decrease if the target neuron was activated perfectly (column 1) or the non-target neuron was silent perfectly (column 2)
    costSourceCompareOnlineB = [(SBCostsOnline(:,1) - SBCostsOnline(:,2)) (SBCostsOnline(:,1) - SBCostsOnline(:,3))];
    
    targetProblemA = diff(costSourceCompareA')' < 0; % The data where problems with the target neuron (too LITTLE activity) were the main contributor to high costs in stimulation A
    nonTargetProblemA = diff(costSourceCompareA')' > 0; % The data where problems with the non-target neuron (too MUCH activity) were the main contributor to high costs in stimulation A
    targetProblemB = diff(costSourceCompareB')' < 0; % The data where problems with the target neuron (too LITTLE activity) were the main contributor to high costs in stimulation B
    nonTargetProblemB = diff(costSourceCompareB')' > 0; % The data where problems with the non-target neuron (too MUCH activity) were the main contributor to high costs in stimulation B
    
    targetProblemOnlineA = diff(costSourceCompareOnlineA')' < 0; % The data where problems with the target neuron (too LITTLE activity) were the main contributor to high costs in stimulation A
    nonTargetProblemOnlineA = diff(costSourceCompareOnlineA')' > 0; % The data where problems with the non-target neuron (too MUCH activity) were the main contributor to high costs in stimulation A
    targetProblemOnlineB = diff(costSourceCompareOnlineB')' < 0; % The data where problems with the target neuron (too LITTLE activity) were the main contributor to high costs in stimulation B
    nonTargetProblemOnlineB = diff(costSourceCompareOnlineB')' > 0; % The data where problems with the non-target neuron (too MUCH activity) were the main contributor to high costs in stimulation B
    
    sumTargetProblemA = sum(targetProblemA); % The total number of data where problems with the target neuron lead to high costs in A
    sumNonTargetProblemA = sum(nonTargetProblemA); % The total number of data where problems with the non-target neuron lead to high costs in A
    sumTargetProblemB = sum(targetProblemB); % The total number of data where problems with the target neuron lead to high costs in B
    sumNonTargetProblemB = sum(nonTargetProblemB); % The total number of data where problems with the non-target neuron lead to high costs in B
    
    sumTargetProblemSameOnlineA = sum(targetProblemA & targetProblemOnlineA); % The total number of data where problems with the target neuron lead to high costs in A for BOTH online AND offline sorting
    sumNonTargetProblemSameOnlineA = sum(nonTargetProblemA & nonTargetProblemOnlineA); % The total number of data where problems with the non-target neuron lead to high costs in A for BOTH online AND offline sorting
    sumTargetProblemSameOnlineB = sum(targetProblemB & targetProblemOnlineB); % The total number of data where problems with the target neuron lead to high costs in B for BOTH online AND offline sorting
    sumNonTargetProblemSameOnlineB = sum(nonTargetProblemB & nonTargetProblemOnlineB); % The total number of data where problems with the non-target neuron lead to high costs in B for BOTH online AND offline sorting
    
    % Get the response fraction difference for each stimulation (A attempts to
    % activate 1, B attempts to activate 2)
    rFDiffA = tPFrac1 - fAFrac2;
    rFDiffB = tPFrac2 - fAFrac1;
    
    difCostA = difCostFn(tPFrac1, fAFrac2);
    difCostB = difCostFn(tPFrac2, fAFrac1);
    
    sqeCostA = sqeCostFn(tPFrac1, fAFrac2);
    sqeCostB = sqeCostFn(tPFrac2, fAFrac1);
    
    % Make a total cost by combining the costs from both electrodes
    difCost = squaredMeanCostCombine(difCostA, difCostB);
    
    sqeCost = squaredMeanCostCombine(sqeCostA, sqeCostB);
    
    % Find the mean difference/cost over all sequences for each set
    rFDiffMeanA = nanmean(rFDiffA,1);
    rFDiffMeanB = nanmean(rFDiffB,1);
    
    difCostMeanTotalA = nanmean(difCostA,1);
    difCostMeanTotalB = nanmean(difCostB,1);
    difCostMeanPost = nanmean(difCost(:,doAdapt),1);
    
    sqeCostMeanTotalA = nanmean(sqeCostA,1);
    sqeCostMeanTotalB = nanmean(sqeCostB,1);
    %     sqeCostMeanPre = nanmean(sqeCost(:,~doAdapt),1);
    sqeCostMeanPost = nanmean(sqeCost(:,doAdapt),1);
    
    rFDiffSTEA = nanstd(rFDiffA,1)/sqrt(seqsPerSet); % Calculate the standard error for the performance according to the response fraction
    rFDiffSTEB = nanstd(rFDiffB,1)/sqrt(seqsPerSet);
    
    difCostMeanTotal = difCostMeanPost;
    sqeCostMeanTotal = sqeCostMeanPost;
    
    doAdaptGrouped = doAdapt;
    
    difCostSTE = nanstd(difCost,1)/sqrt(seqsPerSet);
    
    sqeCostSTE = nanstd(sqeCost,1)/sqrt(seqsPerSet);
    
%     % Find the number of units and pairs that pass minimum controllability
%     numPairsMinimumControllableTotal = sum((rFDiffSTEA > 0) & (rFDiffSTEB > 0));
%     numPairsMinimumControllableA = sum(rFDiffSTEA > 0);
%     numPairsMinimumControllableB = sum(rFDiffSTEB > 0);
    
    % Get the densities of each performance tier
    % First generate the performance tier levels
    rfDiffSorted = sort(min([rFDiffMeanA; rFDiffMeanB], [], 1), 'descend');
    % Declare that we want some fraction of pairs in the top tier, and
    % generate our tier list
    performanceTiers = [-Inf linspace(0, rfDiffSorted(round(numData/numPerformanceTiers)), numPerformanceTiers - 1)];
    
    % First, get the densities according to response fraction difference
    %     tierDensitiesPreAdapt = zeros(numISIs,numPerformanceTiers);
    tierDensitiesPostAdapt = zeros(numISIs,numPerformanceTiers);
    %     assignedPerformancePerSetPre = zeros(numISIs,numPreAdapt/2);
    assignedPerformancePerSet = zeros(numISIs,numPostAdapt);
    for tierNum = numPerformanceTiers:-1:1
        
        thisTierLowerPostA = reshape(performanceTiers(tierNum) + rFDiffSTEA(doAdaptGrouped), numISIs, []);
        thisTierLowerPostB = reshape(performanceTiers(tierNum) + rFDiffSTEB(doAdaptGrouped), numISIs, []);
        
        thisTierIndsPostA = thisTierLowerPostA <= reshape(rFDiffMeanA(doAdaptGrouped), numISIs, []) & ~assignedPerformancePerSet;
        thisTierIndsPostB = thisTierLowerPostB <= reshape(rFDiffMeanB(doAdaptGrouped), numISIs, []) & ~assignedPerformancePerSet;
        
        thisTierIndsPost = thisTierIndsPostA & thisTierIndsPostB;
        tierDensitiesPostAdapt(:, tierNum) = sum(thisTierIndsPost,2);
        assignedPerformancePerSet(thisTierIndsPost) = tierNum;
    end
    
    assignedPerformancePerSetPostGrouped = assignedPerformancePerSet;
    tierDensitiesPostAdaptGrouped = tierDensitiesPostAdapt;
    % Second, get the performance densities according to the cost functions
    
    % Trim both both mean cost vectors to only contain sets that only satisfy
    % the minimum control criteria
    goodPerformanceInds = assignedPerformancePerSet>1; % Find sets whos performance is better than 1 (or 1 partial)
    trimmedNumData = sum(goodPerformanceInds,2);
    
    difTierDensities = cell(numISIs, 1);
    sqeTierDensities = cell(numISIs, 1);
    for i = 1:numISIs
        thisGoodPerformanceInds = goodPerformanceInds(i,:);
        thisNumData = trimmedNumData(i,:);
        thisTrimmedDifCostMean = difCostMeanTotal(i, thisGoodPerformanceInds);
        thisTrimmedSqeCostMean = sqeCostMeanTotal(i, thisGoodPerformanceInds);
        
        % Get a single value for the overall performance of each date by
        % summing the costs of both stimulation types
        difCostBinSize = (max(thisTrimmedDifCostMean) - min(thisTrimmedDifCostMean))/numCostFunctionPerformanceTiers;
        sqeCostBinSize = (max(thisTrimmedSqeCostMean) - min(thisTrimmedSqeCostMean))/numCostFunctionPerformanceTiers;
        
        % Create tiers
        difPerformanceTiers = min(thisTrimmedDifCostMean) + difCostBinSize*(0:numCostFunctionPerformanceTiers);
        sqePerformanceTiers = min(thisTrimmedSqeCostMean) + sqeCostBinSize*(0:numCostFunctionPerformanceTiers);
        
        % Separate the data in to tiers based on the cost calculations
        difTierDensities{i} = zeros(1,numCostFunctionPerformanceTiers); % Should probably all sum to trimmedNumData?
        sqeTierDensities{i} = zeros(1,numCostFunctionPerformanceTiers);
        difDataAlreadyAssigned = zeros(1,thisNumData);
        sqeDataAlreadyAssigned = zeros(1,thisNumData);
        difPlotLabels = cell(1,numCostFunctionPerformanceTiers);
        sqePlotLabels = cell(1,numCostFunctionPerformanceTiers);
        for tierNum = numCostFunctionPerformanceTiers:-1:1
            tierDensitiesInd = tierNum; % The index in tierDensities which will hold the data for this tier
            
            % Lower boundary calculation for tier (and partial tier)
            % TODO: Include standard error calculation in here somehow?
            difTierLower = difPerformanceTiers(tierNum);% + trimmedDifCostSTE;
            difTierLowerPartial = difPerformanceTiers(tierNum);% + trimmedDifCostSTE;
            
            sqeTierLower = sqePerformanceTiers(tierNum);% + trimmedSqeCostSTE;
            sqeTierLowerPartial = sqePerformanceTiers(tierNum);% + trimmedSqeCostSTE;
            
            % Find which sets fall into this tier
            difTierInds = difTierLower <= thisTrimmedDifCostMean & ~difDataAlreadyAssigned;
            sqeTierInds = sqeTierLower <= thisTrimmedSqeCostMean & ~sqeDataAlreadyAssigned;
            
            difTierDensities{i}(tierDensitiesInd) = sum(difTierInds);
            difDataAlreadyAssigned(difTierInds) = tierDensitiesInd;
            sqeTierDensities{i}(tierDensitiesInd) = sum(sqeTierInds);
            sqeDataAlreadyAssigned(sqeTierInds) = tierDensitiesInd;
            
            % Create plot labels
            difPlotLabels{tierDensitiesInd} = sprintf('[%0.2f, %.2f)', difPerformanceTiers(tierNum), difPerformanceTiers(tierNum + 1));
            sqePlotLabels{tierDensitiesInd} = sprintf('[%0.2f, %.2f)', sqePerformanceTiers(tierNum), sqePerformanceTiers(tierNum + 1));
        end
    end
    
    % Get the average total cost for each response fraction tier
    
    averageDifCostPerTierPost = zeros(1,numPerformanceTiers);
    averageSqeCostPerTierPost = zeros(1,numPerformanceTiers);
    
    averageDifCostPerTierPostGrouped = zeros(numISIs,numPerformanceTiers);
    averageSqeCostPerTierPostGrouped = zeros(numISIs,numPerformanceTiers);
    plotLabels = cell(1,2*numPerformanceTiers);
    plotLabelsGrouped = cell(1,numPerformanceTiers);
    for tierNum = 1:numPerformanceTiers
        tierDensitiesInd = tierNum; % The index in tierDensities which will hold the data for this tier
        
        averageDifCostPerTierPost(tierDensitiesInd) = mean(difCostMeanPost(assignedPerformancePerSet == tierDensitiesInd)); % Take the mean of the total cost of dates whose performance falls within this tier
        
        averageSqeCostPerTierPost(tierDensitiesInd) = mean(sqeCostMeanPost(assignedPerformancePerSet == tierDensitiesInd)); % Take the mean of the total cost of dates whose performance falls within this tier
        
        % Get the averaged performance for each set, still split up by ISI
        for isiNum = 1:numISIs
            averageDifCostPerTierPostGrouped(isiNum, tierNum) = mean(difCostMeanPost((find(assignedPerformancePerSetPostGrouped(isiNum, :) == tierNum) - 1)*numISIs + isiNum));
            averageSqeCostPerTierPostGrouped(isiNum, tierNum) = mean(sqeCostMeanPost((find(assignedPerformancePerSetPostGrouped(isiNum, :) == tierNum) - 1)*numISIs + isiNum));
        end
        
        plotLabels{tierDensitiesInd} = sprintf('%d', tierNum - 1);
        plotLabelsGrouped{tierNum} = sprintf('%d', tierNum - 1);
    end
    
    xDataDifGroupedPost = 1:size(tierDensitiesPostAdaptGrouped, 2);
    xDataSqeGroupedPost = xDataDifGroupedPost;
    
    % Clustering parameters
    clusterNumISIs = numISIs;
    numDataPerISI = max(sum(uniqueISInds, 2)); % Max number of data points that exist for an ISI
    numLines = clusterNumISIs*numDataPerISI; % Number of lines for the excel document
    experimentNames = cell(numLines, 1);
    
    if useOfflineSort
        % Extract the desired clustering parameter
        allClustParams = nan(numDataPerISI, clusterNumISIs, numUnits(1));
        numCategories = 6; % Number of categories in the data table
    else
        numCategories = 4; % Number of categories in the data table
    end
    
    dataTable = nan(numLines, numCategories);
    dataTableInd = 1;
    
    % Choose which cost to use
    if useDiff
        totalCostA = difCostMeanTotalA;
        totalCostB = difCostMeanTotalB;
    else
        totalCostA = sqeCostMeanTotalA;
        totalCostB = sqeCostMeanTotalB;
    end
    
    for ISINum = 1:clusterNumISIs
        allThisISIInds = find(uniqueISInds(ISINum, :)); % Find the indices in all data that represent this ISI
        numThisISIData = length(allThisISIInds); % Number of data points that really exist for this ISI
        if useOfflineSort
            headers = {'Experiment Order', 'Data and Time', sprintf('%s N_A', clustParam), sprintf('%s N_B', clustParam), 'Total cost', 'N_A cost', 'N_B cost'}; % Define the headers for the excel document
        else
            headers = {'Experiment Order', 'Data and Time', 'Total cost', 'N_A cost', 'N_B cost'}; % Define the headers for the excel document
        end
        
        ISIStr = sprintf('_%d', uniqueIS(ISINum));
        
        for dataInd = 1:numThisISIData
            dataNum = allThisISIInds(dataInd);
            if useOfflineSort
                for neuronNum = 1:numUnits(1)
                    allClustParams(dataInd, ISINum, neuronNum) = clustStat{dataNum}(neuronNum).L_Ratio.(clustParam); % Extract the clustering parameters
                end
                dataTable(dataTableInd, :) = [dataInd squeeze(allClustParams(dataInd, ISINum, :)) difCostMeanTotal(ISINum,dataInd) totalCostA(ISINum,dataInd) totalCostB(ISINum,dataInd)]; % Organize the data needed to put into the excel document
            else
                dataTable(dataTableInd, :) = [dataInd difCostMeanTotal(ISINum,dataInd) totalCostA(ISINum,dataInd) totalCostB(ISINum,dataInd)]; % Organize the data needed to put into the excel document
            end
            
            experimentNames{dataTableInd} = sprintf('%s_%s%s', dateStrs{dataNum}, timeStrs{dataNum}, ISIStr);
            dataTableInd = dataTableInd + 1;
        end
    end
    dataTable = dataTable(1:(dataTableInd - 1),:);
    experimentNames = experimentNames(1:(dataTableInd - 1));
    finalTable = [headers;[experimentNames, num2cell(dataTable)];[{'';'Mean';'STD'}, [cell(1, numCategories);num2cell([nanmean(dataTable,1);nanstd(dataTable,[],1)])]]];
    
    doAdaptFirstISIInds = 1:numISIs:numData;
    doAdaptFirstISI = doAdapt(doAdaptFirstISIInds);
    if useOfflineSort
        clustParams = allClustParams(doAdaptFirstISI, :, :);
        clustParamsInclude = all(clustParams < lRatioLimit, 3);
    end
    
    
    dateStrsInclude = dateStrs(doAdapt);
    timeStrsInclude = timeStrs(doAdapt);
    trialNumsInclude = trialNums(doAdapt);
    %     end
    
    % Create x axes for both cost per tier plots
    
    xDataDifPost = 1:(numPerformanceTiers); % Double the number of performance tiers, because we are including "partial" performance tiers
    xDataSqePost = xDataDifPost;
    
    
    % Remove NaNs from both x and y data
    
    difRemovedPost = isnan(averageDifCostPerTierPost);
    sqeRemovedPost = isnan(averageSqeCostPerTierPost);
    averageDifCostPerTierPost(difRemovedPost) = [];
    averageSqeCostPerTierPost(sqeRemovedPost) = [];
    xDataDifPost(difRemovedPost) = [];
    xDataSqePost(sqeRemovedPost) = [];
    %
    difRemovedPostGrouped = isnan(averageDifCostPerTierPostGrouped);
    sqeRemovedPostGrouped = isnan(averageSqeCostPerTierPostGrouped);
    averageDifCostPerTierPostGrouped(difRemovedPostGrouped) = [];
    averageSqeCostPerTierPostGrouped(sqeRemovedPostGrouped) = [];
    xDataDifGroupedPost(difRemovedPostGrouped) = [];
    xDataSqeGroupedPost(sqeRemovedPostGrouped) = [];
end

%% Do statistical
if doStatisticalTests
    % First, test for how many pairs was the control significantly
    % effective (and/or how confident are we)?
    rho = 0; % Correlation between CQA and CQB; Used to calculate the distribution of the min of the two stimulations' CQ's
    hdrMass = .95; % The probability mass that must exist in the high density region
    p = .05; % The p value that we are looking to find a confidence interval for
    
    alpha = 1 - p;
    %% Overall control
    % First test if the pair was overall controllable
    
    % Make a couple of useful function shortcuts
    sCDF = @(x) normcdf(x, 0, 1); % Standard normal cdf
    sPDF = @(x) normpdf(x, 0, 1); % Standard normal pdf
    
    % First calculate the means and standard error of the means of the data
    % % TODO: FIX THIS
    % Create matrices that will first accept all of the response fraction
    % differences, indexed by stimulation type (I know, it's annoying...)
    aStimMat = zeros(numData, numAllStims/2); % Half of all stimulations will be A
    bStimMat = zeros(numData, numAllStims/2);
    
    for dataNum = 1:numData
        aStimMat(dataNum, :) = spikeSR1(dataNum, spikeER1(dataNum, :)) - spikeSR2(dataNum, spikeER1(dataNum, :)); % Calculate the response fraction differences per stim, and fill into aMat
        bStimMat(dataNum, :) = spikeSR2(dataNum, spikeER2(dataNum, :)) - spikeSR1(dataNum, spikeER2(dataNum, :)); % Calculate the response fraction differences per stim, and fill into aMat
    end
    
    CQATrue = mean(aStimMat, 2); % numData x 1
    CQBTrue = mean(bStimMat, 2); % numData x 1
    CQTrue = min(CQATrue, CQBTrue); % numData x 1
    
    muA = CQATrue; % Reponse fraction difference A
    muB = CQBTrue; % Reponse fraction difference B
    
    sigmaA = std(aStimMat, 0, 2);
    sigmaB = std(bStimMat, 0, 2);
    steA = sigmaA./sqrt(numAllStims);
    steB = sigmaB./sqrt(numAllStims);
    estPDistMin = @(x, i) pDistMin(x, muA(i), muB(i), steA(i), steB(i), rho); % Get the value in the probability density at x, of the i-th data
%     estPDistMin = @(x, i) pDistMin(x, muA(i), muB(i), sigmaA(i), sigmaB(i), rho); % Get the value in the probability density at x, of the i-th data
    
    % Calculate the distributions of the minimums of the control quality
    % metric (the response fraction difference)
    theta = sqrt(steA.^2 + steB.^2 - 2*rho*steA.*steB); % Used in other expressions
    cqEstMean = muA.*sCDF((muB - muA)./theta) + muB.*sCDF((muA - muB)./theta) - theta.*sPDF((muB - muA)./theta); % The mean of the estimate
    cqEstUncVar = (steA.^2 + muA.^2).*sCDF((muB - muA)./theta) + (steB.^2 + muB.^2).*sCDF((muA - muB)./theta) - (muA + muB).*theta.*sPDF((muB - muA)./theta); % Uncentered variance of the estimate
    cqEstSte = sqrt(cqEstUncVar - cqEstMean.^2);

    %% Determine if we have acheived good control
    %     ciRadius = zeros(numData, 1); % The confidence interval radii
    %     for dataInd = 1:numData
    %         ciRadius(dataInd) = fminunc(@(int) (integral(@(x) estPDistMin(x, dataInd), cqEstMean(dataInd) - int, cqEstMean(dataInd) + int) - alpha).^2, 2*cqEstSte(dataInd)); % Find the radius for our confidence interval (whatever the hell that's called)
    %     end
    %     confidenceIntervals = cqEstMean + ciRadius.*repmat([-1 1], numAllStims, 1); % Each row is the confidence interval for each data set
    
    % Determine if we have acheived good control by finding the high
    % density region
    thisHdrBounds = zeros(numData, 2); % The high density region boundaries
    for dataInd = 1:numData
        initYVal = .8*estPDistMin(cqEstMean(dataInd), dataInd);
        hdrYVal = fminunc(@(yVal) (integral(@(x) estPDistMin(x, dataInd), fmincon(@(xMin) (estPDistMin(xMin, dataInd) - yVal)^2, cqEstMean(dataInd)-3*cqEstSte(dataInd), [], [], [], [], [], cqEstMean(dataInd)), fmincon(@(xMax) (estPDistMin(xMax, dataInd) - yVal)^2, cqEstMean(dataInd)+3*cqEstSte(dataInd), [], [], [], [], cqEstMean(dataInd))) - hdrMass).^2, initYVal);
        thisHdrBounds(dataInd, :) = [fmincon(@(xMin) (estPDistMin(xMin, dataInd) - hdrYVal)^2, cqEstMean(dataInd)-3*cqEstSte(dataInd), [], [], [], [], [], cqEstMean(dataInd)) fmincon(@(xMax) (estPDistMin(xMax, dataInd) - hdrYVal)^2, cqEstMean(dataInd)+3*cqEstSte(dataInd), [], [], [], [], cqEstMean(dataInd))]; % Use the y value to calculate the hdr bounds
    end
    
    %     conclusiveControl = (abs(cqEstMean) - abs(ciRadius)) > 0; % For each data set, is the confidence interval test conclusive (does the interval NOT contain 0?)
    conclusiveTotalControl = (0 < thisHdrBounds(:, 1)) | (thisHdrBounds(:, 2) < 0);
    goodTotalControl = cqEstMean > 0; % Is the mean above 0 (i.e. good control, whether or not it is conclusive)
    goodConclusiveTotalControl = conclusiveTotalControl & goodTotalControl; % If we conclusively have good control
    
    conclusiveAControl = ((CQATrue + 2*steA) < 0) | (0 < (CQATrue - 2*steA)); % The confidence interval is defined as a radius of 2*ste around the mean (http://www.mas.ncl.ac.uk/~njnsm/medfac/docs/se&ci.pdf)
    conclusiveBControl = ((CQBTrue + 2*steB) < 0) | (0 < (CQBTrue - 2*steB));
    goodAControl = CQATrue > 0;
    goodBControl = CQBTrue > 0;
    goodConclusiveAControl = conclusiveAControl & goodAControl;
    goodConclusiveBControl = conclusiveBControl & goodBControl;
    goodConclusiveOneWayControl = goodConclusiveAControl | goodConclusiveBControl;
    
    % Create the plot
    % Generate the y values
    yHDRVals = [1:numData;1:numData];
    yCQMeanVals = 1:numData;
    
    % Generate the x values
    [cqMeansSorted, cqMeanSortInds] = sort(CQTrue);
    hdrBoundsSorted = thisHdrBounds(cqMeanSortInds, :)'; % 2 x numData
    
    if doStatisticalTestingPlot
        hStatsControl = figure;
        a = axes;
        hold(a, 'on');
        plot(a, cqMeansSorted, yCQMeanVals, 'k.', 'markersize', 20); % First plot all of the mean values
        plot(a, hdrBoundsSorted, yHDRVals, 'k-', 'linewidth', 3);
        
        yLims = ylim;
        plot(a, [0 0], yLims, 'k:');
        ylim([0 (numData + 1)]);
        xlim([-1 1]);
        
        xlabel('Control Quality CQ');
        ylabel('Pair Number');
        title('Summary of Mean Control Quality');
    end
    
    % Calculate the p value of the confidence interval that contains 0
    % TODO: UPDATE THIS
    %     ciBounds = sort([0 2*cqEstMean], 2); % Find the boundaries for the
    %     thisAlphaVal = zeros(numData, 1); % The alpha values of the confidence interval that contains 0
    %     for dataInd = 1:numData
    %         thisAlphaVal(dataInd) = integral(@(x) estPDistMin(x, dataInd), ciBounds(dataInd, 1), ciBounds(dataInd, 2));
    %     end
    %     thisPVal = 1 - thisAlphaVal; % The confidence with which we can say we have acheived control
    
    % TODO: NOT useful
    % Finally, calculate the p value of CQ being above 0
%     thisPVal = zeros(numData, 1);
%     for dataInd = 1:numData
%         thisPVal = integral(@(x) estPDistMin(x, dataInd), -Inf, 0);
%     end
    
    %% Differentiability of stimuli
    % Next, make sure that the stimuli actually do different things.  To do
    % this, we are going to shuffle the stimulation A/B labels and see if
    % the data looks significantly different.  Shuffle will be done by
    % stimulation (as opposed to by sequence).
    
    % Perform the shuffle
    erUnsorted = repelem([false; true], numAllStims/2, 1); % A vector of an equal number of 0's and 1's, numAllStims in total, that will be rearranged to generate a new ER vector for each shuffle
    
    CQAShuffled = zeros(numData, numShuffles);
    CQBShuffled = zeros(numData, numShuffles);
    for shuffInd = 1:numShuffles
        thisSpikeER1 = erUnsorted(randperm(numAllStims)); % Generate this shuffle's ER
%         thisSpikeER1 = repmat(thisSpikeER1(:), 1, numData); % (numAllStims x numData) Repeat the matrix for each data set
        thisSpikeER2 = ~thisSpikeER1; % numAllStims x numData
        
        % Find the response fraction difference per stim 
        thisAStimMat = spikeSR1(:, thisSpikeER1) - spikeSR2(:, thisSpikeER1); % (numData x (numAllStims/2) ) For every data set, apply the new ER vector
        thisBStimMat = spikeSR2(:, thisSpikeER2) - spikeSR1(:, thisSpikeER2);
        
        % Take the mean of the RFD's, and put into storage matrix
        CQAShuffled(:, shuffInd) = mean(thisAStimMat, 2);
        CQBShuffled(:, shuffInd) = mean(thisBStimMat, 2);
    end
    CQShuffled = min(cat(3, CQAShuffled, CQBShuffled), [], 3); % (numData x numShuffles) For each shuffle of each data point, take the minimum between CQA and CQB
    
    cqMean = mean(CQShuffled, 2); % (numData x 1)
    cqStd = std(CQShuffled, 0, 2); % (numData x 1)
    
    % Z-score our sample CQ according to this distribution
    sampleZScore = (CQTrue - cqMean)./cqStd; % (numData x 1)
    stimDifferentiabilityPVal = 2*(1 - normcdf(abs(sampleZScore)));
    
%     hdrBounds = nan(numData, 2); % The boundaries for each pair's high density region
%     conclusive = false(numData, 1); % Do the stim labels inform the spike results?
%     allHDRChecks = nan(numData, 1); % A way to check the validity of each HDR, this value will be the fraction of the probability mass between the defined HDR bounds (should be equal to hdrMassFrac)
%     for dataNum = 1:numData
%         [binCounts, binEdges] = histcounts(CQShuffled(dataNum, :)); % Get a histogram of CQShuffled
%         binMids = mean([binEdges(1:(end - 1));binEdges(2:end)], 1);
%         smoothedBinCounts = smooth(binCounts, round(length(binCounts)/10)); % Smooth out the histogram
%         cqPostInterp = @(x) interp1(binMids, smoothedBinCounts, x);
%         yInit = .8*cqPostInterp(cqMean(dataNum));
%         totalMass = integral(cqPostInterp, binMids(1), binMids(end));
%         intAtY = @(yVal, xStartOffset) integral(cqPostInterp, fmincon(@(xMin) (cqPostInterp(xMin) - yVal)^2, cqMean(dataNum) - xStartOffset, [], [], [], [], [], cqMean(dataNum)), fmincon(@(xMax) (cqPostInterp(xMax) - yVal)^2, cqMean(dataNum) + xStartOffset, [], [], [], [], cqMean(dataNum)));
%         
%         % Calculate the high density region
%         cqMeanTestOffsets = cqStd(dataNum)*[.01 .1 .2 .5 1 1.5];
%         numOffsets = length(cqMeanTestOffsets);
%         tolErr = 1e-4;
%         haveHDR = false;
%         thisDataChecks = nan(numOffsets, 1);
%         for offsetNum = 1:numOffsets
%             cqDistYVal = fmincon(@(yVal)(intAtY(yVal, cqMeanTestOffsets(offsetNum)) - totalMass*hdrMassFrac)^2, yInit, [], [], [], [], 0, cqPostInterp(cqMean(dataNum)));
%             thisHdrBounds = [fmincon(@(xMin) (cqPostInterp(xMin) - cqDistYVal)^2, cqMean(dataNum) - cqMeanTestOffsets(offsetNum), [], [], [], [], [], cqMean(dataNum)), fmincon(@(xMax) (cqPostInterp(xMax) - cqDistYVal)^2, cqMean(dataNum) + cqMeanTestOffsets(offsetNum), [], [], [], [], cqMean(dataNum))];
%             check = integral(cqPostInterp, thisHdrBounds(1), thisHdrBounds(2))/totalMass; % This should be approximately equal to hdrMassFrac
%             thisDataChecks(offsetNum) = check;
%             if abs(check - hdrMassFrac) < tolErr
%                 haveHDR = true;
%                 continue;
%             end
%         end
%         
%         if haveHDR
%             conclusive(numData) = (CQTrue(dataNum) < thisHdrBounds(1)) | (thisHdrBounds(2) < CQTrue(dataNum));
%             hdrBounds(numData, :) = thisHdrBounds;
%             allHDRChecks(numData) = check;
%         else
%             %             warning('Was not able to find a good HDR boundary!');
%             conclusive(numData) = stimDifferentiabilityPVal(dataNum) < .025; % P for a two-tailed test
%         end
%         
%         %         if conclusive
%         %             hasInfoStr = 'stim labels inform results';
%         %         else
%         %             hasInfoStr = 'stim labels do not inform results';
%         %         end
%     end
    
    % TODO: Finish this test (perform similar confidence interval test as
    % above?)
    if doStatisticalTestingPlot
       hStatsDiff = figure;
       a = axes;
       histogram(a, sampleZScore);
       xlabel('Z-Score');
       ylabel('Bin Count');
       title('Stimulation Differentiability Z-Scores');
    end
end

%% Analyze epoch-border
if doCostFunctionPlot
    allNumEpochs = cellfun(@(x) length(x), epochStartInds);
    minNumEpochs = min(allNumEpochs); % The minimum number of epochs that a recalibration experiment has
    trimmedEpochStartInds = cellfun(@(x) trimEpochStartIndHelperFun(x, minNumEpochs, numSequences(1)), epochStartInds, 'uniformoutput', false);
    trimmedEpochStartInds = cell2mat(trimmedEpochStartInds);
    
    prevCostPurAll1 = [];
    prevCostPurAll2 = [];
    prevCostPurAllS = [];
    postCostPurAll1 = [];
    postCostPurAll2 = [];
    postCostPurAllS = [];
    
    
    prevCostMixAll1 = [];
    prevCostMixAll2 = [];
    prevCostMixAllS = [];
    postCostMixAll1 = [];
    postCostMixAll2 = [];
    postCostMixAllS = [];
    
    for epochBoundInd = 1:(minNumEpochs - 1)
        prevStartInds = (trimmedEpochStartInds(epochBoundInd,:) - 1)*numStimsPerSequence + 1;
        thisStartInds = (trimmedEpochStartInds(epochBoundInd + 1,:) - 1)*numStimsPerSequence + 1;
        postStartInds = (trimmedEpochStartInds(epochBoundInd + 2,:) - 1)*numStimsPerSequence + 1;
        prevEpochSizes = floor((thisStartInds - prevStartInds)/2);
        postEpochSizes = ceil((postStartInds - thisStartInds)/2);
        
        prevCostPur1 = nan(1, numData);
        prevCostPur2 = nan(1, numData);
        prevCostPurS = nan(1, numData);
        postCostPur1 = nan(1, numData);
        postCostPur2 = nan(1, numData);
        postCostPurS = nan(1, numData);
        
        prevCostMix1 = nan(1, numData);
        prevCostMix2 = nan(1, numData);
        prevCostMixS = nan(1, numData);
        postCostMix1 = nan(1, numData);
        postCostMix2 = nan(1, numData);
        postCostMixS = nan(1, numData);
        
        for dataNum = 1:numData
            prevAllInds = prevStartInds(dataNum):(thisStartInds(dataNum) - 1);
            postAllInds = thisStartInds(dataNum):(postStartInds(dataNum) - 1);
            prevAllInds = prevAllInds((end - prevEpochSizes(dataNum) + 1):end); % Get the indices for ALL stimulations in the second half of the previous epoch
            postAllInds = postAllInds(1:(postEpochSizes(dataNum) - 1)); % Get the indices for ALL stimulations in the first half of the next epoch
            
            % Prepare expected results
            prevER1 = logical(spikeER1(dataNum,prevAllInds));
            postER1 = logical(spikeER1(dataNum,postAllInds));
            
            prevER2 = ~prevER1;
            postER2 = ~postER1;
            
            % Prepare stim results
            prevStimResults1 = spikeSR1(dataNum,prevAllInds);
            prevStimResults2 = spikeSR2(dataNum,prevAllInds);
            postStimResults1 = spikeSR1(dataNum,postAllInds);
            postStimResults2 = spikeSR2(dataNum,postAllInds);
            
            % Get the allowed G's (is this a "pure" comparison, or "mixed"?)
            prevAllG = spikeSS(dataNum, prevAllInds); % Previous epoch's G's
            postAllG = spikeSS(dataNum, postAllInds); % Post epoch's G's
            
            if ~doAdapt(dataNum)
                prevPurG1 = prevAllG(find(prevER1, 1));
                prevPurG2 = prevAllG(find(prevER2, 1));
                postPurG1 = postAllG(find(postER1, 1));
                postPurG2 = postAllG(find(postER2, 1));
                
                prevMixG1 = prevPurG1;
                prevMixG2 = prevPurG2;
                postMixG1 = postPurG1;
                postMixG2 = postPurG2;
            else
                if isempty(adaptiveInds1{dataNum})
                    allPrevMixG1 = prevAllG(prevER1); % Get all instances of all strengths for this unit in this epoch
                    allPrevMixG2 = prevAllG(prevER2); % Get all instances of all strengths for this unit in this epoch
                    allPostMixG1 = postAllG(postER1); % Get all instances of all strengths for this unit in this epoch
                    allPostMixG2 = postAllG(postER2); % Get all instances of all strengths for this unit in this epoch
                    
                    prevMixG1 = unique(allPrevMixG1, 'sorted'); % Get only one instance of each strength in this unit of this epoch (can only be current, search up, and search down)
                    prevMixG2 = unique(allPrevMixG2, 'sorted'); % Get only one instance of each strength in this unit of this epoch (can only be current, search up, and search down)
                    postMixG1 = unique(allPostMixG1, 'sorted'); % Get only one instance of each strength in this unit of this epoch (can only be current, search up, and search down)
                    postMixG2 = unique(allPostMixG2, 'sorted'); % Get only one instance of each strength in this unit of this epoch (can only be current, search up, and search down)
                    
                    if (length(prevMixG1) ~= 3 || length(prevMixG2) ~= 3)
                        if (epochBoundInd == 1)
                            fprintf('Skipping data number %d, badly adjusted ER\n', dataNum);
                        end
                        continue;
                    end
                    prevPurG1 = prevMixG1(2); % Get the middle value in the mixed (which must be the "current" point, because it will be in the middle of the sorted mixed vector)
                    prevPurG2 = prevMixG2(2); % Get the middle value in the mixed (which must be the "current" point, because it will be in the middle of the sorted mixed vector)
                    postPurG1 = postMixG1(2); % Get the middle value in the mixed (which must be the "current" point, because it will be in the middle of the sorted mixed vector)
                    postPurG2 = postMixG2(2); % Get the middle value in the mixed (which must be the "current" point, because it will be in the middle of the sorted mixed vector)
                else
                    a1 = adaptiveInds1{dataNum};
                    a2 = adaptiveInds2{dataNum};
                    prevPurG1 = a1.searchPoints(epochBoundInd, 1, 1);
                    prevPurG2 = a2.searchPoints(epochBoundInd, 1, 1);
                    postPurG1 = a1.searchPoints(epochBoundInd + 1, 1, 1);
                    postPurG2 = a2.searchPoints(epochBoundInd + 1, 1, 1);
                    
                    prevMixG1 = [squeeze(a1.searchPoints(epochBoundInd, 1, :));squeeze(a1.escapePoints(epochBoundInd, 1,:))];
                    prevMixG2 = [squeeze(a2.searchPoints(epochBoundInd, 1, :));squeeze(a2.escapePoints(epochBoundInd, 1,:))];
                    postMixG1 = [squeeze(a1.searchPoints(epochBoundInd + 1, 1, :));squeeze(a1.escapePoints(epochBoundInd + 1, 1,:))];
                    postMixG2 = [squeeze(a2.searchPoints(epochBoundInd + 1, 1, :));squeeze(a2.escapePoints(epochBoundInd + 1, 1,:))];
                end
            end
            prevPurGInds1 = prevAllG == prevPurG1;
            prevPurGInds2 = prevAllG == prevPurG2;
            postPurGInds1 = postAllG == postPurG1;
            postPurGInds2 = postAllG == postPurG2;
            
            prevMixGInds1 = ismember(prevAllG, prevMixG1);
            prevMixGInds2 = ismember(prevAllG, prevMixG2);
            postMixGInds1 = ismember(postAllG, postMixG1);
            postMixGInds2 = ismember(postAllG, postMixG2);
            
            
            % Calculate the response fractions
            prevTP1Log = prevStimResults1 & prevER1;
            prevFA1Log = prevStimResults1 & prevER2;
            prevTP2Log = prevStimResults2 & prevER2;
            prevFA2Log = prevStimResults2 & prevER1;
            
            postTP1Log = postStimResults1 & postER1;
            postFA1Log = postStimResults1 & postER2;
            postTP2Log = postStimResults2 & postER2;
            postFA2Log = postStimResults2 & postER1;
            
            % Average the response fractions
            prevPurTP1 = sum(prevTP1Log(prevPurGInds1))/sum(prevER1(prevPurGInds1));
            prevPurFA1 = sum(prevFA1Log(prevPurGInds2))/sum(prevER2(prevPurGInds2));
            prevPurTP2 = sum(prevTP2Log(prevPurGInds2))/sum(prevER2(prevPurGInds2));
            prevPurFA2 = sum(prevFA2Log(prevPurGInds1))/sum(prevER1(prevPurGInds1));
            
            prevMixTP1 = sum(prevTP1Log(prevMixGInds1))/sum(prevER1(prevMixGInds1));
            prevMixFA1 = sum(prevFA1Log(prevMixGInds2))/sum(prevER2(prevMixGInds2));
            prevMixTP2 = sum(prevTP2Log(prevMixGInds2))/sum(prevER2(prevMixGInds2));
            prevMixFA2 = sum(prevFA2Log(prevMixGInds1))/sum(prevER1(prevMixGInds1));
            
            postPurTP1 = sum(postTP1Log(postPurGInds1))/sum(postER1(postPurGInds1));
            postPurFA1 = sum(postFA1Log(postPurGInds2))/sum(postER2(postPurGInds2));
            postPurTP2 = sum(postTP2Log(postPurGInds2))/sum(postER2(postPurGInds2));
            postPurFA2 = sum(postFA2Log(postPurGInds1))/sum(postER1(postPurGInds1));
            
            postMixTP1 = sum(postTP1Log(postMixGInds1))/sum(postER1(postMixGInds1));
            postMixFA1 = sum(postFA1Log(postMixGInds2))/sum(postER2(postMixGInds2));
            postMixTP2 = sum(postTP2Log(postMixGInds2))/sum(postER2(postMixGInds2));
            postMixFA2 = sum(postFA2Log(postMixGInds1))/sum(postER1(postMixGInds1));
            
            % Calculate the costs
            %         prevCost1(dataNum) = sqeCostFnSc(prevTP1, prevFA2);
            %         prevCost2(dataNum) = sqeCostFnSc(prevTP2, prevFA1);
            %
            %         postCost1(dataNum) = sqeCostFnSc(postTP1, postFA2);
            %         postCost2(dataNum) = sqeCostFnSc(postTP2, postFA1);
            
            %         prevCost1(dataNum) = difCostFnSc(prevTP1, prevFA2);
            %         prevCost2(dataNum) = difCostFnSc(prevTP2, prevFA1);
            %
            %         postCost1(dataNum) = difCostFnSc(postTP1, postFA2);
            %         postCost2(dataNum) = difCostFnSc(postTP2, postFA1);
            
            %             prevCostPur1(dataNum) = difCostFn(prevPurTP1, prevPurFA2);
            %             prevCostPur2(dataNum) = difCostFn(prevPurTP2, prevPurFA1);
            %
            %             postCostPur1(dataNum) = difCostFn(postPurTP1, postPurFA2);
            %             postCostPur2(dataNum) = difCostFn(postPurTP2, postPurFA1);
            
            prevCostPur1(dataNum) = sqeCostFn(prevPurTP1, prevPurFA2);
            prevCostPur2(dataNum) = sqeCostFn(prevPurTP2, prevPurFA1);
            
            postCostPur1(dataNum) = sqeCostFn(postPurTP1, postPurFA2);
            postCostPur2(dataNum) = sqeCostFn(postPurTP2, postPurFA1);
            
            prevCostPurS(dataNum) = (prevCostPur1(dataNum) + prevCostPur2(dataNum))/2;
            postCostPurS(dataNum) = (postCostPur1(dataNum) + postCostPur2(dataNum))/2;
            
            %             prevCostMix1(dataNum) = difCostFn(prevMixTP1, prevMixFA2);
            %             prevCostMix2(dataNum) = difCostFn(prevMixTP2, prevMixFA1);
            %
            %             postCostMix1(dataNum) = difCostFn(postMixTP1, postMixFA2);
            %             postCostMix2(dataNum) = difCostFn(postMixTP2, postMixFA1);
            
            prevCostMix1(dataNum) = sqeCostFn(prevMixTP1, prevMixFA2);
            prevCostMix2(dataNum) = sqeCostFn(prevMixTP2, prevMixFA1);
            
            postCostMix1(dataNum) = sqeCostFn(postMixTP1, postMixFA2);
            postCostMix2(dataNum) = sqeCostFn(postMixTP2, postMixFA1);
            
            prevCostMixS(dataNum) = (prevCostMix1(dataNum) + prevCostMix2(dataNum))/2;
            postCostMixS(dataNum) = (postCostMix1(dataNum) + postCostMix2(dataNum))/2;
        end
        
        %     prevCostAll1 = [prevCostAll1 prevCost1];
        %     prevCostAll2 = [prevCostAll2 prevCost2];
        %
        %     postCostAll1 = [postCostAll1 postCost1];
        %     postCostAll2 = [postCostAll2 postCost2];
        
        prevCostPurAll1 = [prevCostPurAll1 prevCostPur1];
        prevCostPurAll2 = [prevCostPurAll2 prevCostPur2];
        postCostPurAll1 = [postCostPurAll1 postCostPur1];
        postCostPurAll2 = [postCostPurAll2 postCostPur2];
        prevCostPurAllS = [prevCostPurAllS prevCostPurS];
        postCostPurAllS = [postCostPurAllS postCostPurS];
        
        prevCostMixAll1 = [prevCostMixAll1 prevCostMix1];
        prevCostMixAll2 = [prevCostMixAll2 prevCostMix2];
        postCostMixAll1 = [postCostMixAll1 postCostMix1];
        postCostMixAll2 = [postCostMixAll2 postCostMix2];
        prevCostMixAllS = [prevCostMixAllS prevCostMixS];
        postCostMixAllS = [postCostMixAllS postCostMixS];
    end
    plotWidth = .25;
    xAxisBounds = (1:(minNumEpochs - 1)) + [-plotWidth; plotWidth];
    xAxisAllCosts = (1:minNumEpochs)';
    
    %% Get cost of each epoch
    allCostPurAll1 = [];
    allCostPurAll2 = [];
    allCostMixAll1 = [];
    allCostMixAll2 = [];
    for epochInd = 1:minNumEpochs
        thisStartInds = (trimmedEpochStartInds(epochInd,:) - 1)*numStimsPerSequence + 1;
        thisEndInds = (trimmedEpochStartInds(epochInd + 1,:) - 1)*numStimsPerSequence;
        
        thisCostPur1 = nan(1, numData);
        thisCostPur2 = nan(1, numData);
        thisCostMix1 = nan(1, numData);
        thisCostMix2 = nan(1, numData);
        for dataNum = 1:numData
            thisAllInds = thisStartInds(dataNum):thisEndInds(dataNum);
            thisER1 = logical(spikeER1(dataNum, thisAllInds));
            thisER2 = ~thisER1;
            
            thisSR1 = spikeSR1(dataNum, thisAllInds);
            thisSR2 = spikeSR2(dataNum, thisAllInds);
            
            % Get the allowed G's (is the a "pure" comparison, or "mixed"?)
            allG = spikeSS(dataNum, thisAllInds); % Previous epoch's G's
            
            if ~doAdapt(dataNum)
                purG1 = allG(find(thisER1, 1));
                purG2 = allG(find(thisER2, 1));
                
                mixG1 = purG1;
                mixG2 = purG2;
            else
                if isempty(adaptiveInds1{dataNum})
                    allMixG1 = allG(thisER1); % Get all instances of all strengths for this unit in this epoch
                    allMixG2 = allG(thisER2); % Get all instances of all strengths for this unit in this epoch
                    
                    mixG1 = unique(allMixG1, 'sorted'); % Get only one instance of each strength in this unit of this epoch (can only be current, search up, and search down)
                    mixG2 = unique(allMixG2, 'sorted'); % Get only one instance of each strength in this unit of this epoch (can only be current, search up, and search down)
                    
                    if (length(mixG1) ~= 3 || length(mixG2) ~= 3)
                        if (epochInd == 1)
                            fprintf('Skipping data number %d, badly adjusted ER\n', dataNum);
                        end
                        continue;
                    end
                    
                    purG1 = mixG1(2); % Get the middle value in the mixed (which must be the "current" point, because it will be in the middle of the sorted mixed vector)
                    purG2 = mixG2(2); % Get the middle value in the mixed (which must be the "current" point, because it will be in the middle of the sorted mixed vector)
                else
                    a1 = adaptiveInds1{dataNum};
                    a2 = adaptiveInds2{dataNum};
                    purG1 = a1.searchPoints(epochInd, 1, 1);
                    purG2 = a2.searchPoints(epochInd, 1, 1);
                    
                    mixG1 = [squeeze(a1.searchPoints(epochInd, 1, :));squeeze(a1.escapePoints(epochInd, 1,:))];
                    mixG2 = [squeeze(a2.searchPoints(epochInd, 1, :));squeeze(a2.escapePoints(epochInd, 1,:))];
                end
            end
            purGInds1 = allG == purG1;
            purGInds2 = allG == purG2;
            
            mixGInds1 = ismember(allG, mixG1);
            mixGInds2 = ismember(allG, mixG2);
            
            
            thisTP1Log = thisSR1 & thisER1;
            thisFA1Log = thisSR1 & thisER2;
            thisTP2Log = thisSR2 & thisER2;
            thisFA2Log = thisSR2 & thisER1;
            
            thisPurTP1 = sum(thisTP1Log(purGInds1))/sum(thisER1(purGInds1));
            thisPurFA1 = sum(thisFA1Log(purGInds2))/sum(thisER2(purGInds2));
            thisPurTP2 = sum(thisTP2Log(purGInds2))/sum(thisER2(purGInds2));
            thisPurFA2 = sum(thisFA2Log(purGInds1))/sum(thisER1(purGInds1));
            
            thisMixTP1 = sum(thisTP1Log(mixGInds1))/sum(thisER1(mixGInds1));
            thisMixFA1 = sum(thisFA1Log(mixGInds2))/sum(thisER2(mixGInds2));
            thisMixTP2 = sum(thisTP2Log(mixGInds2))/sum(thisER2(mixGInds2));
            thisMixFA2 = sum(thisFA2Log(mixGInds1))/sum(thisER1(mixGInds1));
            
            thisCostPur1(dataNum) = sqeCostFn(thisPurTP1, thisPurFA2);
            thisCostPur2(dataNum) = sqeCostFn(thisPurTP2, thisPurFA1);
            
            thisCostMix1(dataNum) = sqeCostFn(thisMixTP1, thisMixFA2);
            thisCostMix2(dataNum) = sqeCostFn(thisMixTP2, thisMixFA1);
        end
        
        allCostPurAll1 = [allCostPurAll1; thisCostPur1];
        allCostPurAll2 = [allCostPurAll2; thisCostPur2];
        allCostMixAll1 = [allCostMixAll1; thisCostMix1];
        allCostMixAll2 = [allCostMixAll2; thisCostMix2];
    end
    
    % Separate cost data into pre and post adapt
    doAdaptAllEpochBounds = repmat(doAdapt,1,minNumEpochs - 1);
    includeISIEpochBounds = repmat(includedISIInds, 1, minNumEpochs - 1);
    prevCostPurPreAdapt1 = prevCostPurAll1(~doAdaptAllEpochBounds & includeISIEpochBounds);
    prevCostPurPreAdapt2 = prevCostPurAll2(~doAdaptAllEpochBounds & includeISIEpochBounds);
    prevCostPurPreAdaptS = prevCostPurAllS(~doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostPurPreAdapt1 = postCostPurAll1(~doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostPurPreAdapt2 = postCostPurAll2(~doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostPurPreAdaptS = postCostPurAllS(~doAdaptAllEpochBounds & includeISIEpochBounds);
    
    prevCostMixPreAdapt1 = prevCostMixAll1(~doAdaptAllEpochBounds & includeISIEpochBounds);
    prevCostMixPreAdapt2 = prevCostMixAll2(~doAdaptAllEpochBounds & includeISIEpochBounds);
    prevCostMixPreAdaptS = prevCostMixAllS(~doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostMixPreAdapt1 = postCostMixAll1(~doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostMixPreAdapt2 = postCostMixAll2(~doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostMixPreAdaptS = postCostMixAllS(~doAdaptAllEpochBounds & includeISIEpochBounds);
    
    prevCostPurPostAdapt1 = prevCostPurAll1(doAdaptAllEpochBounds & includeISIEpochBounds);
    prevCostPurPostAdapt2 = prevCostPurAll2(doAdaptAllEpochBounds & includeISIEpochBounds);
    prevCostPurPostAdaptS = prevCostPurAllS(doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostPurPostAdapt1 = postCostPurAll1(doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostPurPostAdapt2 = postCostPurAll2(doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostPurPostAdaptS = postCostPurAllS(doAdaptAllEpochBounds & includeISIEpochBounds);
    
    prevCostMixPostAdapt1 = prevCostMixAll1(doAdaptAllEpochBounds & includeISIEpochBounds);
    prevCostMixPostAdapt2 = prevCostMixAll2(doAdaptAllEpochBounds & includeISIEpochBounds);
    prevCostMixPostAdaptS = prevCostMixAllS(doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostMixPostAdapt1 = postCostMixAll1(doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostMixPostAdapt2 = postCostMixAll2(doAdaptAllEpochBounds & includeISIEpochBounds);
    postCostMixPostAdaptS = postCostMixAllS(doAdaptAllEpochBounds & includeISIEpochBounds);
    
    
    allCostPurPreAdapt1 = allCostPurAll1(:, ~doAdapt & includedISIInds);
    allCostPurPreAdapt2 = allCostPurAll2(:, ~doAdapt & includedISIInds);
    
    allCostPurPostAdapt1 = allCostPurAll1(:, doAdapt & includedISIInds);
    allCostPurPostAdapt2 = allCostPurAll2(:, doAdapt & includedISIInds);
    
    allCostMixPreAdapt1 = allCostMixAll1(:, ~doAdapt & includedISIInds);
    allCostMixPreAdapt2 = allCostMixAll2(:, ~doAdapt & includedISIInds);
    
    allCostMixPostAdapt1 = allCostMixAll1(:, doAdapt & includedISIInds);
    allCostMixPostAdapt2 = allCostMixAll2(:, doAdapt & includedISIInds);
    
    deltaCostPurPostAdapt1 = reshape(postCostPurAll1(doAdaptAllEpochBounds & includeISIEpochBounds) - prevCostPurAll1(doAdaptAllEpochBounds & includeISIEpochBounds), numPostAdaptISI, []);
    deltaCostPurPostAdapt2 = reshape(postCostPurAll2(doAdaptAllEpochBounds & includeISIEpochBounds) - prevCostPurAll2(doAdaptAllEpochBounds & includeISIEpochBounds), numPostAdaptISI, []);
    
    deltaCostMixPostAdapt1 = reshape(postCostMixAll1(doAdaptAllEpochBounds & includeISIEpochBounds) - prevCostMixAll1(doAdaptAllEpochBounds & includeISIEpochBounds), numPostAdaptISI, []);
    deltaCostMixPostAdapt2 = reshape(postCostMixAll2(doAdaptAllEpochBounds & includeISIEpochBounds) - prevCostMixAll2(doAdaptAllEpochBounds & includeISIEpochBounds), numPostAdaptISI, []);
    
    clc;
    %     fprintf('Sets before recalibration algorithm: %d\nSets after recalibration algorithm: %d\n\n', numPreAdaptISI, numPostAdaptISI);
    fprintf('Sets after recalibration algorithm: %d\n\n', numPostAdaptISI);
    for epochBoundInd = 1:(minNumEpochs - 1)
        fprintf('Epoch Bound %d:\nPreRecal: %0.3f +/- %0.3f\n\nPostRecal (Pure): %0.3f +/- %0.3f\nPostRecal (Mixed): %0.3f +/- %0.3f\n\n\n', epochBoundInd,mean(deltaCostPurPreAdapt1(:,epochBoundInd)), std(deltaCostPurPreAdapt1(:,epochBoundInd)), nanmean(deltaCostPurPostAdapt1(:,epochBoundInd)), nanstd(deltaCostPurPostAdapt1(:,epochBoundInd)), nanmean(deltaCostMixPostAdapt1(:,epochBoundInd)), nanstd(deltaCostMixPostAdapt1(:,epochBoundInd)));
    end
    
    hCostPur3 = figure;
    aCost21 = subplot(2,1,1);
    aCost22 = subplot(2,1,2);
    % set(groot, 'defaultaxescolororder', lines(numPostAdapt));
    plot(aCost21, repelem(xAxisBounds,1, numPostAdaptISI), [prevCostPurPostAdapt1;postCostPurPostAdapt1], 'b.-', 'markersize', 15);
    title(aCost21, 'Post-Adaptive Boundary Differences E1 (Pure)');
    ylim(aCost21, [0 1]);
    plot(aCost22, repelem(xAxisBounds,1, numPostAdaptISI), [prevCostPurPostAdapt2;postCostPurPostAdapt2], 'b.-', 'markersize', 15);
    title(aCost22, 'Post-Adaptive Boundary Differences E2 (Pure)');
    ylim(aCost22, [0 1]);
    
    hCostPur4 = figure;
    % set(groot, 'defaultaxescolororder', lines(numPostAdapt));
    plot(repelem(xAxisBounds,1, numPostAdaptISI), [prevCostPurPostAdaptS;postCostPurPostAdaptS], 'b.-', 'markersize', 15);
    title('Post-Adaptive Boundary Differences Combined (Pure)');
    ylim([0 2]);
    
    hCostMix3 = figure;
    aCost21 = subplot(2,1,1);
    aCost22 = subplot(2,1,2);
    % set(groot, 'defaultaxescolororder', lines(numPostAdapt));
    plot(aCost21, repelem(xAxisBounds,1, numPostAdaptISI), [prevCostMixPostAdapt1;postCostMixPostAdapt1], 'b.-', 'markersize', 15);
    title(aCost21, 'Post-Adaptive Boundary Differences E1 (Mixed)');
    ylim(aCost21, [0 1]);
    plot(aCost22, repelem(xAxisBounds,1, numPostAdaptISI), [prevCostMixPostAdapt2;postCostMixPostAdapt2], 'b.-', 'markersize', 15);
    title(aCost22, 'Post-Adaptive Boundary Differences E2 (Mixed)');
    ylim(aCost22, [0 1]);
    
    hCostMix4 = figure;
    % set(groot, 'defaultaxescolororder', lines(numPostAdapt));
    plot(repelem(xAxisBounds,1, numPostAdaptISI), [prevCostMixPostAdaptS;postCostMixPostAdaptS], 'b.-', 'markersize', 15);
    title('Post-Adaptive Boundary Differences Combined (Mixed)');
    ylim([0 2]);
    
    hAllCostPur2 = figure;
    aAllCost21 = subplot(2,1,1);
    aAllCost22 = subplot(2,1,2);
    % set(groot, 'defaultaxescolororder', lines(numPostAdapt));
    plot(aAllCost21, repmat(xAxisAllCosts, 1, numPostAdaptISI), allCostPurPostAdapt1, '.-', 'markersize', 15);
    title(aAllCost21, 'Post-Adaptive Epoch Costs E1 (Pure)');
    ylim(aAllCost21, [0 1]);
    plot(aAllCost22, repmat(xAxisAllCosts, 1, numPostAdaptISI), allCostPurPostAdapt2, '.-', 'markersize', 15);
    title(aAllCost22, 'Post-Adaptive Epoch Costs E2 (Pure)');
    ylim(aAllCost22, [0 1]);
    
    hAllCostMix2 = figure;
    aAllCost21 = subplot(2,1,1);
    aAllCost22 = subplot(2,1,2);
    plot(aAllCost21, repmat(xAxisAllCosts, 1, numPostAdaptISI), allCostMixPostAdapt1, '.-', 'markersize', 15);
    title(aAllCost21, 'Post-Adaptive Epoch Costs E1 (Mixed)');
    ylim(aAllCost21, [0 1]);
    plot(aAllCost22, repmat(xAxisAllCosts, 1, numPostAdaptISI), allCostMixPostAdapt2, '.-', 'markersize', 15);
    title(aAllCost22, 'Post-Adaptive Epoch Costs E2 (Mixed)');
    ylim(aAllCost22, [0 1]);
    
end

%% Plot data
if doPerformancePlot
    % Plot the RF-difference based performance for post-adapt
    hHist = figure;
    a2 = axes;
    yyaxis left;
    bar(a2, tierDensitiesPostAdapt);
    text(1:length(tierDensitiesPostAdapt),tierDensitiesPostAdapt,num2str(tierDensitiesPostAdapt'),'vert','bottom','horiz','center');
    title(a2, sprintf('Histogram of Performance with Adaptive (N = %d)', numPostAdapt), 'fontsize', plotFontSize);
    xlim(a2, [0 (length(tierDensitiesPostAdapt) + 1)]);
    ylim([0 1.1*max(tierDensitiesPostAdapt)]);
    a2.XTickLabel = plotLabels;
    a2.XTickLabelRotation = 45;
    xlabel('Performance Tier', 'fontsize', plotFontSize);
    ylabel('Data Sets per Tier', 'fontsize', plotFontSize);
    
    hold(a2, 'on');
    yyaxis right;
    plot(xDataDifPost, averageDifCostPerTierPost, 'x', 'markersize', 15, 'linewidth', 2);
    plot(xDataSqePost, averageSqeCostPerTierPost, 'o', 'markersize', 15, 'linewidth', 2);
    ylabel('Cost Function Value', 'fontsize', plotFontSize);
    legend('Data (binned by performance)', 'Adjusted Difference Cost', 'Square Error Cost');
    set(findall(hHist,'-property','FontSize'),'FontSize',plotFontSize)
    
    
    if useDiff
        thisCostPost = difCostMeanPost;
        thisCostStr = 'Cost (f_D_i_f_f)';
    else
        thisCostPost = sqeCostMeanPost;
        thisCostStr = 'Cost (f_S_S_E)';
    end
    
    %% Apply limits to the data based on its clustering parameters and its performance
    limitStr = '';
    limitStrFileName = '';
    postInclude = true(size(thisCostPost));
    
    if useOfflineSort && useClustParamLimit
        limitStr = sprintf(', %s < %0.2f', clustParam, lRatioLimit);
        limitStrFileName = '_lim';
        postInclude = postInclude(:) & clustParamsInclude(:);
    end
    
    % Apply the limits
    thisCostPost = thisCostPost(postInclude);
    assignedPerformancePerSetPostGrouped = assignedPerformancePerSetPostGrouped(postInclude);
    
    % Sort the costs (needed for CDF), and save the sorting
    [costPostSorted, orderPost] = sort(thisCostPost);
    
    % Find the dates and times for each group of data that's better than
    % the clustering parameter limit, in the correct order.
    dateStrsInclude = dateStrsInclude(postInclude);
    timeStrsInclude = timeStrsInclude(postInclude);
    trialNumsInclude = trialNumsInclude(postInclude);
    dateStrsInclude = dateStrsInclude(orderPost);
    timeStrsInclude = timeStrsInclude(orderPost);
    trialNumsInclude = trialNumsInclude(orderPost);
    
    % Find the data indices of each experiment, sorted according to cost
    allDataInds = allDataInds(postInclude);
    allDataInds = allDataInds(orderPost);
    
    % Sort the assignedPerformance vectors (indicating how well each Set
    % performed) according to sorting quality
    assignedPerformancePerSetPostGroupedSorted = assignedPerformancePerSetPostGrouped(orderPost);
    
    % Create the x-axis of the CDF
    numCostPost = length(thisCostPost);
    postX = (1:numCostPost)/numCostPost;
    
    % Determine which Sets have "acceptable" performance
    minPerformance = 1 - ~doAcceptableControlabilityLabeling;
    accIndsPost = assignedPerformancePerSetPostGroupedSorted > minPerformance; % Inds in costPostSorted that are acceptable (achieve the minimum requirements for control)
    notAccIndsPost = assignedPerformancePerSetPostGroupedSorted == minPerformance; % Inds in costPostSorted that are not acceptable (do not achieve the minimum requirements for control)
    
    % Split the x-axes into acceptable and not acceptable (needed to plot
    % both on the same plot)
    postXAcc = postX(accIndsPost);
    postXNotAcc = postX(notAccIndsPost);
    
    % Split costs according to acceptable and not acceptable
    thisCostPostAcc = costPostSorted(accIndsPost);
    thisCostPostNotAcc = costPostSorted(notAccIndsPost);
    
    numCostPostAcc = length(thisCostPostAcc);
    numCostPostNotAcc = length(thisCostPostNotAcc);
    
    % Choose sample data sets to show the POT of on the plot (low, medium,
    % and high cost)
    allQuartiles = linspace(0, 1, numExamples);
    allExampleInds = nan(numExamples, 1);
    
    if doAcceptableControlabilityLabeling && exampleQuartilesFromControllable
        possibleExampleInds = find(accIndsPost);
        examplePoolCosts = thisCostPostAcc;
    else
        possibleExampleInds = 1:numCostPost;
        examplePoolCosts = costPostSorted;
    end
    
    for exampleNum = 1:numExamples
        [~, examplePoolInd] = min(abs(examplePoolCosts - quantile(examplePoolCosts, allQuartiles(exampleNum))));
        allExampleInds(exampleNum) = possibleExampleInds(examplePoolInd);
    end
    
    %% Plot the probability space plot
    if doProbSpacePlot
        meanTPFrac1 = nanmean(tPFrac1, 1);
        meanFAFrac1 = nanmean(fAFrac1, 1);
        meanTPFrac2 = nanmean(tPFrac2, 1);
        meanFAFrac2 = nanmean(fAFrac2, 1);
        
        hProb = figure;
        a = axes;
        hold(a, 'on');
        
        h = plot(a, meanTPFrac1, meanFAFrac2, 'b.', 'markersize', plotMarkerSize);
        j = plot(a, meanTPFrac2, meanFAFrac1, 'r.', 'markersize', plotMarkerSize);
        
        hAll = [h(1) j(1)];
        lineStr = {'S_A', 'S_B'};
        
        if doConnectPoints
            plot(a, [meanTPFrac1;meanTPFrac2], [meanFAFrac2;meanFAFrac1], 'k:');
        end
        
        if doVectorField && ~doConnectPoints
            costSensTP = @(TP, FA) -1*(-3 + 2*TP); % Each is multiplied by -1 to show the direction of DECREASING cost, instead of INCREASING
            costSensFA = @(TP, FA) -1*(ones(size(FA)));
            [pointX, pointY] = meshgrid(linspace(0, 1, vectorFieldMeshSize));
            vecX = costSensTP(pointX, pointY);
            vecY = costSensFA(pointX, pointY);
            
            k = quiver(a, pointX, pointY, vecX, vecY, vectorScale, 'c');
            hAll = [hAll k(1)];
            lineStr = [lineStr 'Cost Sensitivity'];
        end
        
        ylim(a, [0 1]);
        xlim(a, [0 1]);
        xlabel(a, 'True Positives');
        ylabel(a, 'False Alarms');
        title(a, 'Probability Space summary of results');
        legend(a, hAll, lineStr);
        set(findall(hProb,'-property','FontSize'),'FontSize',plotFontSize)
    end
    
    %% Plot the string plot
    if doStringPlot
        if useDiff
            f = difCostFn;
        else
            f = sqeCostFn;
        end
        allIncludedCostsAOnline = nanmean(f(tPOnlineFrac1, fAOnlineFrac2), 1)';
        allIncludedCostsBOnline = nanmean(f(tPOnlineFrac2, fAOnlineFrac1), 1)';
        
        allIncludedCostsAOffline = nanmean(f(tPOfflineFrac1, fAOfflineFrac2), 1)';
        allIncludedCostsBOffline = nanmean(f(tPOfflineFrac2, fAOfflineFrac1), 1)';
        
        allIncludedCostsTotalOnline = nanmean([allIncludedCostsAOnline, allIncludedCostsBOnline], 2);
        allIncludedCostsTotalOffline = nanmean([allIncludedCostsAOffline, allIncludedCostsBOffline], 2);
        
        xSingle = [0;1];
        x = repmat(xSingle, 1, size(allIncludedCostsTotalOnline,1));
        
        hS = figure;
        a = axes;
        hold(a, 'on');
        h = plot(a, x, [allIncludedCostsTotalOnline allIncludedCostsTotalOffline]', 'b.-', 'markersize', plotMarkerSize);
        j = plot(a, xSingle, median([allIncludedCostsTotalOnline allIncludedCostsTotalOffline], 1), 'rx-', 'linewidth', 2, 'markersize', plotMarkerSize);
        set(a, 'xtick', [0 1], 'xticklabel', {'Online', 'Offline'});
        xlim([-1 2]);
        ylim([0 1]);
        ylabel(a, 'Cost');
        title(a, 'Comparison of Online to Offline sort costs');
        legend([h(1) j], 'All Data Change', 'Median Data Change');
        set(findall(hS,'-property','FontSize'),'FontSize',plotFontSize)
    end
    
    %% Plot the CDF
    maxCost = max(thisCostPost(:))*1.1;
    h1 = figure;
    a1 = axes;
    allH = plot(a1, thisCostPostAcc, postXAcc, 'kx', thisCostPostNotAcc, postXNotAcc, 'ko', 'linewidth', plotLineWidth, 'markersize', plotMarkerSize);
    if doAcceptableControlabilityLabeling
        lineStrs = {sprintf('Adaptive Algorithm Data Set , controlled (n = %d)', numCostPostAcc), sprintf('Adaptive Algorithm Data Set, not controlled (n = %d)', numCostPostNotAcc)};
    else
        lineStrs = {sprintf('Adaptive Algorithm Dataset (n = %d)', numCostPost)};
    end
    hold(a1, 'on');
    plot(a1, costPostSorted, postX, 'k:', 'linewidth', plotLineWidth);
    if doQuartilesPlot
        j = plot(a1, costPostSorted(allExampleInds), postX(allExampleInds), 'ro', 'markersize', plotMarkerSize + 5, 'linewidth', plotLineWidth);
        allH = [allH;j(1)];
        lineStrs = [lineStrs 'Example datasets'];
    end
    ylabel('Fraction of All Sets', 'fontsize', plotFontSize);
    xlabel(thisCostStr, 'fontsize', plotFontSize);
    title(['CDF of Adaptive Algorithm Dataset' isiStr limitStr], 'fontsize', plotFontSize)
    legend(allH, lineStrs, 'location', 'southeast');
    set(findall(h1,'-property','FontSize'),'FontSize',plotFontSize)
    xlim([0 maxCost]);
    
    if doQuartilesPlot
        % Make detail plots
        hPost = nan(numExamples, 1);
        sizeDif = .15;
        %         subAxisSizeDif = .1;
        for exampleNum = 1:numExamples
            % Create figures and axes
            if doWaveformPlots
                hPost(exampleNum) = figure;
                aPost = subplot(2,1,1);
                aPostWav1 = subplot(2,2,3);
                aPostWav2 = subplot(2,2,4);
                
                posMain = get(aPre, 'position');
                posMain([2 4]) = posMain([2 4]) + [-sizeDif, sizeDif];
                set(aPre,'position', posMain);
                set(aPost,'position', posMain);
                
                posSub1 = get(aPreWav1, 'position');
                posSub2 = get(aPreWav2, 'position');
                posSub1(4) = posSub1(4) - sizeDif;
                posSub2(4) = posSub2(4) - sizeDif;
                set(aPreWav1, 'position', posSub1);
                set(aPreWav2, 'position', posSub2);
                set(aPostWav1, 'position', posSub1);
                set(aPostWav2, 'position', posSub2);
            else
                hPost(exampleNum) = figure;
                aPost = axes;
            end
            
            hold(aPost, 'on');
            
            % Get the indices in stimTimesAll/performance vectors (tPFrac1,
            % etc.)
            iPostIncluded = allExampleInds(exampleNum);
            iPost = allDataInds(iPostIncluded);
            
            % Ready the subsampling, used to make relatively few "X" and "O"
            % markers so as to not clutter the screen
            numMarkers = 50;
            subSampInds = round(linspace(1,size(stimTimesAll, 2), numMarkers));
            
            % Plot the lines
            plot(aPost, stimTimesAll(iPost, :), tPFrac1(:,iPost), 'b-', stimTimesAll(iPost, :), fAFrac1(:,iPost), 'b--', stimTimesAll(iPost, :), tPFrac2(:,iPost), 'r-', stimTimesAll(iPost, :), fAFrac2(:,iPost), 'r--', 'linewidth', plotLineWidth);
            
            % Plot the markers
            jPost = plot(aPost, stimTimesAll(iPost, subSampInds), tPFrac1(subSampInds,iPost), 'bx', stimTimesAll(iPost, subSampInds), fAFrac1(subSampInds,iPost), 'bo', stimTimesAll(iPost, subSampInds), tPFrac2(subSampInds,iPost), 'ro', stimTimesAll(iPost, subSampInds), fAFrac2(subSampInds,iPost), 'rx', 'linewidth', plotLineWidth, 'markersize', plotMarkerSize);
            
            xlabel(aPost, 'Time (s)');
            ylabel(aPost, 'Response Fraction');
            title(aPost, sprintf('POT of example #%d, (%s %s %d)', exampleNum, dateStrsInclude{iPostIncluded}, timeStrsInclude{iPostIncluded}, trialNumsInclude(iPostIncluded)));
            legend(aPost, jPost, 'P(A|S_A)', 'P(A|S_B)', 'P(B|S_B)', 'P(B|S_A)');
            ylim(aPost, [0 1]);
            
            % Plot the waveforms
            if doWaveformPlots
                % Find which traces are actually spikes
                spikeIndsPost1 = find(spikeSR1(iPost, :));
                spikeIndsPost2 = find(spikeSR2(iPost, :));
                
                numTraces1Post = length(spikeIndsPost1); % Get the number of traces
                numTraces2Post = length(spikeIndsPost2); % Get the number of traces
                randInds1Post = spikeIndsPost1(randperm(numTraces1Post, min(numTracesMax, numTraces1Post)));
                randInds2Post = spikeIndsPost2(randperm(numTraces2Post, min(numTracesMax, numTraces2Post)));
                
                % Get randomly ordered traces
                
                traces1Post = convToMicroV*cell2mat(vT1{iPost}(randInds1Post))';
                traces2Post = convToMicroV*cell2mat(vT2{iPost}(randInds2Post))';
                
                % Extract 32-sample waveforms to plot
                waveformXAxis = ((-numSamplesPrePeak):(numSamples - numSamplesPrePeak - 1));
                
                waveFormInds1Post = spikeHitInds{iPost}(1, randInds1Post) + waveformXAxis';
                wFInds1PostAdjustmentMask = waveFormInds1Post < 1 | size(traces1Post, 1) < waveFormInds1Post; % All inds that are within scope of the trace (not 0 or below, nor above the size of the trace)
                waveFormInds1Post(wFInds1PostAdjustmentMask) = 1; % Make sure that the traces can be properly indexed (is junk data, will be later masked out with nans)
                waveForms1Post = nan(size(waveFormInds1Post));
                waveForms1Post(:) = traces1Post(sub2ind(size(traces1Post), waveFormInds1Post, repmat(1:size(traces1Post,2), numSamples, 1))); % Get all waveforms
                waveForms1Post(wFInds1PostAdjustmentMask) = nan; % Set to nan any samples that are out of scope
                
                waveFormInds2Post = spikeHitInds{iPost}(2, randInds2Post) + waveformXAxis';
                wFInds2PostAdjustmentMask = waveFormInds2Post < 1 | size(traces2Post, 1) < waveFormInds2Post; % All inds that are within scope of the trace (not 0 or below, nor above the size of the trace)
                waveFormInds2Post(wFInds2PostAdjustmentMask) = 1; % Make sure that the traces can be properly indexed (is junk data, will be later masked out with nans)
                waveForms2Post = nan(size(waveFormInds2Post));
                waveForms2Post(:) = traces2Post(sub2ind(size(traces2Post), waveFormInds2Post, repmat(1:size(traces2Post,2), numSamples, 1))); % Get all waveforms
                waveForms2Post(wFInds2PostAdjustmentMask) = nan; % Set to nan any samples that are out of scope
                
                % Plot the waveforms
                plot(aPostWav1, waveformXAxis, waveForms1Post, 'b:');
                plot(aPostWav2, waveformXAxis, waveForms2Post, 'r:');
                
                xlabel(aPostWav1, 'Sample');
                ylabel(aPostWav1, 'Potential (uV)');
                title(aPostWav1, sprintf('Unit A, n = %d', size(waveForms1Post,2)));
                
                xlabel(aPostWav2, 'Sample');
                title(aPostWav2, sprintf('Unit B, n = %d', size(waveForms2Post,2)));
            end
            set(findall(hPost(exampleNum),'-property','FontSize'),'FontSize',plotFontSize);
        end
    end
    
    if doQualityVsPerformancePlot && useOfflineSort
        % Plot the clustering quality metric vs the performance of each
        % pair for ALL recorded pairs.
        % Using max (worst-case) clustering metric, and average performance
        % between the two stimulations
        
        % Process the clustering parameter to find the worst-case
        clustParamWorst = max(allClustParams, [], 3);
        
        % Extract the selected cost
        if useDiff
            thisCost = difCostMeanTotal;
        else
            thisCost = sqeCostMeanTotal;
        end
        
        % Get columnized parameters
        thisCost = thisCost(:);
        clustParamWorst = clustParamWorst(:);
        
        % Plot the parameters
        hQP = figure;
        a = axes;
        plot(a, clustParamWorst, thisCost, '.', 'markersize', plotMarkerSize);
        xlabel(a, 'L_r_a_t_i_o');
        ylabel(a, thisCostStr);
        title(a, 'Performance as a function of cluster quality');
        curXLim = xlim;
        
        xlim([0 curXLim(2)]);
        ylim([0 1]);
    end
end


if doSummarySpreadsheet
    try
        delete(xlsFileName);
        xlswrite(xlsFileName, finalTable);
    catch e
        fprintf('Could not write to %s.  Failed with error:\n\n*****\n\n%s\n\n*****\n\n', xlsFileName, getReport(e));
    end
end

%% Save Data
if doSaveData && doStatisticalTests
   save(dataSaveFile, 'dateStrs', 'timeStrs', 'CQATrue', 'CQBTrue', 'CQTrue', 'sigmaA', 'sigmaB', 'steA', 'steB', 'cqEstMean', 'cqEstSte', 'thisHdrBounds', 'conclusiveTotalControl', 'goodTotalControl', 'goodConclusiveTotalControl', 'conclusiveAControl', 'conclusiveBControl', 'goodAControl', 'goodBControl', 'goodConclusiveAControl', 'goodConclusiveBControl', 'goodConclusiveOneWayControl', 'cqMeansSorted', 'cqMeanSortInds', 'cqMean', 'cqStd', 'sampleZScore', 'stimDifferentiabilityPVal');
end

%% Save plots
if doSaveFigs
    if ~exist(parentFigureDir, 'dir')
        mkdir(parentFigureDir);
    end
    
    if useOfflineSort
        sortDirName = 'Offline';
    else
        sortDirName = 'Online';
    end
    
    mainFigureDir = sprintf('%s%s%s', parentFigureDir, sortDirName, filesep);
    
    if ~exist(mainFigureDir, 'dir')
        mkdir(mainFigureDir);
    end
    
    if doPerformancePlot
        %
        fName2 = sprintf('%sPerformanceHistogram.', mainFigureDir);
        saveFigure(hHist, [fName2 'png']);
        savefig(hHist, [fName2 'fig']);
        
        fName1 = sprintf('%sDataSetCostCDF%s.', mainFigureDir, limitStrFileName);
        saveFigure(h1, [fName1 'png']);
        savefig(h1, [fName1 'fig']);
        
        if doProbSpacePlot
            fName1 = sprintf('%sProbabilitySpace%s.', mainFigureDir, limitStrFileName);
            saveFigure(hProb, [fName1 'png']);
            savefig(hProb, [fName1 'fig']);
        end
        
        if doStringPlot
            fName1 = sprintf('%sClusteringComparisonString%s.', mainFigureDir, limitStrFileName);
            saveFigure(hS, [fName1 'png']);
            savefig(hS, [fName1 'fig']);
        end
        
        if doQuartilesPlot
            for exampleNum = 1:numExamples
                
                fName2 = sprintf('%sQuartile%dExPost%s.', mainFigureDir, exampleNum - 1, limitStrFileName);
                saveFigure(hPost(exampleNum), [fName2 'png']);
                savefig(hPost(exampleNum), [fName2 'fig']);
            end
        end
        
        if doQualityVsPerformancePlot && useOfflineSort
            fName1 = sprintf('%sQualityVsPerformancePlot.', mainFigureDir);
            saveFigure(hQP, [fName1 'png']);
            savefig(hQP, [fName1 'fig']);
        end
        
    end
    
    if doCostFunctionPlot
        fName3 = sprintf('%sBoundDifPurePostInd.', mainFigureDir);
        saveFigure(hCostPur3, [fName3 'png']);
        savefig(hCostPur3, [fName3 'fig']);
        
        fName4 = sprintf('%sBoundDifPurePostCom.', mainFigureDir);
        saveFigure(hCostPur4, [fName4 'png']);
        savefig(hCostPur4, [fName4 'fig']);
        
        
        fName3 = sprintf('%sBoundDifMixedPostInd.', mainFigureDir);
        saveFigure(hCostMix3, [fName3 'png']);
        savefig(hCostMix3, [fName3 'fig']);
        
        fName4 = sprintf('%sBoundDifMixedPostCom.', mainFigureDir);
        saveFigure(hCostMix4, [fName4 'png']);
        savefig(hCostMix4, [fName4 'fig']);
        
        fName2 = sprintf('%sEpochPurePost.', mainFigureDir);
        saveFigure(hAllCostPur2, [fName2 'png']);
        savefig(hAllCostPur2, [fName2 'fig']);
        
        fName4 = sprintf('%sEpochMixedPost.', mainFigureDir);
        saveFigure(hAllCostMix2, [fName4 'png']);
        savefig(hAllCostMix2, [fName4 'fig']);
    end
    
    if doStatisticalTestingPlot
        fName1 = sprintf('%sControlQualitySummary.', mainFigureDir);
        saveFigure(hStatsControl, [fName1 'png']);
        savefig(hStatsControl, [fName1 'fig']);
        
        fName1 = sprintf('%sStimDifferentiability.', mainFigureDir);
        saveFigure(hStatsDiff, [fName1 'png']);
        savefig(hStatsDiff, [fName1 'fig']);
    end
end

function y = trimEpochStartIndHelperFun(x, minNumEpochs, numTotalSequences)
x = [x numTotalSequences];
y = x(1:minNumEpochs + 1)';
end