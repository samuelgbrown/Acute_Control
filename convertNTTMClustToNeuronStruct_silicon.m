function convertNTTMClustToNeuronStruct_silicon(varargin)
% This function will convert the output clusters from MClust (that were
% input from an .nse file created by convertTankToNSE) into a NEW
% neuronStruct file, whose resultsData reflects the new sorting (to analyze
% the units' performance using offline sorting methods).
%
% Input:
% inputFileName - the date and time name from an .nse file, without the
% electrode deisgnation (e.g. "20170803_1721").  Both .t files, made from
% the A electrode and the B electrode, will be required.
%
% BEFORE RUNNING:
% Create an .nse file from the convertTankToNSE script, keeping the
% resulting expTankNums and expTimeStrs variables intact (if they are
% destroyed, running convertTankToNSE once more will replace them while
% preserving any .t files already made).

% Perform MClust clustering on BOTH units in the original neuronStruct
% file.  Run MClust, loading one .nse file created by convertTankToNSE and
% sort using the desired clustering features (look up MClust documentation
% for more information on clustering).  When this is done, export the
% cluster and save the .t file.  Once one .nse file has been sorted and its
% .t file saved, sort the other.

%% Set parameters
dataFolder = 'D:\samuelgb\Documents\MATLAB\Acute Control\Data\Second Round Data\'; % Where the data should be
numStimsPerExp = 2000; % Number of stims per experiment
timespanLim = 20; % Acquire all threshold crossings between the onset of stimulation and this value (in ms) after
maxTimeRecorded = 50; % Number of ms of data recorded by Matlab after each stim
sampPerRecording = 1000; % The number of samples per recorded waveform in the Matlab program (covers maxTimeRecording ms of data)
indsAfterMClustTimeStart = 0;
indsAfterMClustTimeEnd = 200; % The number of indices after the indicated time in the MClust generated t-file that the script will look for a peak
indsAfterMClustTimeNum = indsAfterMClustTimeEnd - indsAfterMClustTimeStart + 1; % Size of waveform after the t-file timestamps to check


%% Make sure that this file can be created
if evalin('base', 'exist(''expTankNums'', ''var'') && exist(''expTimeStrs'', ''var'') && exist(''expTimeBounds'', ''var'') && exist(''expPenetrationNums'', ''var'')&& exist(''timespanLim'', ''var'')')
    expTankNums = evalin('base', 'expTankNums');
    expTimeStrs = evalin('base', 'expTimeStrs');
    expTimeBounds = evalin('base', 'expTimeBounds');
    expPenetrationNums = evalin('base', 'expPenetrationNums');
    timespanLim = evalin('base', 'timespanLim');
else
    error('Cannot create neuronStruct file without ''expTankNums'', ''expTimeStrs'', ''expPenetrationNums'', and ''expTimeBounds''.  Please run convertTankToNTT.');
end

%% Extract inputs
filePrefix = 'neuronStruct';
inputFileName = [];
if nargin > 1
    error('Too many inputs')
elseif nargin == 1
    inputDateTimeName = varargin{1};
    inputFileName = [filePrefix inputDateTimeName '.mat'];
    
    if ~exist(inputFileName, 'file')
        fprintf('Could not find file %s.\n', inputFileName);
        return;
    end
end

%% Load the original neuronStruct file
% User must select the original neuronStruct file (to serve as a
% "prototype")
haveFileOriginal = false;
while ~haveFileOriginal
    if ~isempty(inputFileName)
        [origFilePath, origFileName] = fileparts(which(inputFileName)); % Get the full path and name from the input file
        inputFileName = [];
    else
        hW = warndlg('Please select the originally recorded neuronStruct file.', 'File selection');
        uiwait(hW);
        [origFileName, origFilePath] = uigetfile('*.mat', '', dataFolder); % Get the file from the user via a getfile UI
    end
    
    if ~ischar(origFileName)
        % See why the user canceled getting the file
        answer = questdlg('Skip this file?', 'Canceled file load', 'Try again', 'Skip', 'Skip');
        if strcmp(answer, 'Skip')
            return;
        else
            continue;
        end
    end
    
    fileNameParts1 = strsplit(origFileName, '.');
    origFileName = fileNameParts1{1}; % Remove the file suffix
    fileNameParts2 = strsplit(origFileName, '_');
    thisTimeStr = fileNameParts2{2};
    isInExpTimeStrs = cellfun(@(a) strcmp(a, thisTimeStr), expTimeStrs);
    if ~any(isInExpTimeStrs)
        hW = warndlg('Please select different file, this one was not a part of the most recently loaded .nse file.', 'File selection');
        uiwait(hW);
        continue;
    end
    
    % Get information about this file and its location in the TDT data
    tankNum = expTankNums(isInExpTimeStrs); % Tank number
    thisExpTimeBounds = expTimeBounds(isInExpTimeStrs, :); % Experiment time boundaries
    thisExpPenetrationNums = expPenetrationNums(isInExpTimeStrs, :); % Experiment penetration numbers
    
    % Construct the full file name for mat file, and load it
    fullFileName = [origFilePath filesep origFileName];
    load(fullFileName);
    
    % If the file has been loaded, then continue
    if exist('neuronStruct', 'var')
        haveFileOriginal = true;
    end
end

%% Load the two t-files
% Guess the name of the first unit's t-file
% guessName = [origPathName filesep filePrefix '_01.t'];

% Navigate to the directory containing the original file, and ask the user
% to pick out the .t file corresponding to the first unit
cd(origFilePath);
haveFileFirst = false;
haveFileSecond = false;

dateTimeStr = regexp(origFileName, '\d{8}_\d{4}', 'match');
dateTimeStr = dateTimeStr{1};

possibleT1FileName = sprintf('%s_A_01.t', dateTimeStr);
possibleT2FileName = sprintf('%s_B_01.t', dateTimeStr);

if exist(possibleT1FileName, 'file') && exist(possibleT2FileName, 'file')
    haveFileFirst = true;
    haveFileSecond = true;
    
    t1FullFileName = possibleT1FileName;
    t2FullFileName = possibleT2FileName;
end

while ~haveFileFirst
    hW = warndlg(sprintf('Please select the MClust .t file corresponding to the first electrode (A) for file %s.', origFileName), 'File selection');
    uiwait(hW);
    [t1FileName, t1PathName] = uigetfile('*.t', '', origFilePath);
    
    if ~ischar(t1FileName)
        % See why the user canceled getting the file
        answer = questdlg('Skip this file?', 'Canceled file load', 'Try again', 'Skip', 'Skip');
        if strcmp(answer, 'Skip')
            return;
        end
    else
        haveFileFirst = true;
        t1FullFileName = [t1PathName t1FileName];
    end
end


% Ask the user to pick out the .t file corresponding to the second unit
% if strcmp(t1FullFileName, guessName)
%     guessName = [origPathName filesep filePrefix '_02.t'];
% else
%     guessName = [origPathName filesep filePrefix '_01.t'];
% end
while ~haveFileSecond
    hW = warndlg(sprintf('Please select the MClust .t file corresponding to the second electrode (B) for file %s.', origFileName), 'File selection');
    uiwait(hW);
    [t2FileName, t2PathName] = uigetfile('*.t', '', origFilePath);
    
    if ~ischar(t2FileName)
        % See why the user canceled getting the file
        answer = questdlg('Skip this file?', 'Canceled file load', 'Try again', 'Skip', 'Skip');
        if strcmp(answer, 'Skip')
            return;
        end
    else
        haveFileSecond = true;
        t2FullFileName = [t2PathName t2FileName];
    end
end


%% Load the two t files
fprintf('Using files for MClust .t files:\nt1: %s\nt2: %s\n\n', t1FullFileName, t2FullFileName);
fprintf('Reading t files...\n');
tData{1} = t_load(t1FullFileName);
tData{2} = t_load(t2FullFileName);

%% Compare these to the TDT file
% Get the timeString
parts = strsplit(dateTimeStr, '_');
dateStr = parts{1};
timeStr = parts{2};

% Determine which tank number this experiment's data is in
tankSubName = sprintf('Acute_%s_%0.2d', dateStr, tankNum);
tankFileName = sprintf('%s%s%s%sTutTank_%s.tsq', origFilePath, filesep, tankSubName, filesep, tankSubName);
if ~exist(tankFileName, 'file')
    fprintf('%s does not exist, cannot load tank to convert to neuron struct.\n', tankFileName);
    return;
end

% Load tank
block = tdt_block_sam(tankFileName); % Create a block object
% fullProbeA = getdata(block, 'Raw1', 'channel', 1, 'interval', expTimeBounds); % Get the filtered data for probe 1
% fullProbeB = getdata(block, 'Raw1', 'channel', 2, 'interval', expTimeBounds); % Get the filtered data for probe 1
pNum = getdata(block, 'TDat', 'channel', 4, 'interval', thisExpTimeBounds); % Get the filtered data for probe 1
% expectedResults = getdata(block, 'TDat', 'channel', 5, 'interval', thisExpTimeBounds); % Get the expected results of each stimulation

numTDTPens = length(thisExpPenetrationNums);

% Find which trialData rows in neuronStruct correspond to each stim in the
% TDT data
NSPenInds = nan(numTDTPens,1);
NSPenNums = nan(numTDTPens,1);
trialDataIsNotEmpty = false(1,length(neuronStruct{1}.trialData));
NSPenNumLoc = 1;
trialStart = 0; % Mechanism to deal with bug in previous code that left many empty entries in trialData
for i = 1:length(neuronStruct{1}.trialData)
    if isempty(neuronStruct{1}.trialData(i).TDTPenetNum)
        trialStart = trialStart + 1;
       continue;
    end
    trialDataIsNotEmpty(i) = true; % This element is not empty
    
    [isInTDTPenNums, locInTDTPenNums] = ismember(neuronStruct{1}.trialData(i).TDTPenetNum, thisExpPenetrationNums);
    
    if any(isInTDTPenNums)
        NSPenNums(NSPenNumLoc) = thisExpPenetrationNums(locInTDTPenNums);
        NSPenInds(NSPenNumLoc) = i - trialStart;
        NSPenNumLoc = NSPenNumLoc + 1;
    end
end

% Adjust for empty entries in trialData (ex. 20170323_1536)
neuronStruct{1}.trialData = neuronStruct{1}.trialData(trialDataIsNotEmpty);
neuronStruct{2}.trialData = neuronStruct{2}.trialData(trialDataIsNotEmpty);

% NSPenNums(1:numTDTPens) = NSPenNums(trialStart:(trialStart + numTDTPens - 1));
% NSPenNums((numTDTPens + 1):end) = [];
% neuronStruct{1}.trialData((numTDTPens + 1):end) = [];
% neuronStruct{2}.trialData((numTDTPens + 1):end) = [];

sizeData = length(neuronStruct{1}.voltageTraces{1});
timeGranularity = maxTimeRecorded/sampPerRecording; % The number of seconds per sample

% Create a blank sortCodes cell array to populate neuronStruct
sortCodesA = zeros(neuronStruct{1}.resultsDataAllocated, sizeData);
sortCodesB = zeros(neuronStruct{2}.resultsDataAllocated, sizeData);

for trialDataInd = 1:length(NSPenInds)
    thisTDTPen = NSPenNums(trialDataInd);
    NSPenInd = NSPenInds(trialDataInd);
    if isnan(thisTDTPen)
       continue;
    end
    
    thisPenInds = find(pNum.vals == thisTDTPen, numStimsPerExp, 'last');
    
    % Determine which stimulation (in TDT) each spike belongs to (if any)
    %     timeBoundsA = neuronStruct{1}.triggerTimespan(neuronStruct{1}.resultsDataCurLoc - 1, :)/1000;
    %     timeBoundsB = neuronStruct{2}.triggerTimespan(neuronStruct{2}.resultsDataCurLoc - 1, :)/1000;
    timeBoundsA = [0 maxTimeRecorded]/1000;
    timeBoundsB = [0 maxTimeRecorded]/1000;
    stimTimes = pNum.times(thisPenInds); % Times of each stimulus
    
    binEdgesMatA = bsxfun(@plus, timeBoundsA,stimTimes)'; % A matrix containing the start and end times of each stim response period, in TDT time
    binEdgesA = sort(binEdgesMatA(:)); % The bin edges of A.  Every other pair of numbers (inds 1 and 2, 3 and 4, 5 and 6, etc.) indicate the beginning and end of the stim response period, in which a spike is considered a stim response.
    
    binEdgesMatB = bsxfun(@plus, timeBoundsB,stimTimes)'; % A matrix containing the start and end times of each stim response period, in TDT time
    binEdgesB = sort(binEdgesMatB(:)); % The bin edges of A.  Every other pair of numbers (inds 1 and 2, 3 and 4, 5 and 6, etc.) indicate the beginning and end of the stim response period, in which a spike is considered a stim response.
    
    [~, binNumsA] = histc(tData{1}, binEdgesA); % Find where each stimulus is
    [~, binNumsB] = histc(tData{2}, binEdgesB); % Find where each stimulus is
    
    
%     binNumsA(binNumsA == 0) = max(binNumsA);
%     binNumsB(binNumsB == 0) = max(binNumsB);
    
    tFileStimTimesA = tData{1};
    tFileStimTimesB = tData{2};
    
    % Remove all spikes that fall into an even numbered bin (those that fell BETWEEN stim response periods, instead of INSIDE of them)
    removeIndsA = mod(binNumsA,2) == 0;
    binNumsA(removeIndsA) = [];
    tFileStimTimesA(removeIndsA) = [];
    
    removeIndsB = mod(binNumsB,2) == 0;
    binNumsB(removeIndsB) = [];
    tFileStimTimesB(removeIndsB) = [];
    
%     % Remove all spikes that are doubles of spikes already in a stim's
%     % response region
%     [binNumsA, iA] = unique(binNumsA);
%     tFileStimTimesA = tFileStimTimesA(iA);
%     
%     [binNumsB, iB] = unique(binNumsB);
%     tFileStimTimesB = tFileStimTimesB(iB);
    
    % Get the stim numbers within this experiment that each spike belongs to
    stimNumsA = (binNumsA + 1)/2; % Convert all odd indices (1, 3, 5...) to counting numbers (1, 2, 3...)
    stimNumsB = (binNumsB + 1)/2; % Convert all odd indices (1, 3, 5...) to counting numbers (1, 2, 3...)
    
    firstStimTime = neuronStruct{1}.trialData(NSPenInd).stimTimes(1); % Time in resultsData of the first stim
    firstStimIndA = find(neuronStruct{1}.resultsData(:,4) == firstStimTime);
    firstStimIndB = find(neuronStruct{2}.resultsData(:,4) == firstStimTime); % Should be idential to firstStimIndA
    
    resultsDataIndsA = (firstStimIndA:(firstStimIndA + numStimsPerExp - 1))'; % Indices in resultsData of ALL stims for this experiment
    resultsDataIndsB = (firstStimIndB:(firstStimIndB + numStimsPerExp - 1))';
    
    resultsDataIndsSpikesA = stimNumsA + firstStimIndA - 1; % Indices in resultsData of all stims for this experiment that threshold crossed on A
    resultsDataIndsSpikesB = stimNumsB + firstStimIndB - 1; % Indices in resultsData of all stims for this experiment that threshold crossed on B
    
    % Retreive all of the required traces
    tracesA = cell2mat(cellfun(@double, neuronStruct{1}.voltageTraces(resultsDataIndsSpikesA), 'uniformoutput', false))';
    tracesB = cell2mat(cellfun(@double, neuronStruct{2}.voltageTraces(resultsDataIndsSpikesB), 'uniformoutput', false))';
    
    % Find the threshold crossing times of each spike
    indsPerMS = sizeData/maxTimeRecorded;
    crossIndsStartA = round((tFileStimTimesA - binEdgesA(binNumsA))*1000*indsPerMS); % The start of the small range in which the maximum of a spike should be found
    crossIndsStartB = round((tFileStimTimesB - binEdgesB(binNumsB))*1000*indsPerMS);
    crossIndsStartA((crossIndsStartA + indsAfterMClustTimeEnd) > size(tracesA,1)) = size(tracesA,1) - indsAfterMClustTimeEnd; % Ensure that indices are bounded by the length of the trace
    crossIndsStartB((crossIndsStartB + indsAfterMClustTimeEnd) > size(tracesB,1)) = size(tracesB,1) - indsAfterMClustTimeEnd;
    crossSubsRangeMatA = (crossIndsStartA + (indsAfterMClustTimeStart:indsAfterMClustTimeEnd))';
    crossSubsRangeMatB = (crossIndsStartB + (indsAfterMClustTimeStart:indsAfterMClustTimeEnd))';
    crossSubsRangeMatA(crossSubsRangeMatA == 0) = 1;
    crossSubsRangeMatB(crossSubsRangeMatB == 0) = 1;
    crossIndsRangeMatA = sub2ind(size(tracesA),crossSubsRangeMatA, repmat(1:size(tracesA,2), indsAfterMClustTimeNum, 1));
    crossIndsRangeMatB = sub2ind(size(tracesB),crossSubsRangeMatB, repmat(1:size(tracesB,2), indsAfterMClustTimeNum, 1));
    [~, crossIndsWithinDataPartA] = max(tracesA(crossIndsRangeMatA), [], 1);
    [~, crossIndsWithinDataPartB] = max(tracesB(crossIndsRangeMatB), [], 1);
    crossIndsA = crossIndsWithinDataPartA' + crossSubsRangeMatA(1,:)' - 1;
    crossIndsB = crossIndsWithinDataPartB' + crossSubsRangeMatB(1,:)' - 1;
    
    % Match each spike from MClust with a stim from the original neuron struct
    resultsA = zeros(numStimsPerExp, 1);
    resultsB = zeros(numStimsPerExp, 1);
    
    % Ensure that all crossInds fall inside of the timespan of spikes being
    % accepted as results of a stimulation (timespanLim ms after stim
    % onset)
    crossIndsRemoveA = crossIndsA > (timespanLim/timeGranularity); % Indices to be removed
    crossIndsRemoveB = crossIndsB > (timespanLim/timeGranularity); % Indices to be removed
    
    resultsDataIndsSpikesA(crossIndsRemoveA) = [];
    crossIndsA(crossIndsRemoveA) = [];
    stimNumsA(crossIndsRemoveA) = [];
    resultsDataIndsSpikesB(crossIndsRemoveB) = [];
    crossIndsB(crossIndsRemoveB) = [];
    stimNumsB(crossIndsRemoveB) = [];
    
    sortCodesA(sub2ind(size(sortCodesA), resultsDataIndsSpikesA, crossIndsA)) = 1;
    sortCodesB(sub2ind(size(sortCodesB), resultsDataIndsSpikesB, crossIndsB)) = 1;
    resultsA(stimNumsA) = 1;
    resultsB(stimNumsB) = 1;
    
    % Insert the results from the t-file into the resultsData from that
    % .mat file
    neuronStruct{1}.resultsData(resultsDataIndsA,3) = resultsA;
    neuronStruct{2}.resultsData(resultsDataIndsB,3) = resultsB;
end

% Populate the pcaSortCodes field
neuronStruct{1}.pcaSortCodes = mat2cell(sortCodesA, ones(neuronStruct{1}.resultsDataAllocated, 1), sizeData);
neuronStruct{2}.pcaSortCodes = mat2cell(sortCodesB, ones(neuronStruct{2}.resultsDataAllocated, 1), sizeData);

%% Save this file as a new neuronStruct file
fprintf('Writing .mat file...\n');
save([origFilePath filesep 'MC_' origFileName '.mat'], 'neuronStruct');
fprintf('Done\n\n');