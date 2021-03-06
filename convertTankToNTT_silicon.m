function convertTankToNTT_silicon(varargin)
% This function converts tank information (generated by the Acute Control
% TDT circuit) into an .ntt file to be read by MClust. It also stores
% metadata in the 'base' workspace.
%
% convertTankToNTT_silicon()
% convertTankToNTT_silicon(dateStr)
% convertTankToNTT_silicon(dateStr, writeNTT)
% convertTankToNTT_silicon(dateStr, writeNTT, threshold)
% convertTankToNTT_silicon(dateStr, writeNTT, threshold, timeSpanLim)Z
% convertTankToNTT_silicon(dateStr, writeNTT, threshold, timeSpanLim, useTetrode)
% convertTankToNTT_silicon(dateStr, writeNTT, threshold, timeSpanLim, useTetrode, numStimsPerExp)
%
% dateStr - A string of numbers representing the date on which the
% experiment was performed (e.g. '20170811').
%
% timeStrs - A cell array of strings of numbers, representing the names of
% the .mat files that contain data corresponding to each trial, in
% chronological order (e.g. {'1302', '1646'}).  (Multiple experiments on
% the same unit pair are ok if in different .mat files.)
%
% threshold - A scalar that represents the threshold in uV for spike
% detection (e.g. 35).
%
% numStimsPerExp - A scalar that represents the number of stimulations that
% are expected during an experiment (equal to stimsPerSequence x
% sequencesPerRun x runsPerExperiment = 5 x 20 x 10 = 1000 for most).
%
% writeNTT - An optional logical that allows the user to not write an NTT
% file during this function call, instead simply populating the base
% workspace with the supplemental variables (i.e. 'expTankNums',
% 'expTimeStrs', etc.) associated with this date's experiment.
%
% timeSpanLim - An optional scalar that represents the maximum number of
% milliseconds after stimulus onset to look for a spike event.

dateStr = '';
writeNTT = true;
threshold = 30;
timespanLim = 30; % Acquire all threshold crossings between the onset of stimulation and this value (in ms) after
useOnlineSortCodesForSpikes = false;
numStimsPerExp = 2000;
postStimTraceLengthInSec = .05; % The length of the post-stim trace in seconds (50ms)
postStimResultsLengthInSec = .03; % The time after which we will accept a sort code to count as a spike

isDebugging = true; % CHANGE ASAP

% Filter design
% The silicon probe circuit saves data as unfiltered traces, so we should
% filter the traces before we analyze them.
order = 8; % The total order of the bandpass filter (4th order lowpass and 4th order highpass)
highCutoff = 6000; % The high cutoff for the bandpass filter
lowCutoff = 600; % The low cutoff for the bandpass filter
samplingFreq = 24414.062500; % The sampling rate in Hz

%%
% TODO: MAKE SURE THIS IS FALSE FOR PRODUCTION
scriptTest = false; % Used for testing the script. Will make it so that an "experiment" made of any number of stimulations will be accepted.
%%

if nargin > 0
    if ischar(varargin{1})
        dateStr = varargin{1};
    else
        error('First parameter must be a string, if used.');
    end
end

if nargin > 1
    if islogical(varargin{2})
        writeNTT = varargin{2};
    else
        error('Second parameter must be a logical, if used.');
    end
end

if nargin > 2
    if isnumeric(varargin{3})
        threshold = varargin{3};
    else
        error('Third parameter must be a double, if used.');
    end
end

if nargin > 3
    if isnumeric(varargin{4})
        timespanLim = varargin{4};
    else
        error('Fourth parameter must be a double, if used.');
    end
end

if nargin > 4
    if islogical(varargin{5})
        useOnlineSortCodesForSpikes = varargin{5};
    else
        error('Fifth parameter must be a logical, if used.');
    end
end

if nargin > 5
    if isnumeric(varargin{6})
        numStimsPerExp = varargin{6};
    else
        error('Sixth parameter must be a double, if used.');
    end
end

if ~writeNTT
    fprintf('Not writing NTT files.\n\n');
end

% Parameters
% threshold = 35; % uV, spiking threshold
peakDistance = 0; % Minimum absolute distance between spikes, in ms
numTetChannels = 32;
numSamples = 32; % Number of samples in each waveform to store
numSamplesPreCrossing = 10; % Number of samples before the threshold crossing to store
convToMicroVolts = 100; % Factor to multiply the TDT data by to convert it to microvolts (data is measured in Volts and has a 10,000x gain)
microVoltsToMilliVolts = 1/1000; % Factor to multiply by to convert microvolts to volts

% Build the full file name
path = 'D:\samuelgb\Documents\MATLAB\Acute Control\Data\Silicon Probe Data';
saveFileName = 'convertNTTResults.mat'; % Name of .mat file that data will be saved to (because this takes a damn long time...)

if strcmp(dateStr, '')
    % The user may want to convert ALL data.  Confirm first, then run
    resp = questdlg('No date was provided.  Convert all data in data directory?', 'All data conversion confirmation', 'Yes', 'No', 'No');
    if strcmp(resp, 'No')
        fprintf('No conversions performed\n');
        return;
    else
        % Go through all data, and convert
        isDataDir = @(a) ~isempty(regexp(a, '^\d{8}$', 'start')); % Is "a" formatted as a data directory (e.g. "20170723")
        cd(path);
        dirContents = dir;
        
        % Go through all data directories
        for dirInd = 1:length(dirContents)
            thisDirName = dirContents(dirInd).name;
            if isDataDir(thisDirName)
                % Check if there are any .ntt files already
                if ~any(size(dir([path thisDirName '/*.ntt' ]),1))
                    % If this is a data directory, run convertTankToNTT
                    thisDateStr = thisDirName;
                    try
                        convertTankToNTT_silicon(thisDateStr);
                    catch
                        fprintf('\n\n**********\nCould not convert %s to tank.  Continuing...\n**********\n\n', thisDateStr);
                    end
                end
            end
        end
        
        return;
    end
end

% Prepare data
cmpPartStr = @(strLng, str, extraChars) (length(strLng) >= (length(str) + extraChars)) && strcmp(strLng(1:length(str)), str);
dirNamePrefix = sprintf('Acute_%s_', dateStr); % Tank name prefix
dirNamePrefixFirstTank = sprintf('Acute_%s', dateStr); % Tank name prefix (if it is the only tank, and therefore does not have any subscripts after it)
% dirNamePrefixLength = length(dirNamePrefix); % Length of the tank name prefix
matFilePrefix = sprintf('neuronStruct%s', dateStr); % First part of he .mat file names
% matFilePrefixLength = length(matFilePrefix); % Length of the mat file name prefix
% numExps = length(timeStrs); % Number of experiments
numSamplesPostCrossing = numSamples - numSamplesPreCrossing;
% peakDistance = peakDistance/1000; % Convert peak distance to milliseconds (TDT data is stored in seconds)

% Determine how many tanks and how many .mat files there are
dataDir = [path filesep dateStr];

% Check existence of data directory
if ~exist(dataDir, 'dir')
    fprintf('Directory does not exist. Check the data directory date and name\n');
    return;
end
timeStrs = cell(1); % Cell array in which time strings will be storeds
timeStrsInd = 1;
cd(dataDir);
dirContents = dir; % Get directory contents
numTanks = 0;
fileExt = '.mat';
for fileNum = 1:length(dirContents)
    thisName = dirContents(fileNum).name;
    if cmpPartStr(thisName, dirNamePrefix, 2)
        numTanks = max(numTanks, str2double(thisName((end - 1):end)));
        continue;
    end
    
    if cmpPartStr(thisName, dirNamePrefixFirstTank, 0)
        numTanks = 1; % This represents the first tank
        newName = [thisName '_01']; % Generate a new, suffixed name
        movefile(thisName, newName); % Rename the file
        
        % Add the suffix to everything in the directory
        cd(newName);
        subDirNamePrefix = sprintf('TutTank_%s', dirNamePrefixFirstTank);
        subDirContents = dir;
        for subFileNum = 1:length(subDirContents)
            thisSubName = subDirContents(subFileNum).name;
            if cmpPartStr(thisSubName, subDirNamePrefix, 4)
                subParts = strsplit(thisSubName, '.');
                subParts{1} = [subParts{1} '_01']; % Generate a new, suffixed name
                newFullName = [subParts{1} '.' subParts{2}];
                movefile(thisSubName, newFullName); % Rename the file
            end
        end
        
        cd('..')
        addpath(newName); % Add this folder name to the path
        continue;
    end
    
    if cmpPartStr(thisName, matFilePrefix, 9) && strcmp(thisName((end - length(fileExt) + 1):end), fileExt)
        % Check that this file actually contains completed trials
        clear neuronStruct
        
        try
            load(thisName);
        catch
            continue;
        end
        
        if ~exist('neuronStruct', 'var')
            continue;
        end
        
        if ~iscell(neuronStruct) || ~isstruct(neuronStruct{1})
            continue;
        end
        
        if ~isfield(neuronStruct{1}, 'trialData')
            continue;
        end
        
        %         numCompletedTrials = 0;
        %         for numTrial = 1:length(neuronStruct{1}.trialData)
        %             if ~isempty(neuronStruct{1}.trialData(numTrial).trialCompleted)
        %                 numCompletedTrials = numCompletedTrials + neuronStruct{1}.trialData(numTrial).trialCompleted;
        %             end
        %         end
        %
        %         if numCompletedTrials == 0 || mod(numCompletedTrials, 2) ~= 0
        %             continue;
        %         end
        
        % If this file has recorded some even number of trials > 2, then
        % record this time string to have its information loaded later
        fileParts1 = strsplit(thisName, '.');
        fileParts2 = strsplit(fileParts1{1}, '_');
        timeStrs{timeStrsInd} = fileParts2{2};
        timeStrsInd = timeStrsInd + 1;
    end
    
end

if numTanks == 0
    fprintf('No data tanks of form %s_xx found in %s\n', dirNamePrefix, pwd);
    return;
end

curExpNum = 1;
usedTetrodeProbe = true;
expTankNums = []; % Will contain information on which tank each experiment came from
expTimeStrs = {}; % Will contain information on what each experiment is named (the timeString) and what number experiment it is for this date
expTimeBounds = []; % Time boundaries for this experimental data
expPenetrationNums = []; % Will contain information on which penetration numbers (in the TDT data) were complete for each experiment
expStimResponses = []; % Will contain the stimulation response traces, 50ms traces that start at the stimulation onset
expExpectedResults = []; % Will contain the expected results from each stimulation (first column is if unit 1 is expected, second column for unit 2)
expResults = []; % Will contain the results of the stimulations judged by the presence of sort-codes
expSortTimes = []; % Will contain the sort code responses as a function of time
expTetNums = []; % Will contain the tetrode numbers of each unit (N x 2, where expTetNums(n, x) gives the tetrode number of unit x (1 or 2) in experiment n)
% expElectrodeTraces = []; 
for curTankNum = 1:numTanks
    tankSubName = sprintf('Acute_%s_%0.2d', dateStr, curTankNum);
    fileName = sprintf('%s%s%s%s%s%sTutTank_%s.tsq', path, filesep, dateStr, filesep, tankSubName, filesep, tankSubName);
    fprintf('Looking for data from %s...\n', dateStr);
    if ~exist(fileName, 'file')
        fprintf('%s does not exist, skipping...\n', fileName);
        continue;
    end
    
    %% Start looking for data in this tank
    tic; % When extraction began
    block = tdt_block_sam(fileName); % Create a block object
    
    %% Determine if a trial was run in this tank
    try
        pNum = getdata(block, 'TDat', 'channel', 4);
    catch
        fprintf('%s does not contain any complete experiments, skipping...\n', fileName);
        toc;
        continue;
    end
    
    %% Determine which experiments were completed, and how many there are
    [numStims, penNums] = hist(pNum.vals, length(unique(pNum.vals)));
    if ~scriptTest
        penNums = penNums(numStims >= numStimsPerExp); % Narrow down the penetration numbers to include only those which represent a completed experiment (the initial GT optimization will create a number of "trials" that have the same id as the main trial; only use the last numStimsPerExp number of stimulations)
    end
    thisNumExp = length(penNums);
    
    %% Extract data from tank
    curPenNumLoc = 1;
    while curPenNumLoc <= thisNumExp
        %% Get the timestamps of the stimulations
        relevantStimInds = find(pNum.vals == penNums(curPenNumLoc), numStimsPerExp, 'last');
        stimTimeStamps = pNum.times(relevantStimInds);
        timeBoundsA = [0 timespanLim]/1000;
        timeBoundsB = [0 timespanLim]/1000;
        
        % Get this penetration number
        [thisNumStims, thisPenNums] = hist(pNum.vals(relevantStimInds), length(unique(pNum.vals(relevantStimInds))));
        if ~scriptTest
            thisPenNums = thisPenNums(thisNumStims == numStimsPerExp); % Narrow down the penetration numbers to include only those which represent a completed experiment
        end
        
        %% Get waveform data from the tank
        timeBounds = [max([0 stimTimeStamps(1) - 1]) (stimTimeStamps(end) + 1)]; % Get the beginning and end times of the experiment (add to the end time, because it indicates the stim time, the neural response takes longer to develop, and buffer second before the start as well)
        %         timeBounds = [0 50]; % For testing
        
        usedTetrodeProbe = true;
        
        % Load in all of the channels
        fprintf('Loading all electrode data...\n');
        tic;
        
        %         [B, A] = butter(order/2, [highCutoff lowCutoff]/(samplingFreq/2)); % Design a filter for the data
        bandPassFilt = designfilt('bandpassiir', 'FilterOrder', order, 'HalfPowerFrequency1', lowCutoff, 'HalfPowerFrequency2', highCutoff, 'SampleRate', samplingFreq);
        testData = getdata(block, 'Mul1', 'channel', 1, 'interval', timeBounds); % Get an example piece of data, so that the electrode data matrix can be shaped properly
        allElectrodeData = zeros(numTetChannels,length(testData.vals));
        for chan = 1:(numTetChannels/2)
            dat1 = getdata(block, 'Mul1', 'channel', chan, 'interval', timeBounds);
            dat2 = getdata(block, 'Mul2', 'channel', chan, 'interval', timeBounds);
            
            allElectrodeData(chan, :) = filtfilt(bandPassFilt, double(dat1.vals));
            allElectrodeData(chan + (numTetChannels/2), :) = filtfilt(bandPassFilt, double(dat2.vals));
            
            if mod(chan, (numTetChannels/8)) == 0
               fprintf('Finished %d channels out of %d...\n', chan, numTetChannels/2); 
            end
        end
        fprintf('Done after %0.2f seconds!\nAnalyzing...\n', toc);
        
        % Extract the channel number(s) that correspond to each unit, so we
        % know which waveforms to use
        if usedTetrodeProbe
            % Assume that the tetrode number doesn't actually get
            % changed over the course of the recording
            unit1TetStruct = getdata(block, 'eNum', 'channel', 1);
            unit2TetStruct = getdata(block, 'eNum', 'channel', 3);
            
            unit1Tet = unit1TetStruct.vals(1);
            unit2Tet = unit2TetStruct.vals(1);
            
            % Get the channel numbers of each tetrode
            unit1Chans = ((unit1Tet - 1)*4 + 1):((unit1Tet - 1)*4 + 4);
            unit2Chans = ((unit2Tet - 1)*4 + 1):((unit2Tet - 1)*4 + 4);
        else
            % TODO: Fill this in when I get the linear probes ready
            % Construct unit1Chans and unit2Chans by getting the singular
            % channel that each unit is defined by
            
        end
        
        % Get tetrode sort codes
        sortCodesA = getdata(block, 'Sort', 'channel', 1, 'interval', timeBounds); % Get the sort codes for unit 1
        sortCodesB = getdata(block, 'Sort', 'channel', 2, 'interval', timeBounds); % Get the sort codes for unit 2
        
        % The sortcodes from TetSort come as a 16-bit mask, indicating if each
        % frame has a waveform in any of the sorting circles (see documentation
        % for more details)
        
        % Convert the values from doubles to uint16's
        sortCodesAInt = uint16(sortCodesA.vals);
        sortCodesBInt = uint16(sortCodesB.vals);
        
        % First, break sortCodes down into bits (make an N x L
        % matrix, where L is the length of the incoming codes,
        % where each value in the matrix represent the nth bit
        % of the lth code)
        N = 16;
        allSortBitsA = bitget(repmat(sortCodesAInt', N, 1), repmat((1:N)', 1, length(sortCodesAInt)));
        allSortBitsB = bitget(repmat(sortCodesBInt', N, 1), repmat((1:N)', 1, length(sortCodesBInt)));
        
        % Find the least significant bit of each code (the
        % "first" sort code that each waveform is sorted into).
        %  We are throwing out some sorting information, but we
        %  don't really need to know EVERY sort code that each
        %  waveform falls into
        [isValidA, finalCodesA] = max(allSortBitsA, [], 1);
        [isValidB, finalCodesB] = max(allSortBitsB, [], 1);
        
        % Make sure that we don't count any zero values
        isValidA = logical(isValidA); % Any value where the max is NOT 0
        isValidB = logical(isValidB);
        finalCodesA(~isValidA) = 0;
        finalCodesB(~isValidB) = 0;
        
        % Store the final sort codes
        thisSortCodesA = double(finalCodesA);
        thisSortCodesB = double(finalCodesB);
        
        % Finally, get the timestamps of the "1" sortcodes
        thisSortTimesA = sortCodesA.times(thisSortCodesA == 1);
        thisSortTimesB = sortCodesB.times(thisSortCodesB == 1);
        
        % Get threshold crossing indices
        if useOnlineSortCodesForSpikes
            % Use the PCA sort codes to find spikes
            % TODO: Process the sort codes first if doing this, because they only represent a 16-bit mask!!!
            crossIndsA = find(sortCodesA.vals == 1);
            crossIndsB = find(sortCodesB.vals == 1);
        else
            % Use peak crossings to find spikes
            timeGranularity = diff(testData.times([1 2])); % Time between two points
            
            crossIndsA = [];
            crossIndsB = [];
            for chanNum = 1:4
                % Go through each channel in the tetrode
                [~, thisCrossIndsA] = findpeaks(allElectrodeData(unit1Chans(chanNum), :)*convToMicroVolts, 'minpeakheight', threshold, 'minpeakdistance', round(peakDistance/timeGranularity));
                [~, thisCrossIndsB] = findpeaks(allElectrodeData(unit2Chans(chanNum), :)*convToMicroVolts, 'minpeakheight', threshold, 'minpeakdistance', round(peakDistance/timeGranularity));
                
                % Find only new threshold crossings
                % TODO: Use uniquetol, to give a tolerance on unique
                % indices?
                crossIndsA = unique([crossIndsA;thisCrossIndsA']);
                crossIndsB = unique([crossIndsB;thisCrossIndsB']);
            end
        end
        
        % Use the crossInds to get the corresponding times in the electrode
        % waveforms
        timeStampsA = testData.times(crossIndsA);
        timeStampsB = testData.times(crossIndsB);
        
        % Find the threshold crossing times of each spike
        binEdgesMatA = bsxfun(@plus, timeBoundsA,stimTimeStamps)'; % A matrix containing the start and end times of each stim response period, in TDT time
        binEdgesMatB = bsxfun(@plus, timeBoundsB,stimTimeStamps)'; % A matrix containing the start and end times of each stim response period, in TDT time
        
        binEdgesA = sort(binEdgesMatA(:)); % The bin edges of A.  Every other pair of numbers (inds 1 and 2, 3 and 4, 5 and 6, etc.) indicate the beginning and end of the stim response period, in which a spike is considered a stim response.
        binEdgesB = sort(binEdgesMatB(:)); % The bin edges of B.  Every other pair of numbers (inds 1 and 2, 3 and 4, 5 and 6, etc.) indicate the beginning and end of the stim response period, in which a spike is considered a stim response.
        
        [~, binNumsA] = histc(timeStampsA, binEdgesA); % Find where each stimulus is
        [~, binNumsB] = histc(timeStampsB, binEdgesB); % Find where each stimulus is
        
        binNumsA(binNumsA == 0) = max(binNumsA);
        binNumsB(binNumsB == 0) = max(binNumsB);
        
        binNumsA(mod(binNumsA,2) == 0) = binNumsA(mod(binNumsA,2) == 0) - 1;
        binNumsB(mod(binNumsB,2) == 0) = binNumsB(mod(binNumsB,2) == 0) - 1;
        
        stimNumsA = ceil(binNumsA/2); % Get the stimulation number for each spike
        stimNumsB = ceil(binNumsB/2); % Get the stimulation number for each spike
        
        spikeDelayA = (timeStampsA - binEdgesA(binNumsA))*1000; % The start of the small range in which the maximum of a spike should be found
        spikeDelayB = (timeStampsB - binEdgesB(binNumsB))*1000;
        
        spikeDelayA(spikeDelayA > timespanLim) = NaN;
        spikeDelayB(spikeDelayB > timespanLim) = NaN;
        
        %% Find the stimulus response traces for each stimulation (not needed in non-tetrode analysis, because Matlab stores one channel's worth of trace information)
        expStimResponses{curExpNum} = cell(numStimsPerExp, 2); %#ok<AGROW> % Initialize a cell for this experiment's stim responses
        expResults{curExpNum, 1} = false(numStimsPerExp, 2); %#ok<AGROW> % Initialize a matrix to hold the results for each stimulation
        for stimNum = 1:numStimsPerExp
            startTime = stimTimeStamps(stimNum);
            endTimeTrace = startTime + postStimTraceLengthInSec;
            endTimeResults = startTime + postStimResultsLengthInSec;
            
            % Get the response traces
            startInd = find(testData.times >= startTime, 1); % Find the first index later than the stimulus onset
            endInd = find(testData.times >= endTimeTrace, 1); % Find the first index greater than postStimTraceLengthInSec seconds later than the stimulus onset
            expStimResponses{curExpNum}{stimNum, 1} = allElectrodeData(:, startInd:endInd)'; % Get the voltage traces for all channels
            expStimResponses{curExpNum}{stimNum, 2} = testData.times(startInd:endInd); % Get the time stamps for all channels
            
            % Get the results (spike or no-spike)
            expResults{curExpNum}(stimNum, :) = [any(startTime <= thisSortTimesA & thisSortTimesA <= endTimeResults) any(startTime <= thisSortTimesB & thisSortTimesB <= endTimeResults)];
        end
        
        % Save the expected results that TDT stimulated based on
        eRStruct = getdata(block, 'TDat', 'channel', 5, 'interval', timeBounds);
        expExpectedResults{curExpNum} = [~eRStruct.vals eRStruct.vals]; %#ok<AGROW> % Initialize a matrix to hold the expected results of each stimulation
        
        %% Generate .ntt files from the timestamps and waveforms
        if curExpNum <= length(timeStrs)
            timeString = timeStrs{curExpNum};
        else
            timeString = fprintf('exp%0.2d', curExpNum);
        end
        
        if writeNTT
            % Get waveform windows
            traceDataStartA = crossIndsA - numSamplesPreCrossing; % Starting point for the trace window
            traceDataEndA = crossIndsA + numSamplesPostCrossing; % Ending point for the trace window
            traceDataStartB = crossIndsB - numSamplesPreCrossing; % Starting point for the trace window
            traceDataEndB = crossIndsB + numSamplesPostCrossing; % Ending point for the trace window
            
            % Correct windows that overlap with the beginning of the experiment
            traceDataStartA(traceDataStartA < 1) = 1; % Set the beginning of the window to 1
            traceDataStartB(traceDataStartB < 1) = 1; % Set the beginning of the window to 1
            
            % Correct windows that overlap with the end of the experiment
            traceDataStartA(traceDataEndA > size(allElectrodeData, 2)) = size(allElectrodeData, 2) - numSamples + 1; % Set the beginning of the window so that it ends at the end of the experiment
            traceDataStartB(traceDataEndB > size(allElectrodeData, 2)) = size(allElectrodeData, 2) - numSamples + 1; % Set the beginning of the window so that it ends at the end of the experiment
            
            % Get indices of the threshold crossing waveforms with respect to
            % the fullProbeX.vals vector
            waveformIndsA = bsxfun(@plus, traceDataStartA(:),0:(numSamples - 1));
            waveformIndsB = bsxfun(@plus, traceDataStartB(:),0:(numSamples - 1));
            
            % Extract the waveforms, and separate them by channel
            waveformsA = [];
            waveformsB = [];
            for channel = 1:4
                thisChan1 = allElectrodeData(unit1Chans(channel), :);
                thisChan2 = allElectrodeData(unit2Chans(channel), :);
                
                %                 waveformsA = cat(3, waveformsA, thisChan1(waveformIndsA)*microVoltsToMilliVolts);
                %                 waveformsB = cat(3, waveformsB, thisChan2(waveformIndsB)*microVoltsToMilliVolts);
                
                % Record the data in 100ths of microvolts (it was too
                % heavily quantized in millivolts or microvolts due to the
                % conversion to int16)
                %                 waveformsA = cat(3, waveformsA, thisChan1(waveformIndsA));
                %                 waveformsB = cat(3, waveformsB, thisChan2(waveformIndsB));
                
                % Record the data in 10ths of microvolts (it seems to be
                % too large if recorded in 100ths of microvolts, and maxes
                % out the int16 value
                waveformsA = cat(3, waveformsA, thisChan1(waveformIndsA)/10);
                waveformsB = cat(3, waveformsB, thisChan2(waveformIndsB)/10);
                
                
            end
            
            % Permute the waveform matrices to (numWaveforms x 4
            % waveformLength)
            waveformsA = permute(waveformsA, [1 3 2]);
            waveformsB = permute(waveformsB, [1 3 2]);
            
            if usedTetrodeProbe
                write_waveforms_ntt(sprintf('%s_%s_A.ntt', dateStr, timeString), timeStampsA, waveformsA);
                write_waveforms_ntt(sprintf('%s_%s_B.ntt', dateStr, timeString), timeStampsB, waveformsB);
            else
                write_waveforms_nse(sprintf('%s_%s_A.nse', dateStr, timeString), timeStampsA, waveformsA);
                write_waveforms_nse(sprintf('%s_%s_B.nse', dateStr, timeString), timeStampsB, waveformsB);
            end
            
            if isDebugging
               save(sprintf('convertTankToNTT_DEBUGTEST_%s.mat', timeString), '-v7.3');
            end
        end
        
        %% Add metadata to the 'base' workspace
        expTankNums(curExpNum) = curTankNum; %#ok<AGROW>
        expTimeStrs{curExpNum} = timeString; %#ok<AGROW>
%         expElectrodeTraces{curExpNum} = allElectrodeData; %#ok<AGROW>
        expPenetrationNums(curExpNum, :) = thisPenNums; %#ok<AGROW> % Narrow down the penetration numbers to include only those which represent a completed experiment
        expTimeBounds(curExpNum, :) = timeBounds; %#ok<AGROW>
        expTetNums(curExpNum, :) = [unit1Tet unit2Tet]; %#ok<AGROW>
        expSortTimes{curExpNum} = {thisSortTimesA, thisSortTimesB}; %#ok<AGROW>
        
        assignin('base', 'triggerVoltage', threshold);
        %         assignin('base', 'spikeDelay', []);
        assignin('base', 'trialTimeBounds', timeBounds)
        assignin('base', sprintf('spikeDelay_%s_A', timeString), spikeDelayA);
        assignin('base', sprintf('spikeDelay_%s_B', timeString), spikeDelayB);
        assignin('base', sprintf('stimNums_%s_A', timeString), stimNumsA);
        assignin('base', sprintf('stimNums_%s_B', timeString), stimNumsB);
        
        %% Iterate the location in penNums and timeStrs
        curPenNumLoc = curPenNumLoc + 1;
        curExpNum = curExpNum + 1;
    end
    
    fprintf('Finished reading from tank %d of %d\n', curTankNum, numTanks);
    toc;
end

%% Add metadata to the 'base' workspace
assignin('base', 'expTankNums', expTankNums);
assignin('base', 'expTimeStrs', expTimeStrs);
assignin('base', 'expTetNums', expTetNums);
assignin('base', 'expStimResponses', expStimResponses);
assignin('base', 'expResults', expResults);
assignin('base', 'expExpectedResults', expExpectedResults);
assignin('base', 'expSortTimes', expSortTimes);
assignin('base', 'expTimeBounds', expTimeBounds);
assignin('base', 'expPenetrationNums', expPenetrationNums);
assignin('base', 'timespanLim', timespanLim);
assignin('base', 'useTetrode', usedTetrodeProbe);
assignin('base', 'isSilicon', true);
assignin('base', 'expDateStr', dateStr);
% assignin('base', 'expElectrodeTraces', expElectrodeTraces);
% assignin('base', 'waveformsA', waveformsA);
% assignin('base', 'waveformsB', waveformsB);

evalin('base', sprintf('save(''%s'', ''expTankNums'', ''expTimeStrs'', ''expTetNums'', ''expStimResponses'', ''expResults'', ''expExpectedResults'', ''expSortTimes'', ''expTimeBounds'', ''expPenetrationNums'', ''timespanLim'', ''useTetrode'', ''isSilicon'', ''expDateStr'', ''-v7.3'');', saveFileName));

% Save all of the spike delays and stimulation numbers
for saveExpInd = 1:(curExpNum - 1)
   evalin('base', sprintf('save(''%s'', ''spikeDelay_%s_A'', ''spikeDelay_%s_B'', ''stimNums_%s_A'', ''stimNums_%s_B'', ''-append'');', saveFileName, expTimeStrs{saveExpInd}, expTimeStrs{saveExpInd}, expTimeStrs{saveExpInd}, expTimeStrs{saveExpInd}));
end

fprintf('\nFinished converting to .ntt\n');

if curExpNum <= length(timeStrs)
    fprintf('Could not find data for %d experiment timestring(s)\n', curExpNum - length(timeStrs) + 1);
end
end