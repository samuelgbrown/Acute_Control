function [allMatCell, timeStrings, dateStrings] = loadTrialData(varargin)
% [allMatCell, timeStrings, dateStrings] = loadTrialData
% This function will load all files that were output by the TDT_Control
% program, and will load each neuronStruct into the allMatCell cell array.
% timeStrings returns a cell array of strings that represent the times at
% which the files were created, and dateStrings returns a cell array of
% strings representing the dates that the files were created.  If called in
% a data directory, all files will be loaded.  If called in a directory in
% which data directories are sub-directories, all data in each data
% directory will be loaded.
%
% [allMatCell, numFiles,timeStrings] = loadTrialData(filesToLoad)
% loadMCFiles is a length = 2 logical vector which indicates whether or not
% the original sorts (filesToLoad(1)) and the offline sorts
% (filesToLoad(2)) should be loaded.  These flags must be true IN ADDITION
% to the metaDataFile indicating that a file should be loaded for it to be
% loaded by this file.
%
% [allMatCell, numFiles,timeStrings] = loadTrialData(filesToLoad, readAnyMatName)
% If the logical readAnyMatName is set to true, then all .mat files will be
% read (this allows the user to rename any files that were output from
% TDT_Control, for example when re-exporting old data.
%
% [allMatCell, numFiles,timeStrings] = loadTrialData(filesToLoad, readAnyMatName, dirsToRead)
% If this function is called in a directory with many data directories in
% it (e.g. "20170723"), it will go through each directory, check if it is a
% data directory, and extract data from it.  dirsToRead is a cell array
% that can be populated with the names of directories as strings that
% contain data to be analyzed.  dirsToRead may also be a 0, indicating that
% all directories should be read.
%
% [allMatCell, numFiles,timeStrings] = loadTrialData(filesToLoad, readAnyMatName, dirsToRead, loadOnlyPairs)
% If loadOnlyPairs is set to true, then data for a given date will only be
% loaded if both the online and the offline sorting data exists.  This
% being true overrules the settings for filesToLoad.

%% Assign default output
allMatCell = cell(1);
timeStrings = cell(1);
dateStrings = cell(1);

%% Sort input
% If there are no inputs, set to defaults (read all standard named
% files)
% filesToLoad = 0;
filesToLoad = false(1,2);
readAnyMatName = false; % Try reading .mat files of any name (user may have renamed the mat file, or consolidated one from earlier versions of the program)
dirsToRead = 0;
loadOnlyPairs = false; % Should only offline/online pairs be loaded?
respectDataToAnalyze = true; % If the dataToAnalyze .txt is in a data directory, should we respect the contents or ignore them?
if nargin >= 1
    filesToLoad = varargin{1};
end
if nargin >= 2
    readAnyMatName = varargin{2};
end
if nargin >= 3
    dirsToRead = varargin{3};
end
if nargin >= 4
    loadOnlyPairs = varargin{4};
end
if nargin >= 5
    respectDataToAnalyze = varargin{4};
end

if nargin >= 6
    error('Too many inputs')
end

%% Prepare analysis
isDataDir = @(a) ~isempty(regexp(a, '^\d{8}$', 'start')); % Is "a" formatted as a data directory (e.g. "20170723")
shouldAnalyzeThisDir = @(a) any(cell2mat(cellfun(@(c) strcmp(a, c), dirsToRead, 'uniformoutput', false))); % Is directory name "a" part of the string cell array "dirsToRead"
metaDataFileName = 'dataToAnalyze.txt'; % The name of the text file that shows which .mat files to load
matName = 'neuronStruct';
readAllDirs = ~iscell(dirsToRead) && dirsToRead == 0; % Should all directories be read? (Is dirsToRead equal to 0?)
loadOnlineSorts = filesToLoad(1);
loadOfflineSorts = filesToLoad(2);

%% Load all files
% Prepare to access all data
allMatCellInd = 1;

originalPath = pwd;
pathParts = strsplit(originalPath, filesep);
dirName = pathParts{end};
thisIsDataDir = isDataDir(dirName);
if thisIsDataDir
    numDirs = 1;
else
    upperDirContents = dir; % Get this directory's contents
    numDirs = length(upperDirContents);
end
for upperDirInd = 1:numDirs
    if ~thisIsDataDir
        fprintf('Reading directory %d out of %d...\n', upperDirInd, numDirs);
        thisDir = upperDirContents(upperDirInd);
    end
    if thisIsDataDir || isDataDir(thisDir.name)
        % This is a data directory (or the upper directory is a data
        % directory
        if ~thisIsDataDir
            thisDirName = sprintf('%s%s%s', thisDir.folder, filesep, thisDir.name);
            cd(thisDirName);
        end
        pathParts = strsplit(pwd, filesep);
        dateStr = pathParts{end};
        if thisIsDataDir || readAllDirs || shouldAnalyzeThisDir(dateStr)
            % If this directory should be analyzed (dirsToRead is 0, this
            % directory is on the whitelist)
            dirContents = dir; % Get directory contents
            
            fullMatName = sprintf('%s%s', matName, dateStr);
            fileNum = 1;
            
            % Read the dataToAnalyze file to determine which files should
            % be read
            dataToAnalyzeFileName = [pwd filesep metaDataFileName]; % The name of the data to analyze file (if it exists)
            %             if ((length(fileName) >= length(metaDataFileName)) && strcmp(fileName(1:length(metaDataFileName)), metaDataFileName))
            if respectDataToAnalyze && exist(dataToAnalyzeFileName, 'file')
                % Load this meta data
                fid = fopen(dataToAnalyzeFileName);
                thisFilesToLoad = cell2mat(cellfun(@str2double, strsplit(fgetl(fid), ','), 'uniformoutput', false));
                thisMCFilesToLoad = thisFilesToLoad; % str2num(fgetl(fid));
                fclose(fid);
                
                if (length(thisFilesToLoad) == 1 && (thisFilesToLoad == -1))% && (length(thisMCFilesToLoad) == 1 && (thisMCFilesToLoad == -1))
                    % No files in this directory should be read, so skip
                    % it.
                    continue;
                end
            else
                thisFilesToLoad = 0; % If there is a metaDataFile in this directory, then use the input filesToLoad
                thisMCFilesToLoad = 0; % If there is a metaDataFile in this directory, then use the input filesToLoad
            end
            
            % Load all of the mat files for this day
            for dirInd = 1:length(dirContents)
                fileName = dirContents(dirInd).name;
%                 fprintf('Reading file %d out of %d...\n', dirInd, length(dirContents));
                
                if (((length(fileName) >= length(fullMatName)) && strcmp(fileName(1:length(fullMatName)), fullMatName) && strcmp(fileName((end - 2):end), 'mat')) || readAnyMatName)
                    % For all mat files that were produced during the trial
                    thisTimeString = fileName((end - 7):(end - 4));
                    idealOnlineFileName = sprintf('neuronStruct%s_%s.mat', dateStr, thisTimeString);
                    if ~(strcmp(idealOnlineFileName, fileName) || readAnyMatName)
                        % If this file name is not in the expected format,
                        % and we are not allowing any mat files that are
                        % not in the expected format, then do not read this
                        % file.
                        continue;
                    end
                    
                    offlineFileName = sprintf('MC_neuronStruct%s_%s.mat', dateStr, thisTimeString); % Build the name of the offline sort file
                    
                    fileName = sprintf('%s%s%s', pwd, filesep, fileName); % Make absolute file path
                    offlineFileName = sprintf('%s%s%s', pwd, filesep, offlineFileName);
                    
                    shouldLoadOnline = (loadOnlyPairs || loadOnlineSorts) && (all(thisFilesToLoad == 0) || any(thisFilesToLoad == fileNum));
                    shouldLoadOffline = (loadOnlyPairs || loadOfflineSorts)  && exist(offlineFileName, 'file') && (all(thisMCFilesToLoad == 0) || any(thisMCFilesToLoad == fileNum));
                    if (loadOnlyPairs && shouldLoadOnline && shouldLoadOffline) || (~loadOnlyPairs && (shouldLoadOnline || shouldLoadOffline))
                        % If only offline/online pairs should be loaded,
                        % then both files are required.  Otherwise, either
                        % the online file or the offline file is required.
                        fprintf('Loading file %d out of %d...\n', dirInd, length(dirContents));
                        
                        if shouldLoadOnline
                            %% Load this .mat file
                            load(fileName);
                            if exist('neuronStruct', 'var')
                                for neuronNum = 1:length(neuronStruct) %#ok<USENS>
                                    allMatCell{allMatCellInd, neuronNum, 1} = neuronStruct{neuronNum};
                                end
                            end
                        end
                        
                        if shouldLoadOffline
                            %% Load the offline sorted version of this file, if it exists
                            % Only processed after the online mat file is
                            % considered/processed because thisMCFilesToLoad
                            % numbering is based on online files
                            if exist(offlineFileName, 'file')
                                % Load the offline sort name, and add it to allMatCell
                                load(offlineFileName);
                                if exist('neuronStruct', 'var')
                                    for neuronNum = 1:length(neuronStruct)
                                        allMatCell{allMatCellInd, neuronNum, 2} = neuronStruct{neuronNum};
                                    end
                                end
                            end
                        end
                        
                        if exist('neuronStruct', 'var')
                            % If either an offline or an online file was loaded
                            dateStrings{allMatCellInd} = dateStr;
                            timeStrings{allMatCellInd} = thisTimeString; % Get the string that represents the time at which the file was saved
                            allMatCellInd = allMatCellInd + 1;
                            clear neuronStruct;
                        end
                    end
                    fileNum = fileNum + 1;
                end
            end
        end
    end
    cd(originalPath);
    fprintf('\n');
end

% % Prepare meta data from trials
% if size(allMatCell,2) == 2
%     numFiles = size(allMatCell,1);
% else
%     numFiles = 0;
% end