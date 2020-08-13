function varargout = TDT_Control_Adaptive(varargin)
%TDT_CONTROL M-file for TDT_Control.fig
%      TDT_CONTROL, by itself, creates a new TDT_CONTROL or raises the existing
%      singleton*.
%
%      H = TDT_CONTROL returns the handle to a new TDT_CONTROL or the handle to
%      the existing singleton*.
%
%      TDT_CONTROL('Property','Value',...) creates a new TDT_CONTROL using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to TDT_Control_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      TDT_CONTROL('CALLBACK') and TDT_CONTROL('CALLBACK',hObject,...) call the
%      local function named CALLBACK in TDT_CONTROL.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TDT_Control

% Last Modified by GUIDE v2.5 03-Oct-2019 11:55:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TDT_Control_OpeningFcn, ...
    'gui_OutputFcn',  @TDT_Control_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before TDT_Control is made visible.
function TDT_Control_OpeningFcn(hObject, eventdata, handles, varargin)
% CHECKED
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Matlab initializations
% Choose default command line output for TDT_Control
handles.output = hObject;

% UIWAIT makes TDT_Control wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Initialize parameters
handles.maxStrength = 5; % The maximum allowed strength for stimulation
handles.maxDuration = 50; % The maximum allowed duration for stimulation
handles.maxIrradianceStrength = handles.maxStrength; % The maximum strength allowable by the irradiance condition (constrains the maximum power exerted by the laser)
handles.maxAllowableIrradiance = 80; % (mW/mm^2) Irradiance boundary of 80mW/mm^2
handles.laserRadius = .1; % Radius of wave guide in mm
handles.powToIrr = 1/((handles.laserRadius^2)*pi*1000); % The factor by which to multiply the fiber power in µW (the laser parameters cubic equation converts the control voltage to µW) to get irradiance in mW/mm^2 (calculated by dividing the power by the area of the fiber)
handles.numNeurons = 2; % The number of neurons that are currently being tracked.  DO NOT CHANGE WITHOUT MODIFYING THE ENTIRE PROGRAM
handles.numControlParams = 2; % The number of control parameters for the neurons, (G/strength and T/duration)
handles.numSlots = 5; % The length of the control sequence (number of "slots", or stimulations)
handles.sizeBuffer = 1000; % The number of samples in the buffers retreived from TDT
handles.allNeuronsStruct = initializeAllNeuronStructs(handles.numNeurons); % The struct which will hold all of the currently known information on both neurons
handles.structsLoaded = false(handles.numNeurons,1); % Is the allNeuronStruct variable loaded with real data?
handles.baseFileName = 'neuronStruct'; % The base name for a file created by this program
handles.structVarName = 'neuronStruct'; % The base name for a variable located in the file created by this program
handles.DA = []; % The activeX controller object to communicate with TDT
handles.TDTConnected = false; % Is the program currently connected to TDT?
handles.TDTRecording = false; % Is TDT currently recording?
handles.TDTUsePCA = false; % Is the TDT system using SpikePAC?
handles.optServeConnected = false; % Is the program currently connected to the optimization server?
handles.traceTime = 50; % The length in ms of each voltage trace (can we get this straight from TDT somehow?)
handles.lastSavedNeuronsStruct = cell(handles.numNeurons,1); % The last versions of the structs that were saved (if they are equal to the current version of the struct, no saving will occur)
handles.startTime = clock; % Record the time at which the GUI was started, all neuron spike times will be measured from this time
handles.TDTSampRate = 24414.1; % Hz, the sampling rate of TDT (used to convert "number of samples" to seconds)
% TODO REMOVE % handles.stimTimeMinMax = [1 3]; % The minimum and maximum amount of time between stimulations.  The earliest time for the next pulse (handles.nextStimTime) will be chosen as the current time plus a random value between these numbers
handles.lastStimTime = handles.startTime; % Records the time of the last stimulation using toc
handles.nextStimTime = handles.startTime; % Records the earliest time that the next stimulation may occur
handles.dataLoc = pwd; % Directory where data is being held
handles.neuronColor = [.3 .75 .93;0 .8392 0]; % Colors (in RGB) of the neurons
% TODO REMOVE % handles.exampleSpikeLength = 32; % The number of samples in an example spike
% TODO REMOVE % handles.exampleSpikeCollection = cell(1, handles.numNeurons); % A collection of spikes that represents the spike shape that the user wants (nxm, where n is the number of spikes and m is the length of each voltage trace)
handles.controlSequence = ones(handles.numSlots, 1); % The order of the neurons to stimulate in the control sequence
handles.controlSequence([2 4]) = 2; % Set the sequence to alternate
handles.lastSequenceResults = false(handles.numSlots, handles.numNeurons); % Logical matrix that represents if a neuron fired during the most recent run sequence command
% handles.hasFiredSequence = false; % Has the program fired a sequence yet?
handles.stoppingSequence = false; % Has the user requested that the sequences stop?
handles.stoppingTrial = false; % Has the user requested that the trial stop?
handles.runningTrial = false; % Is the program running a trial right now?
handles.pausingTrial = false; % Has the user requested that the trial that is currently running be paused?
handles.stimEnabled = false; % Is a new stim allowed right now? (Turned off when a previous stim is currently running, or when TDT is not connected/testing isn't on)
handles.unit1IsUnitA = false; % Records which unit is unit A.  According to convention, unit A prefers short/strong light pulses.  This is determined by the ratio of both GT pairs.
handles.doTestingNeuronDrift = false; % Should the test neurons drift over the course of the experiment?

handles.useParallel = false; % Should parallel processing be used, specifically for the theta/GT optimization?
handles.useOptServe = true; % Should the optimization server be used to speed up optimization

handles.lastDrawSDParams = zeros(handles.numNeurons, handles.allNeuronsStruct{1}.numParameters);
handles.lastDrawSDCurve = cell(handles.numNeurons, 1);

% handles = resetAdaptive(handles); % Initialize parameters related to the adaptive optimizer
% handles.costFun = @(avgProb1, avgProb2) avgProb2 - avgProb1 + (1 - avgProb1).^2; % TODO % Optimization cost function

% Trial information
handles.curExpectedResults = []; % The current expected results matrix for one run. Each row, when combined with the same row from the interstimulation periods matrix, gives the required information for 1 sequence. Updated each time a related parameter gets changed
% TODO REMOVE % handles.curInterStimPeriods= []; % The current interstimulation periods matrix for one run. Each row, when combined with the same row from the expected results matrix, gives the required information for 1 sequence. Updated each time a related parameter gets changed
handles.curRunOrders = []; % The trial order for each run.  Each row represents one run, each column a trial in that run.  The number in that location represents the trial number (the row in the expected results/interstim period matrices) to do at that location in the trial.
handles.curTrialNum = 1; % The next trial to be started (iterates each time a trial is started, so the relevant data can be stored)

% Theta optimization information
% Prepare the estimation function for the behavior of the neurons
handles.FPDataFile = 'logisticFPTable.mat';
load(handles.FPDataFile); % Load in the pre-modeled FP data in the form of a matrix
handles.allAlphas = allAlphas;
handles.allSigmas = allSigmas;
handles.allDurs = allDurs;
handles.PAll = PAll;

% These values must be adjusted later, when new laser parameters are found
handles.allStrs_Raw = allStrs;
handles.b0_Raw = .125; % The beta value which was used to calculate the FP data in the FP data file

% % The FP data was constructed using a control voltage scale, so convert it
% % to a power scale
% irrAtMaxG = 1070*handles.powToIrr; % This was calculated using the testLaserParameters.m script in the TDT tests.  It represents the average laser power at 3V (the max of allStrs at time of coding) across the previous experiments.
% handles.irrRatio = irrAtMaxG/allStrs(end); % The ratio by which allStrs needs to be adjusted, and the inverse of that which b0 must be adjusted
% handles.allStrs = allStrs*handles.irrRatio;
% handles.b0 = b0/handles.irrRatio;

% Initialize the parameters needed to keep track of the control voltage
% given the stimulation strength
handles.fControlVFromIrr = @(irr, laserParameters) fzero(@(x) irr/handles.powToIrr - polyval(laserParameters, x), 2.5);
handles.fIrrFromControlV = @(controlV, laserParameters) handles.powToIrr*polyval(laserParameters, controlV);

handles = setLaserParameters(handles, [-9.9757   75.5625  201.5847  -10.7840]);
handles = setControlVoltages(handles);

% Set parameters for the adaptive optimizer
handles.numSequencesToFit = 40; % 20; % 10; % Changed from 10 to 20 on 1/29/2019, changed from 20 to 40 on 3/14/2019
handles.adaptiveBufferReplacementSize = 10; % Number of new sequences to add to the adaptive buffer each "block" (should divide evenly into sequencesPerRun
handles.gracePeriod = 5; % The number of sequences to wait before initially calculating
handles.lastCalculatedG = []; % A matrix to store the G pair that was computed during the last prepareAdaptive optimiziation
handles.lastCalculatedT = []; % A vector to store the T pair that was computed during the last prepareAdaptive optimiziation
handles.trialStartTime = 0; % The timestamp at which the trial starts

handles.doParallelStimAndCompute = false; % Should the GT pairs be updated asynchronously (i.e. continue to stimulate the brain while Matlab is calculating the next GT pair)?  Alternative: Update sequentially, so the stimulations stop while Matlab is computing the best GT pair

% Create warning sound object
handles.rawWarn = [sin(1:.3:400), sin(1:.46:400), sin(1:.4:400)];
handles.doPlayWarn = get(handles.playWarnCheck, 'Value');
handles.warnVolume = str2double(get(handles.warnVolumeEdit, 'String'))/100;
handles.warnAudio = audioplayer(handles.warnVolume*handles.rawWarn, 22050);

% Initialize the input varData's
% The first value represents the last valid value that the field had.  The
% second represents the lowest allowed value, and the third represents the
% highest allowed value
% Manual Stimulation
% handles.manStrengthVarData = [str2double(get(handles.manStimStrength,'String')) 0 handles.maxStrength];
% handles.manDurationVarData = [str2double(get(handles.manStimDuration,'String')) 0 handles.maxDuration];
handles.manDurationVarData = [str2double(get(handles.manStimDuration,'String')) 0 1000];
% handles.manStrength1VarData = [str2double(get(handles.controlStrength1,'String')) 0 handles.maxStrength];
handles.manDuration1VarData = [str2double(get(handles.controlDuration1,'String')) 0 handles.maxDuration];
% handles.manStrength2VarData = [str2double(get(handles.controlStrength2,'String')) 0 handles.maxStrength];
handles.manDuration2VarData = [str2double(get(handles.controlDuration2,'String')) 0 handles.maxDuration];

% yPlot Controller input
handles.yPlotVarData = [str2double(get(handles.controllerYPlotLevel,'String')) 0 Inf];

% Neuron 1
handles.durations1CurVals = [1 2 5 10 15];
handles.durations1VarDataBounds = [0 handles.maxDuration];

% Neuron 2
handles.durations2CurVals = [1 2 5 10 15];
handles.durations2VarDataBounds = [0 handles.maxDuration];

% Sequences
handles.numRunsVarData = [str2double(get(handles.numRuns, 'String')) 1 Inf];
handles.stimPeriodsVarData = [str2double(get(handles.stimPeriods, 'String')) 0 Inf];
handles.runsPerBlockVarData = [str2double(get(handles.adaptingRunsPerEpoch,'String')) 1 Inf];
% TODO REMOVE % handles.stimsPerSequenceVarData = [str2double(get(handles.stimsPerSequence, 'String')) 1 Inf];
handles.seqPeriodLowVarData = [str2double(get(handles.seqPeriodLow, 'String')) 0 Inf];
handles.seqPeriodHighVarData = [str2double(get(handles.seqPeriodHigh, 'String')) 0 Inf];
% handles.startStepFacVarData = [str2double(get(handles.startStepFac,'String')) .001 .999];

% Size of smoothing window
handles.smoothingSizeVarData = [str2double(get(handles.smoothingSize, 'String')) 2 Inf];

% Volume
handles.warnVolumeVarData = [str2double(get(handles.warnVolumeEdit, 'String')) 0 100];

% Testing
% handles.testNoiseVarData = [str2double(get(handles.testNoise, 'String')) -Inf Inf];

% Initialize neuron parameters (each row represents one neuron)
handles.curExpGT = [str2double(get(handles.controlStrength1, 'String')) str2double(get(handles.controlDuration1, 'String')); str2double(get(handles.controlStrength2, 'String')) str2double(get(handles.controlDuration2, 'String'))];
% handles.curExpGT = [handles.curControlVoltages(1) str2double(get(handles.controlDuration1, 'String')); handles.curControlVoltages(2) str2double(get(handles.controlDuration2, 'String'))];
handles.currentTriggerVoltage = nan(handles.numNeurons, 1); % The current NEW trigger voltage levels that are selected for each neuron (only one will be selected at a time).  These can be APPLIED to stimulations by selecting the "Apply to Selected" or "Apply to All" buttons
handles.currentTriggerTimespan = nan(handles.numNeurons, 2); % The current NEW trigger timespan that are selected for each neuron (only one will be selected at a time).  These can be APPLIED to stimulations by selecting the "Apply to Selected" or "Apply to All" buttons
handles.GUITriggerSettings = repmat(handles.allNeuronsStruct{1}.defaultTriggerSettings,2,1); % Voltage and timespan trigger settings for each electrode (if a neuron is cleared, the trigger settings are carried over)
handles.autoSaveFiles = []; % The file names of the save files that will be updated after new data is acquired
handles.stimsExist = false(handles.numNeurons, 1); % Are there any stimulations loaded into the GUI for this neuron?
handles.selectedStims = cell(handles.numNeurons, 1); % The indices in each the resultsData of each neuron that correspond to the currently selected stimulations/corresponding voltage traces (handles.allNeuronsStruct{X}.voltageTraces)
handles.characterizationDurations = cell(handles.numNeurons, 1); % The vector for each neuron that keeps track of which durations are being tested during the automated characterization phase
for i = 1:handles.numNeurons
    handles.characterizationDurations{i} = handles.(['durations' num2str(i) 'CurVals']);  % Use the durations initialized from the 'durationsXCurVals' variable
end

% Initialize data related to viewing sequence results
handles.sequenceTimeseriesFig = [];
handles.sequencePSpaceFig1 = [];
handles.sequencePSpaceFig2 = [];

% Set colors
handles.yellow = [1 1 0];
handles.red = [1 0 0];
handles.blue = [0 0 1];
handles.green = [0 .5 0];
handles.brightGreen = [0 1 0];
handles.gray = [.941 .941 .941];
handles.neuronPlotColor = {handles.blue, handles.green}; % Colors to be used in the plot routine's line-specifications that represent each neuron

% Set the mouse click callback for each of the 4 figures
set(handles.SDCurve1, 'ButtonDownFcn', @(object, eventData) figureClickCallback(object));
set(handles.SDCurve2, 'ButtonDownFcn', @(object, eventData) figureClickCallback(object));
set(handles.voltageTrace1, 'ButtonDownFcn', @(object, eventData) figureClickCallback(object));
set(handles.voltageTrace2, 'ButtonDownFcn', @(object, eventData) figureClickCallback(object));

% Set the keypress function of all graphics objects
set([handles.figure1; findall(handles.figure1, 'type', 'uicontrol');findall(handles.figure1, 'type', 'panel')], 'KeyPressFcn', @(hObject,eventdata) stimOnKeyPress(eventdata, guidata(hObject)));

% % Initialize the control sequence buttons
% updateControlSequenceGUI(handles);

% Initialize a trial expected results and inter-stim periods matrix
handles = updateTrialMatrices(handles);

% Display the name of the TDT block on screen for the user to copy/paste
curClock = clock;
tdtBlockString = sprintf('Acute_%d%.2d%.2d_xx', curClock(1), curClock(2), curClock(3));
tdtCharBlockString = sprintf('Characterization_%d%.2d%.2d_xx', curClock(1), curClock(2), curClock(3));
handles.tdtBlockString = tdtBlockString;
handles.tdtCharBlockString = tdtCharBlockString;
set(handles.tdtBlockName, 'String', tdtBlockString);
set(handles.tdtCharBlockName, 'String', tdtCharBlockString);

% Store the folder/file name of the controller file that needs to be
% re-written
% TODO SILICON %
handles.sourceController = {'PCSort_Seu', 'TetSort_Te1'}; % Name of the "source" (correct) controller in the csf file
handles.targetController = {'PCSort_Se2', 'TetSort_Te2'}; % Name of the "target" (to be modified) controller in the csf file
handles.targetLineStart = {'STATESTRSTART_', 'NAME='}; % Beginning of the line that is of interest in the csf file
handles.tetFirstLine = 'STATESTRSTARTDefTetrode'; % The first line in the tetrode file to copy.  The last line will be the line before the target controller's segment
handles.controllerFile = {'Two electrode SpikePac', 'Silicon Probe Record'};
handles.yPlotFile= 'Silicon Probe Explore 1';
handles.modCSFIndex = 2;
handles.controllerFilePath = {['C:\TDT\OpenEx\MyProjects\Acute Control PCA\UserFiles\' handles.controllerFile{1} '.csf'], ['C:\TDT\OpenEx\MyProjects\Acute Control PCA\UserFiles\' handles.controllerFile{2} '.csf']};
handles.yPlotFilePath = ['C:\TDT\OpenEx\MyProjects\Acute Control PCA\' handles.yPlotFile '.xpc'];
handles.haveControllerFilePath = [false false];
set(handles.controllerExistingFilePath, 'String', handles.controllerFilePath{handles.modCSFIndex});
set(handles.yPlotExistingFilePath, 'String', handles.yPlotFilePath);
handles.controllerFileExtension = '.csf';
handles.yPlotFileExtension = '.xpc';
handles.controllerFileModified = [false false]; % TODOCSF % Has the controller file been modified so that both SpikePac modules have identical eigenvectors and spheres?
handles = checkControllerFileStatus(handles);

% Display the name of the controller file on screen for the user to copy/paste
controllerFileString = {sprintf('%s_%d%.2d%.2d', handles.controllerFile{1}, curClock(1), curClock(2), curClock(3)) sprintf('%s_%d%.2d%.2d', handles.controllerFile{2}, curClock(1), curClock(2), curClock(3))};
handles.controllerSaveName = controllerFileString;
set(handles.controllerFileName, 'String', controllerFileString{handles.modCSFIndex});

% Set up default conditions for the trial to ensure the user does not run a
% trial with stupid settings
defaults.numRuns = 20;
defaults.stimPeriods = 100;
defaults.stimsPerSequence = 5;
defaults.seqPeriodLow = 1;
defaults.seqPeriodHigh = 3;
defaults.useAdaptive = true;
defaults.useAdaptingEpochs = false;
defaults.adaptingRunsPerEpoch = 2;
defaults.usePCACheck1 = true;
defaults.usePCACheck2 = true;
defaults.TDTConnectButton = {'TDTConnected', true};
defaults.optServeConnect = {'optServeConnected', true};
handles.defaults = defaults;

checkDefaults(handles);

% Set up the parallel pool if needed
if handles.useParallel && isempty(gcp)
   parpool('IdleTimeout', 5*60);
   handles.parConstants = []; % A structure which will later hold the parallel.pool.Contant values needed to perform the optimizations in parallel
end

% Set up the SD plot
handles = drawSDPlot(handles);

% Try connecting to TDT by connecting and then testing the connectivity
% (similar to TDTConnectButton callback
handles = startTDTComm(handles);
% handles = checkTDTConnectivity(hObject, handles);

% Get the optimization server connection ready
handles.optServeT = tcpip('128.197.49.98', 30000, 'NetworkRole', 'client');
handles.optServeT.OutputBufferSize = 2^18;
handles.optServeT.InputBufferSize = 2^18;

guidata(hObject, handles); % Update handles

function handles = setControlVoltages(handles)
% Get the control voltage given a desired power
curManPower = str2double(get(handles.manStimStrength, 'String'));
curPower1 = str2double(get(handles.controlStrength1, 'String'));
curPower2 = str2double(get(handles.controlStrength2, 'String'));

handles.curManControlVoltage = handles.fControlVFromIrr(curManPower, handles.curLaserParameters);
handles.curControlVoltages = [handles.fControlVFromIrr(curPower1, handles.curLaserParameters) handles.fControlVFromIrr(curPower2, handles.curLaserParameters)];

% Update the GUI
set(handles.controlV1, 'String', num2str(handles.curControlVoltages(1)));
set(handles.controlV2, 'String', num2str(handles.curControlVoltages(2)));
set(handles.controlV, 'String', num2str(handles.curManControlVoltage));


function handles = setLaserParameters(handles, varargin)
% This function will set the laser parameters in handles.  If another
% parameter is passed (assumed to be a 1x4 array), it will use the inputs
% as the new laser parameters.  Otherwise, it will get the parameters
% directly from TDT.

if nargin == 1
    % Get the parameters from TDT
    if handles.TDTConnected
        A = handles.DA.GetTargetVal('RZ2.LasEquA');
        B = handles.DA.GetTargetVal('RZ2.LasEquB');
        C = handles.DA.GetTargetVal('RZ2.LasEquC');
        D = handles.DA.GetTargetVal('RZ2.LasEquD');
        laserEq = [A B C D];
    else
        return;
    end
elseif nargin == 2
    % Use the parameter passed in as the new laser equation parameters
    laserEq = varargin{1};
end

% Save the new laser parameters
handles.curLaserParameters = laserEq;

% Recalculate the upper limit of power that can be used given these
% parameters
handles.maxIrradiance = handles.fIrrFromControlV(handles.maxStrength, laserEq);
handles.manStrengthVarData = [str2double(get(handles.manStimStrength,'String')) 0 handles.maxIrradiance];
handles.manStrength1VarData = [str2double(get(handles.controlStrength1,'String')) 0 handles.maxIrradiance];
handles.manStrength2VarData = [str2double(get(handles.controlStrength2,'String')) 0 handles.maxIrradiance];

handles.strength1MaxVarData = [str2double(get(handles.strength1Min,'String')) 0 handles.maxIrradiance];
handles.strength1MinVarData = [str2double(get(handles.strength1Max,'String')) 0 handles.maxIrradiance];
% handles.strength2MaxVarData = [str2double(get(handles.strength2Min,'String')) 0 handles.maxIrradiance];
% handles.strength2MinVarData = [str2double(get(handles.strength2Max,'String')) 0 handles.maxIrradiance];

% Scale the upper limit of the SD curves according to the new maximum power
set(handles.SDCurve1, 'Ylim', [0 handles.maxIrradiance]);
set(handles.SDCurve2, 'Ylim', [0 handles.maxIrradiance]);

% Adjust the FP lookup table values so that it is scaled to the current
% highest power
irrAtMaxG = handles.fIrrFromControlV(handles.allStrs_Raw(end), laserEq); % The laser irradiance at 3V (the max of allStrs at time of coding) across the previous experiments.
handles.irrRatio = irrAtMaxG/handles.allStrs_Raw(end); % The ratio by which allStrs needs to be adjusted, and the inverse of that which b0 must be adjusted
handles.allStrs = handles.allStrs_Raw*handles.irrRatio;
handles.b0 = handles.b0_Raw/handles.irrRatio;

% Write the laser parameters to the GUI
set(handles.laserParamA, 'String', num2str(laserEq(1)));
set(handles.laserParamB, 'String', num2str(laserEq(2)));
set(handles.laserParamC, 'String', num2str(laserEq(3)));
set(handles.laserParamD, 'String', num2str(laserEq(4)));
set(handles.maxControlV, 'String', num2str(handles.maxIrradiance, '%0.1f'));

function isDefault = checkDefaults(handles)
% Go through the defaults structure, and ensure that the default values
% are currently set
defaults = handles.defaults;
allUIs = fields(defaults);
isDefault = true; % Is the GUI currently set to default settings?

for fieldInd = 1:length(allUIs)
    getUIField = true;
    field = allUIs{fieldInd};
    val = defaults.(field);
    type = 'String';
    if isnumeric(val)
        type = 'String';
        fEq = @strcmp;
        fVal = @num2str;
    else
        fEq = @eq;
        fVal = @(x) x;
        if islogical(val)
            type = 'Value';
        elseif iscell(val)
            getUIField = false;
            compField = val{1};
            val = val{2};
        end
    end
    
    
    if (getUIField && fEq(get(handles.(field), type), fVal(val))) || (~getUIField && fEq(handles.(compField), fVal(val)))
        set(handles.(field), 'BackgroundColor', handles.gray);
    else
        set(handles.(field), 'BackgroundColor', handles.yellow);
        isDefault = false;
    end
end

function figureClickCallback(hObject)
% CHECKED
% This function waits for a double click on the selected window. When the
% window is double clicked, a new figure will be created that has an exact
% copy of the data on the figure, so that it may be explored easier

if strcmp(get(gcf,'selectiontype'), 'open')
    % Create a new figure with axes (or focus on the figure that was previously
    % spawned from this callback)
    figureNum = 1;
    hFigure = createFig(figureNum);
    hAxes = axes('Parent', hFigure); % Make sure that the new axis is the child of the new figure
    
    % Copy data from hObject over to the new figure
    copyobj(allchild(hObject), hAxes);
    
    set(hAxes,'Xlim', get(hObject, 'xlim'), 'YLim', get(hObject, 'ylim'));
end

% --- Executes on button press in openStructFromFile.
function openStructFromFile_Callback(hObject, eventdata, handles)
% CHECKED
% hObject    handle to openStructFromFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
openStructCallback(hObject, handles);

function openStructCallback(hObject, handles)
% CHECKED
% Extract needed information from the object that called this function
% Prepare the file selecting loop
readyToLoad = false; % Has the user chosen a file?
success = false; % Is the loading operation successful?

while ~readyToLoad
    % Allow the user to select a file
    [fName, fPath, fInd] = uigetfile('.mat', 'Select neuron file', handles.dataLoc);
    
    if fInd == 1
        % If the file is a .mat file, then continue with the loading
        fileToLoad = [fPath fName];
        
        % Record the path as the current data directory
        handles.dataLoc = fPath;
        
        readyToLoad = 1;
    else
        switch fInd
            case 0
                % If the dialog box was canceled, break out of this loop,
                % and do not attempt to load a file
                break;
            otherwise
                % If a file of the wrong type was loaded, allow the user to
                % try again.
                choice = questdlg('This file is not valid.  Try again?','Invalid Selection','Yes','No','Yes');
                
                if strcmp(choice, 'No')
                    % If the user doesn't want to try again, then break out
                    % of the loop. Just can't make up their mind, can they?
                    break;
                end
        end
    end
end

% If a good file name was selected
if readyToLoad
    % Load the file
    [handles, success] = loadStructFile(handles, fileToLoad);
end

% If loading the file was not successful, notify the user
if ~success && readyToLoad
    % If a proper file name was found, but it was still unsuccessful, then
    % something went wrong during loading
    errordlg(sprintf('Neuron file was not properly loaded. Check that the correct file was chosen, and that it contains the variable %s.', handles.structVarName),'Open File Error')
end

% Record changes to handles
guidata(hObject, handles)

function strength1Min_Callback(hObject, eventdata, handles)
% hObject    handle to strength1Min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strength1Min as text
%        str2double(get(hObject,'String')) returns contents of strength1Min as a double
handles = checkNewInput(hObject, handles, 'strength1MinVarData');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function strength1Min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strength1Min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function strength1Max_Callback(hObject, eventdata, handles)
% hObject    handle to strength1Max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strength1Max as text
%        str2double(get(hObject,'String')) returns contents of strength1Max as a double
handles = checkNewInput(hObject, handles, 'strength1MaxVarData');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function strength1Max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strength1Max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function durations1_Callback(hObject, eventdata, handles)
% hObject    handle to durations1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of durations1 as text
%        str2double(get(hObject,'String')) returns contents of durations1 as a double

handles = checkNewVectorInput(hObject, handles, 'durations2CurVals', 'durations2VarDataBounds');

% Record the new values
neuronNum = get(hObject, 'UserData');
handles.characterizationDurations{neuronNum} = handles.durations2CurVals;

guidata(hObject, handles);


function string = vec2String(vector)
% This function converts an input vector of numbers into a formatted string
% to display on the GUI
numString = []; % The string that will contain all of the numbers in order, separated by spaces
for i = 1:length(vector)
    numString = [numString sprintf('%d ', vector(i))]; % Add each number followed by a space
end
numString(end) = []; % Remove the last space from the number string
string = ['[' numString ']']; % Put bracket around the number string to finish the formatting

% --- Executes during object creation, after setting all properties.
function durations1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to durations1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [handles, success] = loadStructFile(handles, fileToLoad)
% CHECKED
% This function will load a new neuronstruct, designated by file, into the
% program

success = true; % Was the file properly loaded?

% Check if file exists
if ~exist(fileToLoad, 'file')
    %     disp('File not found');
    return;
end

% Attempt to extract the allNeuronStructs variable from the file
try
    fileToLoadObject = matfile(fileToLoad);
    newObject = fileToLoadObject.(handles.structVarName); % Attempt to load the allNeuronStructs varaible into the local workspace
catch e
    disp(e); % Display the error
    return; % Return without setting the structsLoaded variable to true
end

% Determine what type of information is contained in the file, and check
% that the information is properly formatted
if iscell(newObject)
    % The object is a cell, so it is a new save file, where each cell in
    % the array represents a struct.  Load it into both slots in the
    % program.
    newNeuronStruct = cell(1, length(newObject));
    for neuronNum = 1:length(newObject)
        % Make sure that the structure is a proper neuron structure (and bring it
        % up to compatability with the latest version number
        newNeuronStruct{neuronNum} = makeProperNeuronStruct(newObject{neuronNum});
        
        if ~isNeuronStruct(newNeuronStruct{neuronNum})
            success = false;
            return; % If the variable that is loaded is not in fact a neuron structure, return from this function without loading it
        end
    end
    
    % If both structs are ok, then load them in
    handles = loadNeuronStruct(handles, newNeuronStruct); % Load the neuronStruct into the GUI
    
elseif isstruct(newObject)
    % The object is a struct, so it is probably a single neuron's structure
    
    % Make sure that the structure is a proper neuron structure (and bring it
    % up to compatability with the latest version number
    newNeuronStruct = makeProperNeuronStruct(newObject);
    
    if ~isNeuronStruct(newNeuronStruct)
        success = false;
        return; % If the variable that is loaded is not in fact a neuron structure, return from this function without loading it
    end
    
    % If only one object is being loaded, ask the user what should be done
    % with it
    c1 = 'Neuron 1 (Blue)'; % Choice 1 (neuron 1)
    c2 = 'Neuron 2 (Green)'; % Choice 2 (neuron 2)
    choice = questdlg('Which slot should this be loaded into?', 'Neuron Slot Choice', c1, c2, c1);
    if strcmp(choice, c1)
        neuronNum = 1;
    elseif strcmp(choice, c2)
        neuronNum = 2;
    else
        success = false;
        return;
    end
    
    newNeuronStruct = handles.allNeuronsStruct;
    newNeuronStruct{neuronNum} = newObject;
    
    % Load the structure
    handles = loadNeuronStruct(handles, newNeuronStruct);
end

function isProper = isNeuronStruct(allegedStruct)
% CHECKED
% This function will return true if a neuron struct is proper, either from
% the most recent version, or of a previous version, and will simply add on
% initialized versions of the missing fields

isProper = ~isempty(makeProperNeuronStruct(allegedStruct)); % If the makeProperNeuronStruct returns an empty variable, then the struct was not a proper neuronStruct
return;

function newStruct = makeProperNeuronStruct(oldStruct)
% CHECKED
% If oldStruct is formatted like an neuronStruct, this function will make
% sure that it has all of the fields necessary to be called a neuronStruct
% in the lastest version (some of the fields may be improperly initialized,
% but if the struct has all of the required fields, it is safe to assume
% that it was created by 'initializeAllNeuronStructs', and is therefore
% properly initialized). It will return the new structure.

newStruct = [];

% Is the variable a struct?
if (~isstruct(oldStruct))
    return;
end

% See if the alleged struct contains all of the fields that would be expected from a properly initialized
% neuronStruct, and ONLY those fields
if isequal(sort(fieldnames(oldStruct)), sort(fieldnames(initializeNeuronStruct)))
    newStruct = oldStruct; % The old structure is correct, and can be used
    return;
end

% If it is a struct, but the fieldnames do not match perfectly, try
% modifying the struct so that it has initialized fieldnames from more
% recent versions.  Versions are tested IN REVERSE ORDER so as to not
% overwrite data
highestVer = 2; % Highest supported version at the moment

for ver = highestVer:-1:2
    switch ver
        case 2
            % Add version 2 fields
            oldStruct.useIrradiance = false; % The neuron struct being modified probably does not use irradiance control
    end
    
    if isequal(sort(fieldnames(oldStruct)), sort(fieldnames(initializeNeuronStruct)))
        newStruct = oldStruct; % The newly modified structure is correct, and can be used
        return
    end
end

% If it still has not passed the tests, then it is not proper, and should
% not be loaded
return;

function handles = characterizeNeuronCallback(handles)
% This function is the callback to the "Characterize Neuron" buttons for
% each neuron.  This goes through the process of testing each laser pulse
% duration requested by the user, and determining what laser strength is
% required to make the neuron spike 50% of the time at each duration.  This
% strength duration curve will then be used to control the two neurons
% individually of each other.

%% Should the characterization only consider data from this time after, or should is consider all data?
ignorePreviousResults = true;

% Delete previous estimates
% handles.allNeuronsStruct{neuronNum}.characterizedDurs = [];

% Get time of characterization
if ignorePreviousResults
    charTime = etime(clock, handles.startTime); % The time at which the characterization was started.  If ignoring previous results, the SD Curve will only be built with data acquired after this time.
else
    charTime = 0;
end

% Set up the parallel pool if needed (this is done now because it's going
% to take a minute, so it should be done BEFORE we annihiliate the neuron but within thirty minutes of the actual trial)
if handles.useParallel && isempty(gcp)
   parpool('IdleTimeout', 5*60);
end

% Build the Strength-Duration curve
[handles, characterizationCanceled] = performStimulationBattery(handles, ignorePreviousResults);

% After the Strength-Duration curve has been built (if it has not been
% canceled), calculate the parameters of the neuron
if ~characterizationCanceled
    handles = fitSDCurveCallback(handles, charTime);
end

function [handles, canceled] = performStimulationBattery(handles, varargin)
% This function is used to test a real neuron to determine what strength(s)
% are required to make the neuron spike with a probability determined by
% targets.  It will do this via a modified bisection algorithm with an
% x-tolerance-based end condition (TODO use logistic regression somehow?
% Likelihood?).  After the end condition has been met, all of the gathered
% data up to this point will be fit to a sigmoid curve to smooth the data,
% at which point the target values will be found.

% Settings for the battery
useSigmoid = false;

% Parse inputs
ignorePreviousResults = false; % Should the algorithm ignore results that were found before starting the characterization? (Will probably replace with learning later)
if nargin > 2
    ignorePreviousResults = varargin{1};
end

% Set up parameters
characterizationStartInd = zeros(handles.numNeurons, 1);
durations = [];
for neuronNum = 1:handles.numNeurons
    characterizationStartInd(neuronNum) = handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc; % The current location in resultsData (used to set all of the time points to be the same, because it is assumed that during the characterization, the system does not change)
    durations = [durations handles.characterizationDurations{neuronNum}'];
end
uniqueDurs = unique(durations); % Get the durations that will be tested to build the curve

canceled = false;
setCharacterizationTimeSame = true; % Should the program set the times of all stimulations during the neuron characterization to be the same, so that they all get the same weight during the SD curve fitting?
targets = .5; % Get the probability value that the strength-duration curve represents
minStrength = 0; % The lowest strength tested will always be 0.

lowMaxStrength = str2double(get(handles.strength1Min, 'String'));
highMaxStrength = str2double(get(handles.strength1Max, 'String'));
% minStr2 = str2double(get(handles.strength2Min, 'String'));
% maxStr2 = str2double(get(handles.strength2Max, 'String'));

% lowMaxStrength = min([minStr1 minStr2]);
% highMaxStrength = max([maxStr1 maxStr2]);

% Program control
% Additional program control: 'ignorePreviousResults' above (see comment
% for more information)
maxNumIterations = 30; % Maximum number of iterations to go through for the probability finding loop (simplest end point)
minNumIterations = 15; % Minimum number of iterations to go through for the probability finding loop (simplest continue signal)
% useSmoothedFinalValue = true; % Should the final value that the program uses be determined by using a sigmoid curve to smooth all of the data
flankWidth = 0; % If set to 0, no flanks will be used.  If set to some other positive value (preferably less than .5), then additional NON-TARGET strengths will be tested on either side of the next target strength each iteration, to increase stability
granularity = .01; % This value determines how finely the durations can change (so, for a granularity of .1, a calculated x value of .46 will round to .5, which will then be tested)
xThreshRatio = .05; % The percentage of the maxStrength that will be the xThreshold (if maxStrength = 5 and xThreshRatio = .1, then the xThresh will be .5)

% Prepare algorithm parameters
% Prepare parameters for all durations
numTargets = size(targets,2); % Number of targets
useFlanks = flankWidth > 0; % Will flanks be used in the algorithm?
numNewStimsPerIteration = 1 + 2*useFlanks; % Number of new values each iteration (1 for the target, and 2 more for the flanks on either side of it)
minDur = min(uniqueDurs);
maxDur = max(uniqueDurs);

% Prepare sigmoid problem parameters
% Define fitting options (modified 5/15/19)
s=fitoptions('Method','NonlinearLeastSquares','Startpoint',[0 0],'Lower',[0 -Inf]);
sigmoidFit=fittype('(1/2)*(1+tanh(a*x+b))','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b'},'options',s);
sigmoidFun = @(x, params) (1/2).*(1 + tanh(params(1).*x + params(2))); % - .25)*1.5; % This function produces a pure sigmoid curve, which is stretched beyond 0 and 1, so that low and high values will either NEVER or ALWAYS return a 1, respectively.
% s=fitoptions('Method','NonlinearLeastSquares','Startpoint',[1 0 0],'Lower',[0 0 -Inf]);
% sigmoidFit=fittype('a/2*(1+tanh(b*x+c))','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c'},'options',s);
% sigmoidFun = @(x, params) (params(1)/2).*(1 + tanh(params(2).*x + params(3))); % - .25)*1.5; % This function produces a pure sigmoid curve, which is stretched beyond 0 and 1, so that low and high values will either NEVER or ALWAYS return a 1, respectively.


% % Check all durations, make sure that a duration that has already been
% % added to characterizedDurs will not be tested here
% remDurInds = ismember(durations, handles.allNeuronsStruct{neuronNum}.characterizedDurs); % This function returns the indexes in durations which are already represented in characterizedDurs (and should therefore be removed from durations)
% durations(remDurInds) = []; % Remove durations that have already been added to characterizedDurs (and have therefore already been characterized)

numDurations = size(uniqueDurs,1); % Get the number of durations that will be tested

% Go through each duration that will be tested
h = waitbar(0, 'Starting...', 'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0);

% Make sure that the waitbar is always set to be on top:
frames = java.awt.Frame.getFrames();
frames(end).setAlwaysOnTop(1);

for durInd = 1:numDurations
    % Get the current duration
    curDuration = uniqueDurs(durInd); % Get the current duration that will be tested
    unitsToTest = any(durations == curDuration, 1); % The units that should be tested with this duration
    numUnitsToTest = sum(unitsToTest); % The total number of units that should be tested with this duration
    
    % Get all possible strength values, taking into account the max
    % strength for this duration and the granularity
    ratio = (curDuration - minDur)/(maxDur - minDur); % The relative value of this duration between the minimum and the maxiumum duration
    maxStrength = lowMaxStrength + (1 - ratio)*(highMaxStrength - lowMaxStrength); % The relative value of the max strength to be used between the minimum and the maxiumum strength, which is inversely proportional to the magnitude of the duration (this allows low durations, such as 1ms, to have a high max strength, while long durations, such as 15ms, can have a low max-strength without overstimulating the cells)
    possibleXs = 0:granularity:maxStrength; % These are all of the possible values of X to be tested
    numPossibleXs = size(possibleXs, 2); % Total number of possible x's
    xThresh = max(xThreshRatio*maxStrength, granularity); %.15; % If the predicted x value changes by less than this amount, then the optimization will be completed
    
    % Update the waitbar
    waitbar(durInd/numDurations,h,sprintf('Testing %dms (%d out of %d durations)', curDuration, durInd, numDurations));
    
    % Start measurement by confirming the boundaries (should this be done?)
    % Go through the min and max strengths
    % Test Min
    % Check if the characterization is being canceled
    if getappdata(h,'canceling')
        canceled = true;
        break;
    end
    
    %% Test min and max
    % Test min
    waitBeforeStim(handles); % Make sure that some random amount of time has passed since the last stimulation
    [handles, minStrengthResults] = stimulate(handles, minStrength, curDuration);
    
    % Check if the characterization is being canceled
    if getappdata(h,'canceling')
        canceled = true;
        break;
    end
    
    % Test Max
    waitBeforeStim(handles); % Make sure that some random amount of time has passed since the last stimulation
    [handles, maxStrengthResults] = stimulate(handles, maxStrength, curDuration);
    
    usedMinStrengthResults = minStrengthResults(unitsToTest);
    usedMaxStrengthResults = maxStrengthResults(unitsToTest);
    
    usedMinStrengthResults = usedMinStrengthResults(:)';
    usedMaxStrengthResults = usedMaxStrengthResults(:)';
    
    % check boundaries
    minProducedSpike = all(usedMinStrengthResults == 1); % Determine if the minimum strength produced no spikes
    maxProducedSpike = all(usedMaxStrengthResults == 1); % Determine if the maximum strength produced a spike
    
    % Check if the boundaries are satisfactory
    if (minProducedSpike || ~maxProducedSpike)
        % Display options to user
        s = '';
        if minProducedSpike
            s = 'The lower strength fired an action potential.  ';
        end
        if ~maxProducedSpike
            s = [s 'The upper strength did not fire an action potential.  '];
        end
        
        % See if user wants to continue, or try again later
        playWarn(handles);
        choice = questdlg([s 'Continue with characterization? ("Yes" to ignore, "No" to adjust boundaries)'], 'Boundary error: The boundaries may need to be expanded.', 'Yes', 'No', 'Yes');
        
        switch choice
            case 'No'
                canceled = true;
                break;
        end
    end
    
    initialResults = [usedMinStrengthResults;usedMaxStrengthResults];
    initialStimStrengths = [minStrength;maxStrength];
    
    % If using a sigmoid, a third initial stimulation is needed
    if useSigmoid
        % Stimulate at the mid-point between the two strengths
        midStrength = mean([maxStrength minStrength]);
        waitBeforeStim(handles);
        [handles, midStrengthResults] = stimulate(handles, midStrength, curDuration);
        
        usedMidStrengthResults = midStrengthResults(unitsToTest);
        usedMidStrengthResults = usedMidStrengthResults(:)';
        
        initialResults = [initialResults;usedMidStrengthResults];
        initialStimStrengths = [initialStimStrengths;midStrength];
    end
    
    % Initialize results and strength tracking variables
    if ignorePreviousResults
        numInitialResults = size(initialResults, 1);
        curResultsAllCell = cell(maxNumIterations, numUnitsToTest);
        curStrengthsAllCell = cell(maxNumIterations, numUnitsToTest);
        curResultsAllCell(1:numInitialResults) = mat2cell(initialResults, 2, ones(numUnitsToTest));
        curStrengthsAllCell(1:numInitialResults) = mat2cell(repmat(initialStimStrengths, 1, numUnitsToTest), 2, ones(numUnitsToTest));
        curResultsLoc = 3; % Location of the next value in curResults and curStrengths
    end
    
    %% Start search for target values
    % Prepare memory for variables needed during the algorithm
    finalVals = zeros(size(targets)); % Vector to store the final values outputed by the algorithm
    numRequiredIterations = zeros(size(targets)); % Vector to store the number of required iterations to acheive the target (or to give up)
    lastStrengths = zeros(numUnitsToTest, numTargets*(1 + 2*useFlanks));
    thisDurUnits = find(unitsToTest); % The neurons being tested in this duration
    done = false(numUnitsToTest, numTargets); % Done optimizing yet? Indexed according to the neurons that are actually being characterized this round
    curIteration = 1; % Algorithm iteration counter
    unitToTestInd = 1;
    unitsContinueTesting = any(~done', 1); % A boolean matrix in the same form of unitsToTest, keeps track of whch units need to continue being tested
    while any(~done(:))
        unitToTestInd = mod(unitToTestInd, sum(unitsContinueTesting)) + 1; % Find the index in THE TRUE (1) VALUES OF unitsContinueTesting that should be tested this round (e.g. unitsContinueTesting == [0 1 0 1], this will go between 1 and 2, to refer to units in thisDurUnits 2 and 4) ("+1" to alternate which unit it being tested)
        potentialUnitsToTestThisLoop = find(unitsContinueTesting); % Get the indices in unitsToTset of all units that are up for testing this loop
        thisDurUnitsInd = potentialUnitsToTestThisLoop(unitToTestInd); % Get the index in thisDurUnits that will be tested
        neuronNum = thisDurUnits(thisDurUnitsInd); % Get the neuron number that this refers to
        thisDone = done(thisDurUnitsInd, :);
        
        % Update the GUI (Let the user know which unit is being tested)
        set(handles.characterizationNeuron, 'String', num2str(neuronNum));
        
        if ignorePreviousResults
            % Prepare the previous data
            curResultsAll = curResultsAllCell{thisDurUnitsInd};
            curStrengthsAll = curStrengthsAllCell{thisDurUnitsInd};
        end
        
        % Prepare for this iteration
        trialTargetsDone = reshape(repmat(thisDone, numNewStimsPerIteration, 1), [numTargets*numNewStimsPerIteration, 1]); % A matrix to show if each target is done, reshaped to include the status of the flanks around each target, as well.
        numNewVals = sum(~trialTargetsDone); % The number of values being added (how many targets AND flanks are not yet done)
        
        % Calculate what the next x value to be tested will be
        if any(~thisDone)
            % Set up flanks if being used
            if ~useFlanks
                thisTrialTargets = targets; % thisTrialTargets represent the values that will be targeted this iteration, which includes all targets as well as the flanks on either side of them
            else
                thisTrialTargetsMatrix = [targets - flankWidth;targets;targets + flankWidth]; % Find all values that will be targeted in this trial, which includes all targets as well as flanks on either side
                thisTrialTargets = thisTrialTargetsMatrix(:)'; % Turn the matrix into a row vector, in the correct order
            end
            
            % Calculate corresponding strengths to test
            
            % First create the x and y data
            if ignorePreviousResults
                currentDataInds = ~isnan(curResultsAll); % Get the indices in curResults to use
                curResults = curResultsAll(currentDataInds); % Get only the results that are ready
                curStrengths = curStrengthsAll(currentDataInds); % Get only the strengths that are ready
            else
                % If previous results are going to be used, then get the
                % data straight from resultsData
                % Use results that are 1) associated only with this duration,
                % 2) not NaN's, and 3) gathered after this characterization
                % began
                allCharInds = characterizationStartInd(neuronNum):(handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc - 1);
                currentDataInds = ~isnan(handles.allNeuronsStruct{neuronNum}.resultsData(allCharInds,1)) & (handles.allNeuronsStruct{neuronNum}.resultsData(allCharInds,2) == curDuration); % Indices for all resultsData rows that should be used for this fit (have the correct duration, and are not NaN's)
                currentDataInds = [false(characterizationStartInd(neuronNum) - 1, 1); currentDataInds];
                curResults = handles.allNeuronsStruct{neuronNum}.resultsData(currentDataInds, 3); % Results associated with the used strengths
                curStrengths = handles.allNeuronsStruct{neuronNum}.resultsData(currentDataInds, 1); % Strengths associated with this duration
                
                allResults = zeros(size(curResults, 1), handles.numNeurons);
                for thisNN = 1:handles.numNeurons
                   allResults(:, thisNN) = handles.allNeuronsStruct{thisNN}.resultsData(currentDataInds, 3); 
                end
            end
            
            % Replace the min/max results with 0's and 1's (always assume
            % that the boundaries are correct, because the user already
            % confirmed wanting to use them)
            curResults(curStrengths == minStrength) = false; % Set lower boundary results to 0
            curResults(curStrengths == maxStrength) = true; % Set upper boundary results to 1
            
            thisTarget = thisTrialTargets(~trialTargetsDone); % The target SD curve isoline to find (usually 50% firing rate)
            if useSigmoid
                % Perform fitting calculation on data
                params = coeffvalues(fit(curStrengths, curResults, sigmoidFit));
                
                if useFlanks
                    thisNumTargets = length(thisTarget);
                   calculatedXVals = zeros(1, thisNumTargets); 
                   for targetInd = 1:thisNumTargets
                      thisTargetScalar = thisTarget(targetInd);
                      sigmoidZeroFindingFunction = @(x) sigmoidFun(x, params) - thisTargetScalar; % A function to use to find the next strength using a sigmoid
                      calculatedXVals(targetInd) = fzero(sigmoidZeroFindingFunction, 0);
                   end
                else
                    sigmoidZeroFindingFunction = @(x) sigmoidFun(x, params) - thisTarget; % A function to use to find the next strength using a sigmoid
                    calculatedXVals = fzero(sigmoidZeroFindingFunction, 0);
                end
            else
                % Calculate a weighted linear regression for the known data
                numCurrentData = sum(currentDataInds); % Number of data points to be interpolated now
                
                W = eye(numCurrentData); % Weights
                X = [curStrengths ones(numCurrentData,1)]; % X Values (strengths)
                Y = curResults; % Y Values (results)
                b = (X'*W*X)\(X'*W*Y); % Calculate the slope and intercept
                
                % Calculate the next most-likely position for the target, only for the
                % ones that are not done yet
                calculatedXVals = (thisTarget - b(2))/b(1);
                
            end
            % Ensure that this value is within the searching boundaries
            nextXTrue = max(min(calculatedXVals, maxStrength), minStrength);
            
            % Next, "snap" it to the nearest possibleXs value
            [~, possibleXInds] = min(abs(repmat(possibleXs,numNewVals,1) - repmat(nextXTrue', 1, numPossibleXs)), [], 2); % Get the index in possibleXs that corresponds to the nextX for each target
            nextStrengths = possibleXs(possibleXInds); % Get the rounded values of the next x's
            
            % Check finishing conditions
            % Store which targets were optimized before checking ending
            % conditions (to determine which are newly finished, so the final
            % value can be recorded
            previousDone = thisDone;
            
            % Check how close the function is to finding the target, by
            % checking the next X-value to be tested (if it is very close
            % to the last value, then no new stimulation needs to be made,
            % because the ending condition has already been satisfied).
            if (curIteration > 1)
                % Check the delta x threshold
                if (any(~thisDone) && (xThresh ~= 0))
                    lastStrengthsThisUnit = lastStrengths(unitToTestInd, :);
                    reachedThreshold = reshape(abs(nextStrengths(~thisDone) - lastStrengthsThisUnit(~thisDone)) < xThresh, [1, sum(~thisDone)]); % Is the newest strength less than xThresh away from the last one?
                    thisDone(~thisDone) = (curIteration > maxNumIterations) | ((minNumIterations <= curIteration) & reachedThreshold); % Is the current iteration count satisfying the minNumIterations count, the maxNumIterations count, and the reachedThreshold boundary?
                end
            end
            
            % Find which targets just finished optimizing (usually only 1
            % target, but code exists for multiple targets to be tested
            % simultaneously).
            justFinished = logical(previousDone ~= thisDone);
            
            % Store this neuron's results
            done(thisDurUnitsInd, :) = thisDone;
            if ignorePreviousResults
                curResultsAllCell{thisDurUnitsInd} = curResultsAll;
                curStrengthsAllCell{thisDurUnitsInd} = curStrengthsAll;
            end
            
            % Check if this unit is done
            unitsContinueTesting = any(~done', 1);
            
            % Measure using the next X values
            % Go through each strength to be measured, and sample the system
            newVals = nan(length(nextStrengths), 1); % Vector that contains results from these stimulation
            for strInd = 1:length(nextStrengths)
                % Check if the characterization is being canceled
                if getappdata(h,'canceling')
                    waitbar(durInd/numDurations,h, 'Canceling...');
                    canceled = true;
                    break; % Break out of strength-testing loop
                end
                
                % Get sample from neuron
                waitBeforeStim(handles); % Make sure that some random amount of time has passed since the last stimulation
                [handles, thisStimResults] = stimulate(handles, nextStrengths(strInd), curDuration);
                
                % Store results
                newVals(strInd) = thisStimResults(neuronNum); % Store in newVals, used to get the real final value from the characterization
                
                if ignorePreviousResults
                    curResultsAll(curResultsLoc) = thisStimResults(neuronNum); % Store in curResults, used to keep track of the results data within this characterization
                    curStrengthsAll(curResultsLoc) = nextStrengths(strInd); % Store the strength used to get the most recent result
                    curResultsLoc = curResultsLoc + 1; % Iterate the index of curResults
                end
            end
            
            if getappdata(h,'canceling')
                waitbar(durInd/numDurations,h, 'Canceling...');
                canceled = true;
                break; % Break out of while(~done) loop
            end
            newVals = newVals(1:numNewStimsPerIteration:end); % Use only the results from the target stimulations, not the flanks
            %             if any(logical(newVals))
            %                 durFiredSpike = true;
            %             end
        end
        
        if getappdata(h,'canceling')
            waitbar(durInd/numDurations,h, 'Canceling...');
            canceled = true;
            break; % Break out of multiple durations loop
        end
        
        % Clean up
        % If any of the targets have finished, store their values in the
        % finalVals vector
        finalVals(justFinished) = newVals(justFinished);
        numRequiredIterations(justFinished) = curIteration;
        
        % Remember the strength values that were tested this iteration, to
        % test the next iteration's end condition
        lastStrengths(unitToTestInd, :) = nextStrengths;
        
        % Iterate
        curIteration = curIteration + 1; % Iterate once for this loop
        %         neuronNum = neuronNum + 1; % Change which neuron is being tested
    end
    
    % Clear the GUI
    set(handles.characterizationNeuron, 'String', '');
    
    if getappdata(h,'canceling')
        waitbar(durInd/numDurations,h, 'Canceling...');
        canceled = true;
        break; % Break out of multiple durations loop
    end
    
    %     % Check that this duration actually caused a spike
    %     if durFiredSpike
    %         % Check that this duration isn't already in characterizedDurs
    %         if all(~ismember(handles.allNeuronsStruct{neuronNum}.characterizedDurs, curDuration))
    %             % Add this duration to the characterizedDurs vector
    %             handles.allNeuronsStruct{neuronNum}.characterizedDurs = [handles.allNeuronsStruct{neuronNum}.characterizedDurs;curDuration]; % Add this duration
    %         end
    %     else
    %         sprintf('Duration %d not characterized for neuron %d, because no action potential fired during characterization.\n', curDuration, neuronNum);
    %     end
end

% Ensure that, if the user wants it, all times for this characterization
% are the same (no difference in weighting between them)
if setCharacterizationTimeSame
    % If the user wants all weights from this characterization to be the
    % same, use the current time
    for neuronNum = 1:handles.numNeurons
        handles.allNeuronsStruct{neuronNum}.resultsData(characterizationStartInd(neuronNum):handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc - 1, 4) = etime(clock, handles.startTime);
    end
end

% Auto save the data (characterized durations)
handles = autosave(handles);

% Reset the app data to no longer canceling
setappdata(gcbf,'canceling',0);
delete(h);

function waitBeforeStim(handles)
% Before triggering the stimulation, wait until enough time has passed
% since the last stimulation

% Calculate important times
tTotal = etime(handles.nextStimTime, handles.lastStimTime); % Total amount of time between the pulses this time around
tRemain = etime(handles.nextStimTime, clock); % Total time remaining

if ~get(handles.testing, 'Value')
    % If the GUI is not in testing mode, do not worry about
    % waiting between each pulses
    % Prepare for giving the user graphical feedback
    if (tRemain > 0)
        % Make the slider visible
        set(handles.stimStatusSlider, 'Value', 0);
        set(handles.stimStatusSlider, 'Visible', 'on');
        
        % Wait for the remaining time to go down to zero
        while (tRemain > 0)
            % While there is still positive time until the next stim time,
            % update the GUI to show the user how much time remains
            set(handles.stimStatusSlider, 'Value', (tTotal - tRemain)/tTotal);
            set(handles.stimStatusText, 'String', sprintf('%0.1f', tRemain));
            
            drawnow;
            pause(.1);
            tRemain = etime(handles.nextStimTime,clock); % Total time remaining
        end
        
        % Reset GUI
        set(handles.stimStatusSlider, 'Visible', 'off');
        set(handles.stimStatusText, 'String', '');
    end
end



% --- Executes on button press in setTriggerVoltage1.
function setTriggerVoltage1_Callback(hObject, eventdata, handles)
% hObject    handle to setTriggerVoltage1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = setTriggerVoltage(hObject, handles); % Call the routine to set the trigger voltage, for neuron 1
guidata(hObject, handles); % Update handles

function handles = setTriggerVoltage(hObject, handles)
neuronNum = get(hObject, 'UserData');

if handles.stimsExist(neuronNum)
    neuronNumStr = num2str(neuronNum); % Neuron number as a string
    axes(handles.(['voltageTrace' neuronNumStr])); % Set the voltage trace axis as the current axis
    [~, newVoltage, button] = ginput(1); % Get the new trigger voltage, and set it to the trigger voltage of all selected spikes
    
    if button ~= 1
        % If the user didn't click the mouse, and instead hit a key, assume
        % that they want to cancel the selection
        return;
    end
    handles.currentTriggerVoltage(neuronNum) = newVoltage;
    handles = drawVoltageTraces(handles, neuronNum); % Draw the plot
end

% Update handles
guidata(hObject, handles);


function handles = drawVoltageTraces(handles, neuronNum)
% This function draws the currently selected (single or multiple)
% stimulations, selectable via the neuronXStimsListbox GUI interface. It
% will draw these traces, as well as the trigger voltage and timespan, and
% whether or not a spike was triggered.

neuronNumStr = num2str(neuronNum); % Make a string from the neuron's number
doAlignSpikes = get(handles.(['alignSpikes' num2str(neuronNum)]), 'Value'); % Does the user want the spikes to be aligned to the same position?
doPlotRaw = get(handles.(sprintf('viewRaw%d', neuronNum)), 'Value'); % Should the raw trace be plotted (as opposed to the filtered one)?
doPCA = get(handles.(sprintf('usePCACheck%d', neuronNum)), 'Value'); % Should the spikes be detected using PCA (or using a voltage threshold)
maxPCACode = 5; % The maximum pca code number (will color up to this number of codes differently on the plot, and everything after will be the same color)
axisHandle = handles.(['voltageTrace' neuronNumStr]);
doFullTrace = logical(get(handles.(['fullTrace' neuronNumStr]), 'value'));
axes(axisHandle); % Select the current neuron's voltage trace axis
cla(axisHandle); % Clear the axis

% Make a color map for displaying pca codes
pcaCMap = parula(maxPCACode - 1);

% Check if any stimulations are selected for this neuron
if isempty(handles.allNeuronsStruct) || isempty(handles.selectedStims{neuronNum})
    % The selected spikes matrix is empty (or the allNeuronStructs variable
    % is empty), clear the axis, and update the stimFlag for this neuron
    set(handles.(['stimFlag' neuronNumStr]), 'BackgroundColor', handles.gray, 'String', 'No Stim');
    
    return % End draw function
end

% Prepare for plotting the traces
numPlottedTraces = length(handles.selectedStims{neuronNum}); % Number of traces to plot
cMap = hsv(numPlottedTraces); % Colormap to use to color the traces
thisCollectionResults = false(numPlottedTraces,1); % The spike results from this collection of stimulation
minTime = 0; % The minimum value of the time axis

if doFullTrace
    maxTime = handles.traceTime; % The maximum value of the time axis, or the right-most time
else
    maxTime = 20; % The maximum value of the time axis, or the right-most time
end

% Determine what the y limits (minimum/maximum voltage) should be
maxVoltage = handles.currentTriggerVoltage(neuronNum); % The maximum voltage value of all of the traces.  Will be used to scale the axis later. Initialized to current trigger voltage of this neuron.
minVoltage = handles.currentTriggerVoltage(neuronNum);
for curStimInd = 1:numPlottedTraces
    traceInd = handles.selectedStims{neuronNum}(curStimInd); % Index in voltageTraces
    if doPlotRaw
        % If plotting the raw trace
        trace = handles.allNeuronsStruct{neuronNum}.voltageTracesRaw{traceInd}; % Trace to plot
    else
        % If plotting the filtered trace
        trace = handles.allNeuronsStruct{neuronNum}.voltageTraces{traceInd}; % Trace to plot
    end
    times = linspace(0, handles.traceTime, length(trace));
    selectedTimesMask = minTime < times & times < maxTime; % Only use the plotted parts of the trace in the calculation of min/max
    
    voltageThresh = handles.allNeuronsStruct{neuronNum}.triggerVoltage(traceInd); % Spike threshold for this trace
    % Find max, and compare it to running maximum
    maxVoltage = nanmax([maxVoltage max(trace(selectedTimesMask)) voltageThresh]); % Find the max between the running max, the maximum value of the trace, and the threshold voltage
    minVoltage = nanmin([minVoltage min(trace(selectedTimesMask)) voltageThresh]); % Find the min between the running min, the minimum value of the trace, and the threshold voltage
end
maxY = maxVoltage*1.1; % Give a small amount of headroom on the top of the plot
minY = min([0 minVoltage*1.1]);
yLims = [minY maxY]; % The limits of the y axis on the plot

if doAlignSpikes
    minTime = inf;
    maxTime = -inf;
end

% If there are stimulations, plot them, plot in bold if a spike occurred
for curStimInd = 1:numPlottedTraces
    % Get index in terms of voltageTraces index (because curSelectedStimInd
    % is an index in selectedStims{neuronNum}, which refers to an intex in
    % voltageTraces
    traceInd = handles.selectedStims{neuronNum}(curStimInd);
    pcaSortCodes = handles.allNeuronsStruct{neuronNum}.pcaSortCodes{traceInd};
    
    % Create x axis
    if doPlotRaw
        % If plotting the raw trace
        trace = handles.allNeuronsStruct{neuronNum}.voltageTracesRaw{traceInd}; % Trace to plot
    else
        % If plotting the filtered trace
        trace = handles.allNeuronsStruct{neuronNum}.voltageTraces{traceInd}; % Trace to plot
    end
    filteredTrace = handles.allNeuronsStruct{neuronNum}.voltageTraces{traceInd}; % Trace to perform calculations on (such as spike detection)
    times = linspace(0, handles.traceTime, length(trace));
    
    % Get the trigger voltage and trigger timespan
    voltageThresh = handles.allNeuronsStruct{neuronNum}.triggerVoltage(traceInd);
    voltageTimespan = handles.allNeuronsStruct{neuronNum}.triggerTimespan(traceInd,:);
    
    % Did a spike occur for this stim?
    if doPCA
        % Use PCA values
        [spikeOccurred, spikeHitInds, validResults] = detectSpikesFromTraces(filteredTrace, handles.traceTime, [handles.allNeuronsStruct{neuronNum}.triggerVoltage(traceInd) handles.allNeuronsStruct{neuronNum}.triggerTimespan(traceInd,:)], pcaSortCodes, 1, get(handles.testing, 'Value')); % Find the index that the trace has a result
    else
        % Do not use the PCA values, use only the voltage threshold
        [spikeOccurred, spikeHitInds, validResults] = detectSpikesFromTraces(filteredTrace, handles.traceTime, [handles.allNeuronsStruct{neuronNum}.triggerVoltage(traceInd) handles.allNeuronsStruct{neuronNum}.triggerTimespan(traceInd,:)], [], 1, get(handles.testing, 'Value')); % Find the index that the trace has a result
    end
    handles.allNeuronsStruct{neuronNum}.resultsData(traceInd,3) = spikeOccurred; % Update the resultsData structure
    thisCollectionResults(curStimInd) = spikeOccurred == 1; % If a spike occurred, record it in thisCollectionResults
    
    % If a spike occurred, align the spikes if the user wants
    if spikeOccurred && doAlignSpikes
        validSpikeHitInds = spikeHitInds(validResults); % Get the indices of the spikes within the timespan
        if ~isempty(validSpikeHitInds)
            alignmentCorrection = times(validSpikeHitInds(1)); % Align to the first valid spike
            times = times - alignmentCorrection; % Align all times to 0
            voltageTimespan = voltageTimespan - alignmentCorrection; % Align the voltage timespan values according to this alignment correction
        end
        
        % Update the minimum and maximum time points
        minTime = min([minTime times(1)]);
        maxTime = max([maxTime times(end)]);
    end
    
    % Plot the trace, using the proper hsv color from the colormap, onto
    % this neuron's voltageTrace axis, bolding the line if there was a
    % spike
    plot(axisHandle, times, trace, 'Color', cMap(curStimInd,:), 'LineWidth', 1 + (spikeOccurred && ~doAlignSpikes)); % Make the line thicker if a spike occcurred, and spikes are NOT being aligned
    hold on;
    
    % Plot the timespan and voltage threshold markers
    if ~doAlignSpikes
        % Only plot the threshold markers if they are not being aligned.
        % Otherwise, it gets way too busy.
        
        % If there is only one trace, make it cyan.  If there are multiple,
        % make them match
        if numPlottedTraces == 1
            spikeBoxColor = 'c';
        else
            spikeBoxColor = cMap(curStimInd,:);
        end
        
        % Plot the voltage timespan markers
        plot(axisHandle,[voltageTimespan(1) voltageTimespan(1)], yLims, '--', [voltageTimespan(2) voltageTimespan(2)], yLims, '--', 'LineWidth', 1.5, 'Color', spikeBoxColor);
        
        % If this trace fired a PCA-detected spike (even if it is outside
        % of the timespan), indicate that to the user by displaying the
        % location of the spike
        if doPCA && ~isempty(pcaSortCodes)
            % If the user wants PCA as the spike detection choice, and PCA
            % information is available for this spike
            % START HERE: MODIFY DETECTSPIKESFROMTRACES SO THAT IT OUTPUTS
            % THE SORTCODE AT EACH SPIKE LOCATION
            for pcaCode = 1:maxPCACode
                if pcaCode ~= 1
                    [~, spikeHitInds, validResults] = detectSpikesFromTraces(filteredTrace, handles.traceTime, [handles.allNeuronsStruct{neuronNum}.triggerVoltage(traceInd) handles.allNeuronsStruct{neuronNum}.triggerTimespan(traceInd,:)], pcaSortCodes, pcaCode, get(handles.testing, 'Value')); % Find the index that the trace has a result
                    c = pcaCMap(pcaCode - 1,:);
                    valPlotStr = 'o';
                    plotWidth = 2;
                    plotSize = 12;
                else
                    c = [0 0 0];
                    valPlotStr = 'x';
                    plotWidth = 3;
                    plotSize = 15;
                end
                validSpikeHitInds = spikeHitInds(validResults);
                invalidSpikeHitInds = spikeHitInds(~validResults);
                if ~isempty(spikeHitInds)
                    if numPlottedTraces == 1
                        for spikeHitInd = 1:length(validSpikeHitInds)
                            plot(axisHandle, times(validSpikeHitInds(spikeHitInd)), 0, valPlotStr, 'linewidth', plotWidth, 'markerSize', plotSize, 'color', c); %, 'ButtonDownFcn', @(object, eventData) spikeMarkerCallback(object, handles), 'UserData', neuronNum);
                        end
                        for spikeHitInd = 1:length(invalidSpikeHitInds)
                            plot(axisHandle, times(invalidSpikeHitInds(spikeHitInd)), 0, 'o', 'linewidth', plotWidth, 'markerSize', plotSize, 'color', c); %, 'ButtonDownFcn', @(object, eventData) spikeMarkerCallback(object, handles), 'UserData', neuronNum);
                        end
                    else
                        %                         plot(axisHandle, times(spikeHitInds), 0, 'kx', 'linewidth', 3, 'markerSize', 15);
                        plot(axisHandle, times(validSpikeHitInds), zeros(1, length(validSpikeHitInds)), valPlotStr, 'linewidth', plotWidth, 'markerSize', plotSize, 'color', spikeBoxColor);
                        plot(axisHandle, times(invalidSpikeHitInds), zeros(1, length(invalidSpikeHitInds)), 'o', 'linewidth', plotWidth, 'markerSize', plotSize, 'color', spikeBoxColor);
                    end
                end
            end
        else
            % Plot the voltage threshold marker
            plot(axisHandle, [voltageTimespan(1) voltageTimespan(2)], [voltageThresh voltageThresh], '--', 'LineWidth', 1.5, 'Color', spikeBoxColor);
        end
    end
end

% Correct for not calculating the min/max time
if ~isfinite(minTime)
    minTime = 0;
end

if ~isfinite(maxTime)
    maxTime = handles.traceTime;
end

% If the spikes are aligned, plot a marker showing where they align to
if doAlignSpikes
    plot([0 0], [minY maxY], 'k:');
end

% Plot the currently selected trigger voltage and timespan (if the current
% value has not been applied yet)
voltageThresh = handles.currentTriggerVoltage(neuronNum); % Selected threshold
voltageTimespan = handles.currentTriggerTimespan(neuronNum,:); % Selected timespan

% Plot new voltage threshold
if all(~isnan(voltageThresh)) && ~doAlignSpikes
    % If the voltage threshold (selected via the "Set Trigger Voltage"
    % button) is not NaN (indicating that the new value has not yet been
    % applied), then plot it
    plot(axisHandle, [minTime maxTime], [voltageThresh voltageThresh], 'b--', 'LineWidth', 3);
end

% Plot new timespans
if ~isnan(voltageTimespan(1)) && ~doAlignSpikes
    % If the selected trigger parameters (selected via the "Set Timespan"
    % button) are not NaN (indicating that the new values have not yet been
    % applied), then plot them
    plot(axisHandle, [voltageTimespan(1) voltageTimespan(1)], yLims, 'b--', 'LineWidth', 3);
end

if ~isnan(voltageTimespan(2)) && ~doAlignSpikes
    % If the selected trigger parameters (selected via the "Set Timespan"
    % button) are not NaN (indicating that the new values have not yet been
    % applied), then plot them
    plot(axisHandle, [voltageTimespan(2) voltageTimespan(2)], yLims, 'b--', 'LineWidth', 3);
end

% Format the axis size
axis(axisHandle, [minTime maxTime minY maxY]);

% Go through the thisCollectionResults and see a) no spikes occurred, b)
% all stimulations lead to spikes, or c) there are some spikes and some no-
% spikes
if (all(~thisCollectionResults))
    % No spikes occurred
    color = handles.red; % Red warning, no spike occured
    s = 'No Spike';
    if numPlottedTraces > 1
        s = [s 's']; % Pluralize, because I'm fancy like that
    end
elseif all(thisCollectionResults)
    % All stimulations lead to spikes
    color = handles.brightGreen; % Green highlight, spikes were found
    s = 'Spike';
    if numPlottedTraces > 1
        s = [s 's']; % Pluralize, because I'm fancy like that
    end
else
    % Some spikes and some no-spikes
    color = handles.yellow;
    s = 'Mixed';
end

% Reset the mouse click callback for this axis and its children, because it
% was changed during the plotting routines
allLineHandles = [axisHandle; allchild(axisHandle)];
blankCallbacks = strcmp(get(allLineHandles, 'buttondownfcn'), '');
set(allLineHandles(blankCallbacks), 'ButtonDownFcn', @(object, eventData) figureClickCallback(axisHandle));

% Set color and string to the stimFlag for this neuron
set(handles.(['stimFlag' neuronNumStr]), 'BackgroundColor', color, 'String', s);

function waitLoop(timeInMS)
% This function waits a given time in milliseconds
tic;
while (toc < timeInMS/1000)
end


function handles = drawSDPlot(handles)
% TODO: Give indication of noise/width of curve?  I.e. plot the 50%
% probability line, as well as the 25% and 75% as dotted lines? %
useFPForPlotting = true; % Should we use the FP data to generate this curve?
isoLineVal = .5; % What percentage isoline on the SD surface should we plot?

% This function draws the Strenth-Duration curve of neuron number
% neuronNum, if it has already been characterized.
S1 = str2double(get(handles.controlStrength1, 'String'));
D1 = str2double(get(handles.controlDuration1, 'String'));
S2 = str2double(get(handles.controlStrength2, 'String'));
D2 = str2double(get(handles.controlDuration2, 'String'));

if useFPForPlotting
    allStrs = handles.allStrs;
    allDurs = handles.allDurs;
    allAlphas = handles.allAlphas;
    allSigmas = handles.allSigmas;
    PAll = handles.PAll;
    b0 = handles.b0;
    approxFun = @(G, T, t) lininterpn(allStrs, allDurs, allAlphas, allSigmas, PAll, G*(t(2)/b0), T, t(1), t(3));
    
    plotDurs = logspace(log10(1), log10(15), 6);
    GLB = zeros(size(plotDurs));
    GUB = max(handles.allStrs)*ones(size(plotDurs));
    GOpt = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'none');
    G0 = .1*ones(size(plotDurs));
else
    fSDCurve = @(p, g) (g>((p(1)*p(3))/(p(2)*(1 - p(3))))).*abs(((-1./(p(1) + g.*p(2))).*log(1 - (((p(1) + g.*p(2)).*p(3))./(g.*p(2)))))) + (g<=((p(1)*p(3))/(p(2)*(1 - p(3))))).*0; % TODO Does this need to change? % Function of the strength-duration curve
end

bigMarkerSize = 15;
smallMarkerSize = 5;
axisMinDur = 15; % Minimum duration to fit the axis to
axisMinStr = handles.maxIrradiance; % Minimum strength to fit the axis to
hold off;
for neuronNum = 1:handles.numNeurons
    % Get other neuron's number
    otherNeuronNum = mod(neuronNum, 2) + 1; % Number of the OTHER neuron
    
    % Get axis handle
    axisHandle = handles.(['SDCurve' num2str(neuronNum)]); % Get the handle of the axis to be modified
    cla(axisHandle); % Clear axis
    
    % Draw raw results data on the plot
    hold(axisHandle, 'on'); % Put this axis into hold on state
    spikeInds = handles.allNeuronsStruct{neuronNum}.resultsData(:,3) == 1; % Indices in resultsData that there was a spike
    noSpikeInds = handles.allNeuronsStruct{neuronNum}.resultsData(:,3) == 0; % Indices in resultsData that there was a spike
    selectedInds = false(size(handles.allNeuronsStruct{neuronNum}.resultsData, 1),1); % Start the creation of selectedInds
    selectedInds(handles.selectedStims{neuronNum}) = true; % Make a logical array of the stimulations that are currently selected in the listbox
    
    % Draw the SD curve and associated data of this and/or the other neuron
    if ~isempty(handles.allNeuronsStruct{neuronNum}.parameters)
        SDParams = handles.allNeuronsStruct{neuronNum}.parameters(end, :); % Parameters for the curve
        
        % If the neuron has been characterized
        precalcData = handles.lastDrawSDCurve{neuronNum};
        if isempty(precalcData) || ~all(handles.lastDrawSDParams(neuronNum, :) == SDParams)
            % Get needed data
            %         durations = handles.allNeuronsStruct{neuronNum}.SDCurvePoints(:,2); % All durations that have been characterized up to this point
            %         strengths = handles.allNeuronsStruct{neuronNum}.SDCurvePoints(:,1); % All strengths that have been characterized up to this point
            
            % Prepare and plot the curve
            if useFPForPlotting
                plotStrs = zeros(size(plotDurs));
%                 tic;
                for strInd = 1:length(plotDurs)
                    plotStrs(strInd) = fzero(@(S) approxFun(S, plotDurs(strInd), SDParams) - isoLineVal, G0(strInd), optimset('Display', 'off'));
                end
%                 toc;
            else
                gRheo = (SDParams(1)*SDParams(3))/(SDParams(2)*(1 - SDParams(3))); % Rheobase of this neuron
                plotStrs = linspace(gRheo*1.000000001, handles.maxStrength, 1000); % Strengths that will be used (to make a smooth curve)
                plotDurs = fSDCurve(SDParams, plotStrs);
            end
            
            % Save the parameters of the curve that were plotted, so we don't need to
            % waste time calculating the strengths/durations again next time
            handles.lastDrawSDParams(neuronNum, :) = SDParams;
            handles.lastDrawSDCurve{neuronNum} = [plotStrs;plotDurs];
        else
            plotStrs = precalcData(1, :);
            plotDurs = precalcData(2, :);
        end
        
        a = plot(axisHandle, plotDurs, plotStrs, 'linewidth', 2);
        set(a, 'color', handles.neuronPlotColor{neuronNum});
    end
    
    % Plot the other neuron's SD curve data
    if ~isempty(handles.allNeuronsStruct{otherNeuronNum}.parameters)
        SDParamsOther = handles.allNeuronsStruct{otherNeuronNum}.parameters(end, :); % Parameters for the other curve
        
        precalcData = handles.lastDrawSDCurve{otherNeuronNum};
        if isempty(precalcData) || ~all(handles.lastDrawSDParams(otherNeuronNum, :) == SDParamsOther)
            % If the OTHER neuron's SDCurvePoints is not empty (meaning it has a
            % curve plotted now), plot the other neuron's curve onto this plot, for
            % easy comparison

            % Prepare and plot the curve
            if useFPForPlotting
                plotStrs = zeros(size(plotDurs));
                for strInd = 1:length(plotDurs)
                    plotStrs(strInd) = fzero(@(S) approxFun(S, plotDurs(strInd), SDParamsOther) - isoLineVal, G0(strInd), optimset('Display', 'off'));
                end
            else
                gRheo = (SDParamsOther(1)*SDParamsOther(3))/(SDParamsOther(2)*(1 - SDParamsOther(3))); % Rheobase of this neuron
                plotStrs = linspace(gRheo*1.000000001, handles.maxStrength, 1000); % Strengths that will be used (to make a smooth curve)
                plotDurs = fSDCurve(SDParamsOther, plotStrs);
            end
            
            % Save the parameters of the curve that were plotted, so we don't need to
            % waste time calculating the strengths/durations again next time
            handles.lastDrawSDParams(otherNeuronNum, :) = SDParamsOther;
            handles.lastDrawSDCurve{otherNeuronNum} = [plotStrs;plotDurs];
        else
            plotStrs = precalcData(1, :);
            plotDurs = precalcData(2, :);
        end
        a = plot(axisHandle, plotDurs, plotStrs, 'linewidth', 2);
        set(a, 'color', handles.neuronPlotColor{otherNeuronNum});
    end
    
    % Plot the control parameters for both neurons
    plot(axisHandle, D1, S1, '+w', 'markersize', bigMarkerSize + 2, 'linewidth', 4);
    plot(axisHandle, D2, S2, '+w', 'markersize', bigMarkerSize + 2, 'linewidth', 4);
    plot(axisHandle, D1, S1, '+', 'markersize', bigMarkerSize, 'linewidth', 2, 'color', handles.neuronPlotColor{1});
    plot(axisHandle, D2, S2, '+', 'markersize', bigMarkerSize, 'linewidth', 2, 'color', handles.neuronPlotColor{2});
    
    % Create all four logical arrays that will be used to plot the four
    % different types of data
    spikeAndSelectedInds = spikeInds & selectedInds; % Large black X's
    noSpikeAndSelectedInds = noSpikeInds & selectedInds; % Large red squares
    spikeAndNotSelectedInds = spikeInds & ~selectedInds; % Small black X's
    noSpikeAndNotSelectedInds = noSpikeInds & ~selectedInds; % Small red squares
    
    % Plot all four groups of stimulations (grouped by the size of the
    % marker)
    plot(axisHandle, handles.allNeuronsStruct{neuronNum}.resultsData(spikeAndNotSelectedInds,2),handles.allNeuronsStruct{neuronNum}.resultsData(spikeAndNotSelectedInds,1), 'xk', handles.allNeuronsStruct{neuronNum}.resultsData(noSpikeAndNotSelectedInds,2), handles.allNeuronsStruct{neuronNum}.resultsData(noSpikeAndNotSelectedInds,1), 'sr', 'markersize', smallMarkerSize); % Plot all not selected stimulations
    plot(axisHandle, handles.allNeuronsStruct{neuronNum}.resultsData(spikeAndSelectedInds, 2), handles.allNeuronsStruct{neuronNum}.resultsData(spikeAndSelectedInds, 1), 'xk', handles.allNeuronsStruct{neuronNum}.resultsData(noSpikeAndSelectedInds, 2), handles.allNeuronsStruct{neuronNum}.resultsData(noSpikeAndSelectedInds, 1), 'sr', 'markersize', bigMarkerSize, 'linewidth', 2); % Plot all selected stimulations
    
    % Get the min and max for both strength and duration
    minMaxDur = [0 1.1*max([axisMinDur;handles.allNeuronsStruct{neuronNum}.resultsData(:,2)])];
    minMaxStr = [0 1.1*max([axisMinStr;handles.allNeuronsStruct{neuronNum}.resultsData(:,1)])];
    
    % Display the irradiance boundary, if relevant
    if (handles.maxIrradianceStrength ~= handles.maxStrength)
        % If the max strength allowed from the irradiance boundary is
        % different than the normal maximum allowed strength, then show
        % this in the SD plot
        plot(axisHandle, minMaxDur, handles.maxIrradianceStrength([1 1]), 'r:', 'linewidth', .5);
    end
    
    % Set axis to size required to display the raw data
    axis(axisHandle, [minMaxDur minMaxStr]); % Make sure to fit all durations/strengths, but make the plot no smaller than [0 15 0 5]
    
    hold(axisHandle, 'off'); % Remove hold
    
    % Reset the mouse click callback for this axis and its children, because it
    % was changed during the plotting routines
    set([axisHandle; allchild(axisHandle)], 'ButtonDownFcn', @(object, eventData) figureClickCallback(axisHandle));
end

drawnow nocallbacks



function alpha1_Callback(hObject, eventdata, handles)
% hObject    handle to alpha1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha1 as text
%        str2double(get(hObject,'String')) returns contents of alpha1 as a double


% --- Executes during object creation, after setting all properties.
function alpha1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta1_Callback(hObject, eventdata, handles)
% hObject    handle to beta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta1 as text
%        str2double(get(hObject,'String')) returns contents of beta1 as a double

% --- Executes during object creation, after setting all properties.
function beta1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w1_Callback(hObject, eventdata, ~)
% hObject    handle to w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w1 as text
%        str2double(get(hObject,'String')) returns contents of w1 as a double

% --- Executes during object creation, after setting all properties.
function w1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = fitSDCurveCallback(handles, varargin)
% Fit thetas to all neurons
if nargin > 1
    charTime = varargin{1}; % The time of the most recent characterization.  Only data after this time will be used in the fit
    handles = calculateOptimalThetas(handles, false, charTime);
else
    handles = calculateOptimalThetas(handles, false);
end

% Save these new GT's to the GUI
for neuronNum = 1:handles.numNeurons
    set(handles.(sprintf('controlStrength%d', neuronNum)), 'String', num2str(handles.curExpGT(neuronNum, 1)));
    set(handles.(sprintf('controlDuration%d', neuronNum)), 'String', num2str(handles.curExpGT(neuronNum, 2)));
end
handles = setControlVoltages(handles);

% Redraw the SD curves
handles = drawSDPlot(handles);

% Auto save
handles = autosave(handles);

% Notify the user that the fitting is done
playWarn(handles);

% --- Executes on button press in TDTConnectButton.
function TDTConnectButton_Callback(hObject, eventdata, handles)
% hObject    handle to TDTConnectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable the button
set(handles.TDTConnectButton, 'Enable', 'off');

% Toggle the state of the connection
if (handles.TDTConnected)
    % If there is a connection, end it
    handles = endTDTComm(handles);
else
    % If there is no connection, start one
    handles = startTDTComm(handles);
end

% Check the TDT connection and update the GUI
handles = checkTDTConnectivity(handles);

% Re-enable the button
set(handles.TDTConnectButton, 'Enable', 'on');

guidata(hObject, handles);

% --- Executes on button press in setTriggerTimespan1.
function setTriggerTimespan1_Callback(hObject, eventdata, handles)
% hObject    handle to setTriggerTimespan1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = setTriggerTimespan(hObject, handles);
guidata(hObject, handles);


function handles = setTriggerTimespan(hObject, handles)
neuronNum = get(hObject, 'UserData');

if handles.stimsExist(neuronNum)
    neuronNumStr = num2str(neuronNum); % Neuron number as a string
    axisHandle = handles.(['voltageTrace' neuronNumStr]); % Get the handle of the axis to be modified
    axes(axisHandle); % Set the voltage trace axis as the current axis
    newTimespan = zeros(1,2); % The variable to which the timespan will be temporarily assigned
    
    [newTimespan(1), ~, button] = ginput(1); % Get the first trigger timespan, and set it temporary variable
    if button ~= 1
        % If the user didn't click the mouse, and instead hit a key, assume
        % that they want to cancel the selection
        return;
    end
    hold on;
    plot(axisHandle, [newTimespan(1) newTimespan(1)], get(axisHandle, 'YLim'), 'm:', 'LineWidth', 3); % Plot the first part of the timespan onto the axis, as a visual aid
    
    [newTimespan(2), ~, button] = ginput(1); % Get the second trigger timespan, and set it temporary variable
    if button ~= 1
        % If the user didn't click the mouse, and instead hit a key, assume
        % that they want to cancel the selection
        handles = drawVoltageTraces(handles, neuronNum); % Clear the previously added line
        return;
    end
    
    newTimespan = sort(newTimespan); % Sort the timespan, so that the lower number is always first
    handles.currentTriggerTimespan(neuronNum,:) = newTimespan; % Assign the new timespan to be the currently selected timespan
    handles = drawVoltageTraces(handles, neuronNum); % Draw the plot
end

updateTimespanBoxes(handles); % Record the new timespan boxes in the GUI

% Update handles
guidata(hObject, handles);

% --- Executes on selection change in neuron1StimsListbox.
function neuron1StimsListbox_Callback(hObject, eventdata, handles)
% hObject    handle to neuron1StimsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns neuron1StimsListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from neuron1StimsListbox

handles = stimsListboxCallback(handles, get(hObject, 'Value'));
guidata(hObject, handles);

function handles = stimsListboxCallback(handles, varargin)
% CHECKED %

inputtedValues = false;
if nargin > 1
   values = varargin{1};
   inputtedValues = true;
end

for neuronNum = 1:handles.numNeurons
    listBoxObject = handles.(sprintf('neuron%dStimsListbox', neuronNum));
    % Assign the currently selected values to handles.selectedStims
    if ~inputtedValues
        values = get(listBoxObject, 'Value'); % The indices in the listbox that are currently selected, which correspond to the rows of resultsData and other corresponding data structures
    end
    
    if handles.stimsExist(neuronNum)
        % If some stimulations exist
        handles.selectedStims{neuronNum} = values; % Assign the selected values
        set(listBoxObject, 'Value', values);
    else
        % If no stimulations exist, then the listbox is only populated by the "No Stims" flag
        handles.selectedStims{neuronNum} = [];
    end
    
    % Draw the voltage trace of the stimulation
    handles = drawVoltageTraces(handles, neuronNum);
end

% Finally, draw the SD plot (to see the newly selected stimulations on the
% SD graph)
handles = drawSDPlot(handles); % Done last, because it requires both neuron structs to be up-to-date

updateTimespanBoxes(handles);

% --- Executes during object creation, after setting all properties.
function neuron1StimsListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron1StimsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in applyTrigWindowSelected1.
function applyTrigWindowSelected1_Callback(hObject, eventdata, handles)
% hObject    handle to applyTrigWindowSelected1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

neuronNum = get(hObject, 'UserData'); % The number of the selected neuron
applyTrigParamsCallback(hObject, handles, handles.selectedStims{neuronNum}) % Apply the parameters to only the selected neurons

function applyTrigParamsCallback(hObject, handles, spikeInds)
% CHECKED %

% This function takes the trigger parameters (voltage threshold and
% timespan) and applies them to the stimulations given by spikeInds.  It
% will ALSO apply these parameters to any future stimulations. spikeInds is
% every index that will acquire the current trigger voltage and current
% timespan as their own, to be recorded in the neuronStruct. The "Apply to
% Selected" button applies them only to the stimulations selected in the
% stimsListbox.  The "Apply to All" button applies them to all stimulations

neuronNum = get(hObject, 'UserData'); % Get selected neuron number
numSelectedStims = length(spikeInds); % Number of stims

curResultsLoc = handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc - 1;

if curResultsLoc >= 1
    % For the currently selected trigger voltage
    if ~any(isnan(handles.currentTriggerVoltage(neuronNum)))
        % Get the new trigger voltage
        newVoltage = handles.currentTriggerVoltage(neuronNum);
        
        % Apply the trigger voltage to all selected responses
        handles.allNeuronsStruct{neuronNum}.triggerVoltage(spikeInds) = newVoltage;
        
        % Apply to the neuronStruct (for future responses)
        handles.allNeuronsStruct{neuronNum}.defaultTriggerSettings(1) = newVoltage;
        
        % Apply to GUITriggerSetings (for future neurons)
        handles.GUITriggerSettings(neuronNum,1) = newVoltage;
    end
    
    % For the currently selected trigger timespan
    if any(~isnan(handles.currentTriggerTimespan(neuronNum, :)))
        % Get the new trigger timespan
        newTimespan = handles.currentTriggerTimespan(neuronNum,:);
        
        % If any of the new timespan values are NaN, then extract the currently
        % used value from the corresponding neuronStruct
        isNaNInd = isnan(newTimespan);
        newTimespan(isNaNInd) = handles.allNeuronsStruct{neuronNum}.triggerTimespan(curResultsLoc, isNaNInd);
        
        % Apply the trigger timespan to all selected responses
        handles.allNeuronsStruct{neuronNum}.triggerTimespan(spikeInds,:) = repmat(newTimespan, numSelectedStims, 1);
        
        % Apply to the neuronStruct (for future responses)
        handles.allNeuronsStruct{neuronNum}.defaultTriggerSettings(2:3) = newTimespan;
        
        % Apply to GUITriggerSetings (for future neurons)
        handles.GUITriggerSettings(neuronNum,2:3) = newTimespan;
    end
    
    % Clear the currently selected trigger voltage and timespans
    handles = clearSelectedTriggerParams(handles, neuronNum);
    
    % Check the results of the selected stimulations using the new parameters
    handles = recheckResults(handles, neuronNum, spikeInds);
    
    % Update the timespan boxes
    updateTimespanBoxes(handles);
    
    % Update handles
    guidata(hObject, handles);
end

function handles = recheckResults(handles, neuronNum, spikeInds)
% CHECKED %
% This function goes through the select stimulation indicies
% (given by spikeInds) and reevaluates whether the stimulation resulted in
% a spike or not

% Prepare for looping through all selected stimulations
numSelectedStims = length(spikeInds); % Number of stims

% Revaluate the voltage traces of the selected spikes
% Go through each of the currently selected spikes, and evaluate if there
% is a spike or not
for curSelectedStimIndInd = 1:numSelectedStims
    curStimInd = spikeInds(curSelectedStimIndInd); % Get the index of this selected spike relative to the resultsData matrix
    
    % Re-run the detect spikes routine on the saved voltage trace data,
    % using the new voltage threshold and timespan
    if get(handles.(sprintf('usePCACheck%d', neuronNum)), 'Value')
        handles.allNeuronsStruct{neuronNum}.resultsData(curStimInd, 3) = detectSpikesFromTraces(handles.allNeuronsStruct{neuronNum}.voltageTraces{curStimInd}, handles.traceTime, [handles.allNeuronsStruct{neuronNum}.triggerVoltage(curStimInd) handles.allNeuronsStruct{neuronNum}.triggerTimespan(curStimInd,:)], handles.allNeuronsStruct{neuronNum}.pcaSortCodes{curStimInd}, 1, get(handles.testing, 'Value'));
    else
        handles.allNeuronsStruct{neuronNum}.resultsData(curStimInd, 3) = detectSpikesFromTraces(handles.allNeuronsStruct{neuronNum}.voltageTraces{curStimInd}, handles.traceTime, [handles.allNeuronsStruct{neuronNum}.triggerVoltage(curStimInd) handles.allNeuronsStruct{neuronNum}.triggerTimespan(curStimInd,:)], [], 1, get(handles.testing, 'Value'));
    end
end

% Draw the voltage traces of the selected stimulations
handles = drawVoltageTraces(handles, neuronNum);

% Draw the newly interpretted data (and SD curve, if it has been calculated, even though it may be outdated at this point)
handles = drawSDPlot(handles);

function handles = clearSelectedTriggerParams(handles, neuronNum)
% CHECKED %
% Set the trigger parameters for this neuron to be nan, indicating that
% they should not be plotted
handles.currentTriggerVoltage(neuronNum) = nan;
handles.currentTriggerTimespan(neuronNum,:) = nan(1,2);

updateTimespanBoxes(handles);

function updateTimespanBoxes(handles)
% Update the GUI's text boxes
% Plot the actual trigger timespan for each text box, unless there is a
% suggested trigger timespan

if ~isempty(handles.selectedStims{1})
    % If any stimulations exist
    thisLoc1 = handles.selectedStims{1}(end);
    thisLoc2 = handles.selectedStims{2}(end);
    
    normalC = [1 1 1];
    alteredC = [.75 .75 .05];
    
    if thisLoc1 >= 1
        % If some data has actually be gathered
        if ~isnan(handles.currentTriggerTimespan(1, 1))
            set(handles.tsTextL_1, 'String', num2str(handles.currentTriggerTimespan(1, 1)), 'BackgroundColor', alteredC);
        else
            set(handles.tsTextL_1, 'String', num2str(handles.allNeuronsStruct{1}.triggerTimespan(thisLoc1, 1)), 'BackgroundColor', normalC);
        end
        
        if ~isnan(handles.currentTriggerTimespan(1, 2))
            set(handles.tsTextH_1, 'String', num2str(handles.currentTriggerTimespan(1, 2)), 'BackgroundColor', alteredC);
        else
            set(handles.tsTextH_1, 'String', num2str(handles.allNeuronsStruct{1}.triggerTimespan(thisLoc1, 2)), 'BackgroundColor', normalC);
        end
    end
    
    if thisLoc2 >= 1
        if ~isnan(handles.currentTriggerTimespan(2, 1))
            set(handles.tsTextL_2, 'String', num2str(handles.currentTriggerTimespan(2, 1)), 'BackgroundColor', alteredC);
        else
            set(handles.tsTextL_2, 'String', num2str(handles.allNeuronsStruct{2}.triggerTimespan(thisLoc2, 1)), 'BackgroundColor', normalC);
        end
        
        if ~isnan(handles.currentTriggerTimespan(2, 2))
            set(handles.tsTextH_2, 'String', num2str(handles.currentTriggerTimespan(2, 2)), 'BackgroundColor', alteredC);
        else
            set(handles.tsTextH_2, 'String', num2str(handles.allNeuronsStruct{2}.triggerTimespan(thisLoc2, 2)), 'BackgroundColor', normalC);
        end
    end
end

% --- Executes on button press in applyTrigWindowAll1.
function applyTrigWindowAll1_Callback(hObject, eventdata, handles)
% hObject    handle to applyTrigWindowAll1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

neuronNum = get(hObject, 'UserData'); % The number of the selected neuron
numStims = handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc - 1; % The number of stimulations that are currently being tracked for this neuron
applyTrigParamsCallback(hObject, handles, 1:numStims) % Apply the parameters to only the selected neurons


function manStimStrength_Callback(hObject, eventdata, handles)
% hObject    handle to manStimStrength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of manStimStrength as text
%        str2double(get(hObject,'String')) returns contents of manStimStrength as a double

% Check the input values
handles = checkNewInput(hObject, handles, 'manStrengthVarData');

handles = setControlVoltages(handles);
set(handles.controlV, 'String', num2str(handles.curManControlVoltage));

% If TDT is connected, send this value to the system
if (handles.TDTConnected)
    handles.DA.SetTargetVal('RZ2.Matlab Strength', str2double(get(hObject, 'String')));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function manStimStrength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to manStimStrength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function manStimDuration_Callback(hObject, eventdata, handles)
% hObject    handle to manStimDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of manStimDuration as text
%        str2double(get(hObject,'String')) returns contents of manStimDuration as a double

% Check the input values
handles = checkNewInput(hObject, handles, 'manDurationVarData');

% If TDT is connected, send this value to the system
if (handles.TDTConnected)
    handles.DA.SetTargetVal('RZ2.Matlab Duration', str2double(get(hObject, 'String')));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function manStimDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to manStimDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in manStim.
function manStim_Callback(hObject, eventdata, handles)
% hObject    handle to manStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get strength and duration from the GUI
strength = str2double(get(handles.manStimStrength, 'String'));
duration = str2double(get(handles.manStimDuration, 'String'));

manStim(strength, duration, hObject, handles);


function manStim(strength, duration, hObject, handles)
% CHECKED - Can't think of anything to change...right? %
% Disable the GUI
disableGUI(handles);

% Stimulate
[handles, ~] = stimulate(handles, strength, duration); % Stimulate, but do NOT wait some amount of time

% Autosave
handles = autosave(handles);

% Re-enable the GUI
handles = updateGUI(handles);

% Update handles
guidata(hObject, handles);

function handles = autosave(handles)
% This function saves the data in the neuronsStruct object automatically
% when called, as long as 1) the autosave checkbox for the corresponding
% unit is checked, and 2) there already exists a save file for the unit
% (generated by manually saving with the "Save" button).

if get(handles.autoSave, 'Value') && ~isempty(handles.autoSaveFiles)
    % If the autosave feature is on, and the user has already manually saved with this unit
    
    % Save this unit's data
    handles = saveNeuronStruct(handles, handles.autoSaveFiles);
    
else
    checkSaveBannerState(handles) % See if the user should be notified that a change has been made to the neuronStruct
end

function results = solveIAFNeuronMultiple(handles, neuronParams, strength, duration, tF, numDataPoints)
% This function will call a subfunction to solve an IAF modeled neuron with
% the provided parameters. neuronParams is a cell array that contains the
% [alpha, beta, vThresh, E] parameters for some number of neurons. strength
% is the amplitude at which the neuron will be stimulated, and duration is
% the length of time of the stimulation.  tF is the entire length that will
% be displayed on the GUI. numDataPoints is the number of data points
% requested in the results vector for the signal. results is a cell array
% of the same length as neuronParams that contains the voltage traces from
% the neurons. USED ONLY FOR TESTING.

useFPModel = true; % Should the Fokker-Planck model be used to simulate an extra-cellular spike trace?

% Set up parameters for the simulation
newTimePoints = linspace(0,tF, numDataPoints); % Time points to interpolate the voltage to
TDTScale = 1.1; % Scale at which the TDT outputs data (the maximum value that test neuron data will have)
numNeurons = length(neuronParams); % Number of neurons to test
results = cell(numNeurons, 1); % The final results cell array
for neuronNum = 1:numNeurons
    % Set up the voltage scale
    scaleFactor = TDTScale/neuronParams{neuronNum}(4); % Factor by which the output from solveIAFNeuron must be scaled to match the output from TDT
    
    % Sample from the model neuron
    if useFPModel
        %% Use the FP lookup table
        sigma = .05; % Intrinsic noise to each neuron
        noise = .1; % Amount of noise to include with the simulation
        spikeLocFac = .8; % The spike should be located around this fraction of the way through the duration of the stimulation
        spikeLocSTD = .2; % This amount of variation on the spike location
        
        % Create a baseline voltage trace
        spikeLength = 32;
        thisNeuronVoltages = noise*rand(1, numDataPoints);
        fSpikeWaveform = @(x) [linspace(0, 1, round(x/4)) linspace(1, -.3, x/4) linspace(-.3, 0, x/2)];
        thisSpikeWaveform = fSpikeWaveform(spikeLength);
        
        % Should the neuron spike, or not?
        allStrs = handles.allStrs;
        allDurs = handles.allDurs;
        allAlphas = handles.allAlphas;
        allSigmas = handles.allSigmas;
        PAll = handles.PAll;
        b0 = handles.b0;
        t = [neuronParams{neuronNum}(1:2) sigma];
        thisP = lininterpn(allStrs, allDurs, allAlphas, allSigmas, PAll, strength*(t(2)/b0), duration, t(1), t(3));
        
        % Flip weighted coin to determine if neuron spiked or not
        thisSpike = rand < thisP;
        
        %%%% TEST
%         thisSpike = strength > 14;
        %%%% TEST
        
        if thisSpike
            spikeLoc = min(numDataPoints - spikeLength + 1, max(1,  round(duration*(spikeLocFac + 2*spikeLocSTD*(rand - .5))))); % Ensure that the spike occurs around spikeLocFac, but still within the bounds of the voltage trace
            thisNeuronVoltages(spikeLoc:(spikeLoc + spikeLength - 1)) = thisSpikeWaveform + thisNeuronVoltages(spikeLoc:(spikeLoc + spikeLength - 1)); % Add the spike waveform
        end
        
        % Scale up the voltage trace
        results{neuronNum} = scaleFactor*thisNeuronVoltages;
        
    else
        %% Use the deterministic IAF model
        gFun = @(t) (t < duration)*strength; % The function that solveIAFNeuron will use to solve the neuron integration (make a pulse starting at 0 of length "duration")
        
        [thisNeuronTimes, thisNeuronVoltages] = solveIAFNeuronOutputResults(tF, gFun, 0, neuronParams{neuronNum});
        
        % Remove any repeated time points
        repeatedPoints = find(diff(sort(thisNeuronTimes)) == 0); % Find all points where the time axis is repeated
        thisNeuronTimes(repeatedPoints) = []; % Remove these points
        thisNeuronVoltages(repeatedPoints) = []; % Remove these points
        
        % Resample the data so that it is evenly sampled
        results{neuronNum} = scaleFactor*interp1(thisNeuronTimes, thisNeuronVoltages, newTimePoints);
    end
end

function disableGUI(handles)
% Enables or disables buttons on the GUI that cause a stimulation, or
% that create a large calculation time
itemsToToggle = {'characterizeNeurons', 'fitSDCurve', 'manStim', 'neuron1StimsListbox', 'neuron2StimsListbox', 'TDTConnectButton', 'manStim1', 'manStim2', 'trialStart', 'swapGTs'};
for i = 1:length(itemsToToggle)
    set(handles.(itemsToToggle{i}), 'Enable', 'off');
end

function [handles, results] = stimulate(handles, strength, duration)
% TODO TDT: Get exact time from spike from TDT (relative to tank timing),
% if possible %

% This function is used to stimulate the neurons being recorded.  The
% strength and duration used will be that found in the manStimStrength and
% manStimDuration edit boxes.

results = zeros(handles.numNeurons,1);

% Extract information from TDT
% Set parameter needed for reading the buffer
sizeBuffer = handles.sizeBuffer;

% Preallocate spike data
electrodeData = cell(handles.numNeurons, 1); % Preallocate electrodeData, for the filtered data
electrodeDataRaw = cell(handles.numNeurons, 1); % Preallocate electrodeDataRaw, for the raw (unfiltered) data
sortCodes = cell(handles.numNeurons, 1); % Preallocate sortCodes, for the PCA data

isExploreCircuit = false; % Assume that we are working with a recording circuit (i.e. we should be expecting waveforms back from TDT after stimulating)

if get(handles.testing, 'Value')
    % If the user is just testing, only return a result based on the strength
    % and duration
    
    [electrodeData, electrodeDataRaw, sortCodes] = simulatedNeuron(strength, duration, handles);
else
    % If the user is NOT testing, send a pulse via TDT
    % Only stimulate if the TDT system is connected
    if (~handles.TDTConnected)
        error('TDT Not Connected');
    end
    
    DA = handles.DA; % Assign the shorthand variable for the TDT connection
    
    % TODO TDT: Check the TDT interactions %
    % Set variables to TDT
    % Get the control voltage that corresponds to the inputed power
    controlV = handles.fControlVFromIrr(strength, handles.curLaserParameters); % Find the control voltage that will produce the given power from the laser
    
    % Assign the strength and durations
    DA.SetTargetVal('RZ2.Matlab Strength', controlV);
    DA.SetTargetVal('RZ2.Matlab Duration', duration);
    
    % Clear any current results in TDT (yes, it's redundant, but it's
    % an easy fix)
    DA.SetTargetVal('RZ2.MatResultsRead', 1);
    DA.SetTargetVal('RZ2.MatResultsRead', 0);
    
    % Trigger the stimulation pulse for TDT
    DA.SetTargetVal('RZ2.Matlab Pulse', 1);
    DA.SetTargetVal('RZ2.Matlab Pulse', 0);
    
    % See if the TDT circuit being interfaced with right now is a silicon
    % probe circuit, and if it is the Explore circuit.
    isExploreCircuit = DA.GetTargetVal('RZ2.silExplore');
    isSilCircuit = DA.GetTargetVal('RZ2.isSil');
    isTetCircuit = isSilCircuit && DA.GetTargetVal('RZ2.isTet');
    
    % TODO: Testing
    sortCodesAsMask = true; % Should the tetrode circuit sort codes be interpretted as a mask? (Only for testing)
    
    % If this is not the explore circuit, wait for a response from TDT
    if ~isExploreCircuit
        % Wait for TDT to be done
        while (~DA.GetTargetVal('RZ2.ResultsReady'))
        end
        
        % Get voltage data from both electrodes
        electrodeData{1} = DA.ReadTargetVEX('RZ2.WaveForm1',0,sizeBuffer,'F32','F32');
        electrodeData{2} = DA.ReadTargetVEX('RZ2.WaveForm2',0,sizeBuffer,'F32','F32');
        electrodeDataRaw{1} = DA.ReadTargetVEX('RZ2.WaveFormRaw1',0,sizeBuffer,'F32','F32');
        electrodeDataRaw{2} = DA.ReadTargetVEX('RZ2.WaveFormRaw2',0,sizeBuffer,'F32','F32');
        
        % TODO SILICON: Set up for tetrode, too? %
        % Get PCA sort codes if they are being used
        if handles.TDTUsePCA
            for neuronNum = 1:handles.numNeurons
                incomingCodes = DA.ReadTargetVEX(sprintf('RZ2.sortCodes%d', neuronNum),0,sizeBuffer,'F32','F32');
                if isTetCircuit && sortCodesAsMask
                    %                     incomingCodes = DA.ReadTargetVEX(sprintf('RZ2.sortCodes%d', neuronNum),0,sizeBuffer,'I16','I16');
                    incomingCodes = uint16(incomingCodes); % Convert the float values to uin16, so they can be interpretted as a mask
                    
                    % Do post-processing on sort codes from silicon probe circuit (due to
                    % tetrodes)
                    
                    % The sortcodes from TetSort come as a 16-bit mask, indicating if each
                    % frame has a waveform in any of the sorting circles (see documentation
                    % for more details)
                    
                    % First, break sortCodes down into bits (make an N x L
                    % matrix, where L is the length of the incoming codes,
                    % where each value in the matrix represent the nth bit
                    % of the lth code)
                    N = 16;
                    allSortBits = bitget(repmat(incomingCodes, N, 1), repmat((1:N)', 1, length(incomingCodes)));
                    
                    % Find the least significant bit of each code (the
                    % "first" sort code that each waveform is sorted into).
                    %  We are throwing out some sorting information, but we
                    %  don't really need to know EVERY sort code that each
                    %  waveform falls into
                    [isValid, finalCodes] = max(allSortBits, [], 1);
                    % Make sure that we don't count any zero values
                    isValid = logical(isValid);
                    finalCodes(~isValid) = 0;
                    
                    % Store the final sort codes
                    sortCodes{neuronNum} = double(finalCodes);
                else
                    sortCodes{neuronNum} = incomingCodes;
                end
            end
        end
        
        % Tell TDT that the result has been read, so the next result will be ready
        DA.SetTargetVal('RZ2.MatResultsRead', 1);
        DA.SetTargetVal('RZ2.MatResultsRead', 0);
    end
end

% Record when this stimulation was done, and when the next one can be done
handles.lastStimTime = clock;

if strength == 0
    % If the stimulation had no strength (as is seen during characterization as a low boundary), do not pause after, as the laser did not fire at any strength
    handles.nextStimTime = handles.lastStimTime;
else
    %     stimTimeMinMax = [str2double(get(handles.seqPeriodLow, 'String')) str2double(get(handles.seqPeriodHigh, 'String'))];
    stimTimeMinMax = [.5 1];
    % If the stimulation had strength (as in most cases), add a random number of seconds between the two values in handles.stimMinMax to the last stim time to get the earliest next stim time
    handles.nextStimTime = handles.lastStimTime + [zeros(1,5) (stimTimeMinMax(1) + diff(stimTimeMinMax)*rand)];
end

% If we are not recording the waveforms from this stimulation, then we can
% safely exit this function
if isExploreCircuit
   return; 
end


% Record the time at which the data was read from TDT (relative to the
% reference start time in neuronStimStructure)
stimTime = etime(handles.lastStimTime, handles.startTime);

% Process information by TDT
% Can do threshold comparisons, perhaps within some time window.

for neuronNum = 1:handles.numNeurons
    if get(handles.(sprintf('usePCACheck%d', neuronNum)), 'Value')
        % If this unit should be using PCA, then determine the results
        % using the sort codes from TDT
        results(neuronNum) = detectSpikesFromTraces(electrodeData{neuronNum}, handles.traceTime, handles.allNeuronsStruct{neuronNum}.defaultTriggerSettings, sortCodes{neuronNum}, 1, get(handles.testing, 'Value'));
    else
        % If this unit should not be using PCA, then determine the results
        % using a voltage threshold
        results(neuronNum) = detectSpikesFromTraces(electrodeData{neuronNum}, handles.traceTime, handles.allNeuronsStruct{neuronNum}.defaultTriggerSettings, [], 1, get(handles.testing, 'Value'));
    end
end

% Package information to be output-ready
handles.allNeuronsStruct = addResultToNeuronStruct(handles.allNeuronsStruct, strength, duration, results, electrodeData, electrodeDataRaw, stimTime, sortCodes);

% Update the timespan textboxes for this latest stimulation
updateTimespanBoxes(handles);

% Update GUI with new neuron structs
handles = loadStimulationsIntoGUI(handles);

handles = drawSDPlot(handles); % Draw all new data (and SD curve, if it has been calculated)

function [electrodeData, electrodeDataRaw, sortCodes] = simulatedNeuron(strength, duration, handles)
% TODO SILICON: Include (possible all zeros) tetrode channels 2-4? %
sizeBuffer = handles.sizeBuffer;

% Set alpha, beta, vThreshold, and E for both neurons

% Close parameters:
% Alpha
alpha1_c = .5;
alpha2_c = .05;

% Beta
beta1_c = .5/handles.irrRatio;
beta2_c = .24/handles.irrRatio;

% Far parameters:
% Alpha
alpha1_f = .3;
alpha2_f = .05;

% Beta
% beta1 = .2; % Adjust the beta so that it accepts inputs in terms of laser power instead of control voltage
% beta2 = .12;
beta1_f = .3/handles.irrRatio; % Adjust the beta so that it accepts inputs in terms of laser power instead of control voltage
beta2_f = .24/handles.irrRatio;

if handles.runningTrial && handles.doTestingNeuronDrift
    maxNumRuns = str2double(get(handles.numRuns, 'String'));
    maxNumBlocks = (maxNumRuns*size(handles.curExpectedResults, 1))/handles.adaptiveBufferReplacementSize;
    ratio = (handles.blockNum - 1)/maxNumBlocks;
    
    alpha1 = interp1([0 1], [alpha1_c alpha1_f], ratio);
    alpha2 = interp1([0 1], [alpha2_c alpha2_f], ratio);
    beta1 = interp1([0 1], [beta1_c beta1_f], ratio);
    beta2 = interp1([0 1], [beta2_c beta2_f], ratio);
else
    % If there is no experiment running, then simply send the standard
    % neuron parameters
    alpha1 = alpha1_c;
    alpha2 = alpha2_c;
    beta1 = beta1_c;
    beta2 = beta2_c;
end

% VThreshold
vThresh1 = .5;
vThresh2 = .5;

% E
E1 = 1;
E2 = 1;
neuronParams{1} = [alpha1, beta1, vThresh1, E1];
neuronParams{2} = [alpha2, beta2, vThresh2, E2];

% Get simulated data
electrodeData = solveIAFNeuronMultiple(handles, neuronParams, strength, duration, handles.traceTime, sizeBuffer);
electrodeDataRaw = electrodeData; % For the simulated case, just take the straight simulated data as the "raw" data, as well as the "filtered" data
sortCodes = cellfun(@(x) [0 double(diff(x > 1)>0)], electrodeData, 'uniformOutput', 0); % Get the sortCodes by finding where the spike jump occurs

function newNeuronStruct = addResultToNeuronStruct(allNeuronsStruct, strength, duration, result, voltageTrace, voltageTraceRaw, time, varargin)
% This function is used to add a result to the next available position in
% the neuronStruct structure, modifying the resultsData field and iterating
% the counters associated with the field.  It will also reallocate space in
% the resultsData matrix as needed.
numNeurons = length(allNeuronsStruct); % The number of structures to add data to
pcaSortCodes = cell(1,numNeurons); % Initialize the sort codes to empty
if nargin > 7
    pcaSortCodes = varargin{1};
end

for neuronNum = 1:numNeurons
    % Check if structure has been initialized
    if ~isNeuronStruct(allNeuronsStruct{neuronNum})
        allNeuronsStruct{neuronNum} = initializeNeuronStruct; % Use the current size to make a properly initialized structure
    end
    
    % Check if resultsData needs to be reallocated
    if (allNeuronsStruct{neuronNum}.resultsDataCurLoc > allNeuronsStruct{neuronNum}.resultsDataAllocated)
        allNeuronsStruct{neuronNum}.resultsDataAllocated = allNeuronsStruct{neuronNum}.resultsDataAllocated*2; % Double the amount of allocated space
        
        % Copy old data
        oldResultsData = allNeuronsStruct{neuronNum}.resultsData;
        oldVoltageTraces = allNeuronsStruct{neuronNum}.voltageTraces;
        oldVoltageTracesRaw = allNeuronsStruct{neuronNum}.voltageTracesRaw;
        oldPCASortCodes = allNeuronsStruct{neuronNum}.pcaSortCodes;
        oldTriggerVoltage = allNeuronsStruct{neuronNum}.triggerVoltage;
        oldTriggerTimespan = allNeuronsStruct{neuronNum}.triggerTimespan;
        
        % Reallocate new memory
        allNeuronsStruct{neuronNum}.resultsData = nan(allNeuronsStruct{neuronNum}.resultsDataAllocated, 4);
        allNeuronsStruct{neuronNum}.voltageTraces = cell(allNeuronsStruct{neuronNum}.resultsDataAllocated,1);
        allNeuronsStruct{neuronNum}.voltageTracesRaw = cell(allNeuronsStruct{neuronNum}.resultsDataAllocated,1);
        allNeuronsStruct{neuronNum}.pcaSortCodes = cell(allNeuronsStruct{neuronNum}.resultsDataAllocated,1);
        allNeuronsStruct{neuronNum}.triggerVoltage = allNeuronsStruct{neuronNum}.defaultTriggerSettings(1)*ones(allNeuronsStruct{neuronNum}.resultsDataAllocated, 1);
        allNeuronsStruct{neuronNum}.triggerTimespan = repmat(allNeuronsStruct{neuronNum}.defaultTriggerSettings(2:3),allNeuronsStruct{neuronNum}.resultsDataAllocated, 1).*ones(allNeuronsStruct{neuronNum}.resultsDataAllocated,2);
        
        % Copy old data back into the new memory
        allNeuronsStruct{neuronNum}.resultsData(1:size(oldResultsData,1), :) = oldResultsData;
        allNeuronsStruct{neuronNum}.voltageTraces(1:size(oldVoltageTraces,1), :) = oldVoltageTraces;
        allNeuronsStruct{neuronNum}.voltageTracesRaw(1:size(oldVoltageTraces,1), :) = oldVoltageTracesRaw;
        allNeuronsStruct{neuronNum}.pcaSortCodes(1:size(oldVoltageTraces,1), :) = oldPCASortCodes;
        allNeuronsStruct{neuronNum}.triggerVoltage(1:size(oldTriggerVoltage,1), :) = oldTriggerVoltage;
        allNeuronsStruct{neuronNum}.triggerTimespan(1:size(oldTriggerTimespan,1), :) = oldTriggerTimespan;
    end
    
    % Add data to resultsData
    allNeuronsStruct{neuronNum}.resultsData(allNeuronsStruct{neuronNum}.resultsDataCurLoc, :) = [strength, duration, result(neuronNum), time]; % Add properly formatted results information
    allNeuronsStruct{neuronNum}.voltageTraces{allNeuronsStruct{neuronNum}.resultsDataCurLoc} = voltageTrace{neuronNum}; % Add the voltage trace data
    allNeuronsStruct{neuronNum}.voltageTracesRaw{allNeuronsStruct{neuronNum}.resultsDataCurLoc} = voltageTraceRaw{neuronNum}; % Add the voltage trace data
    allNeuronsStruct{neuronNum}.pcaSortCodes{allNeuronsStruct{neuronNum}.resultsDataCurLoc} = pcaSortCodes{neuronNum}; % Add the voltage trace data
    allNeuronsStruct{neuronNum}.triggerVoltage(allNeuronsStruct{neuronNum}.resultsDataCurLoc) = allNeuronsStruct{neuronNum}.defaultTriggerSettings(1); % Add the trigger voltage
    allNeuronsStruct{neuronNum}.triggerTimespan(allNeuronsStruct{neuronNum}.resultsDataCurLoc, :) = allNeuronsStruct{neuronNum}.defaultTriggerSettings(2:3); % Add the trigger timespan
    allNeuronsStruct{neuronNum}.resultsDataCurLoc = allNeuronsStruct{neuronNum}.resultsDataCurLoc + 1; % Iterate the current data location index
end

% Output data
newNeuronStruct = allNeuronsStruct;


function [result, resultLocs, validResults] = detectSpikesFromTraces(traceData, traceTime, defaultTriggerSettings, sortCodes, varargin)
% TODO SILICON: Must also include spike indicator for tetrode sorting?  Maybe can group tetrode and PCA sorting markers?%

% This function will detect spikes from the trace in traceData, and report
% back if they occurred, as well as where it occurred.
%
% TraceData is a length = m vector which represents m samples of voltage
% recording from an electrode.
%
% result will be a logical scalar, represeting whether or not the trace
% produced a spike (0 for no spike, 1 for spike). If there is an error in
% the spike detection, the result will be a nan.
%
% resultLocs will be a vector that represents the indices in traceData that
% spikes (threshold crossings or sortcodes) occurred.
%
% validResults will be a vector of the same size as resultLocs that is true
% for spikes that are within the timespan and false for spikes that are
% outside of it.
%
% TraceTime represents how many milliseconds the trace lasts (the timestamp
% of the last sample).
%
% DefaultTriggerSettings is a 1x3 vectors that holds the voltage threshold,
% and the lower and upper limits of the timespan in which a spike may
% occur.
%
% SortCodes is either empty, or a length = m vector that contains the
% sort-codes produced by PCA analysis for each given time step (0's
% otherwise). If sortCodes is an empty matrix, then voltage crossings will
% be used to determine spikes.  If it is a length = m vector (as described
% above), it assumed that these sort-codes are the results of PCA analysis,
% and that they should be used to determine if/when spikes occurred.
%
% An optional input is the pca sort code that the function should be
% looking for. (default: 1)

% Parse inputs
if nargin > 4
    pcaCode = varargin{1};
else
    pcaCode = 1;
end

isTesting = false;
if nargin > 5
    isTesting = varargin{2}; % Is the user testing the program now or no?
end

usePCA = false;
if ~isempty(sortCodes)
    usePCA = true; % Because the sortCodes vector is full, do not use voltage thresholding
end

% TEST
% figure;
% plot(sortCodes); % TEST

spikeThreshold = defaultTriggerSettings(1); % Get the spike threshold
spikeTimespan = defaultTriggerSettings(2:3); % Get the spike timespan

% Prepare parameters
if isTesting
    pcaCodeDelay = 0;
    %     vThreshLower = spikeThreshold*.7;
    vThreshLower = spikeThreshold*.95;
else
    pcaCodeDelay = 3.2; %5.8; % The delay in ms between a spike's peak and the sort-code being produced (found via experimentation)
end

% Create time axis for this trace, and find those inds within the
% spikeTimespan
times = linspace(0, traceTime, length(traceData));
timeGranularity = diff(times([1 2])); % Time between two points
timesWithinTimespan = (spikeTimespan(1) < times) & (times < spikeTimespan(2)); % Indices of voltage trace samples within the timespan

% if isTesting
%     [maxData, resultLocs] = max(traceData); % Get the max and location of the max
%     %     validResults = timesWithinTimespan(resultLocs);
%     %     PFiring = max(min((maxData - vThreshLower)/(spikeThreshold - vThreshLower),1),0); % Get the probability of firing, bounded between 0 and 1
%     %     result = PFiring > rand; % Determine if a random variable is larger than the probability of spiking to estimate spiking liklihood
%     result = maxData > vThreshLower;
% else
if usePCA
    sortCodeTimeSpan = spikeTimespan + pcaCodeDelay; % The pca code is delayed relative to the spike, so the timespan is translated
    sortCodeTimesWithinTimespan = (sortCodeTimeSpan(1) < times) & (times < sortCodeTimeSpan(2)); % Indices of sort-code samples within the timespan
    anyCodes = sortCodes == pcaCode; % Get locations of any places that sortCodes is equal to 1
    codesWithinTimespan = sortCodes(sortCodeTimesWithinTimespan) == pcaCode; % Get only codes that are equal to 1 (there could be other artifacts that are sorted into different codes)
    result = any(codesWithinTimespan); % If any of the sort-codes are non zero, then there was a result
    resultLocs = find(anyCodes) - round(pcaCodeDelay/timeGranularity); % Where is the first index of a sort code? (taking into account the pca code delay)
else
    % Determine if a threshold crossings occurred
    if (spikeThreshold > 0)
        % For positive spike thresholds, the voltage must go above it to spike
        hitsWithinTimespan = traceData(timesWithinTimespan) > spikeThreshold;
        resultLocs = find(traceData > spikeThreshold); % Where are the threshold crossings?
    else
        % For negative spike thresholds, the voltage must go below it to spike
        hitsWithinTimespan = traceData(timesWithinTimespan) < spikeThreshold;
        resultLocs = find(traceData < spikeThreshold); % Where are the threshold crossings?
    end
    
    result = any(hitsWithinTimespan); % Check if any threshold crossings occurred
end
% end

% Finalize resultLocs and calculate validResults
if ~isempty(resultLocs)
    outOfBoundsInds = (1 > resultLocs) | (resultLocs > length(traceData));
    resultLocs(outOfBoundsInds) = []; % Delete any references to spikes that are not within the trace (which may show up due to lags between the trace data and the pca code)
    
    % Determine which resultLocs fall within the timespan
    %     potentiallyValidResultLocs = resultLocs(~outOfBoundsInds); % Only see if in-bounds resultLocs are valid
    allSpikeLocs = false(size(times));
    allSpikeLocs(resultLocs) = true; % Create a logical array the size of times, with true's at the location of each spike
    spikeLocsWithinTimespan = allSpikeLocs & timesWithinTimespan; % True for spikes that are within a timespan
    validResults = spikeLocsWithinTimespan(resultLocs); % For each result location, return true if it is within the timespan, and false if it is outside of it
else
    validResults = [];
end

%% Additional GUI Functions

function handles = checkNewInput(hObject, handles, varDataName, varargin)
% CHECKED %
% This function reads the current values of the GUI object pointed to by
% handle (namely editboxes), and ensures that it is a real, finite value
% between the values allowed by varData.
% hObject is the handle to the calling object, handles is the GUI's handle
% data, and varDataName is a string that corresponds to the name of the
% handles field that represents this GUI object's varData
% An optional 4th input allows the user to constrain the input to integers

% See if user wants to constrain input to an integer value
isInteger = false;
if nargin > 3
    isInteger = logical(varargin{1});
end

% Get new values from GUI
newVal = str2double(get(hObject,'String'));

% Get the lower bound, upper bound, and current value for the GUI object
varData = handles.(varDataName);

% Check if they are finite, and if so, constrain them to the min and max
% defined in the varData
% Check strength
newVal = returnValidValue(newVal, varData, isInteger);

% Assign this new (now valid) value back to the varData vector, to be
% recorded as the most recent valid entries into the GUI
varData(1) = newVal; % Update the local version of varData
handles.(varDataName) = varData; % Update the handles version of varData
guidata(hObject, handles) % Update handles (so that the varData change takes effect)

% After checking the values, reassign these new values to the GUI
set(hObject, 'String', num2str(newVal));

% Check the default values of all important parameters of the trial
checkDefaults(handles);

function newVal = returnValidValue(oldVal, varData, varargin)
% CHECKED %

% This function compares oldVal to the varData vector, and returns a newVal
% that is either equal to the oldVal, or is equal to varData(1) if oldVal
% does not satisfy the conditions (finite, real, between bounds set by
% varData(2) and varData(3). Used when checking GUI input data, to ensure
% that the new values are valid
% Optional third parameter allows user to constrain the input to an integer
% (if true)

% See if user wants to constrain input to an integer value
isInteger = false;
if nargin > 2
    isInteger = logical(varargin{1});
end

if isfinite(oldVal) && isreal(oldVal) && (oldVal >= varData(2)) && (oldVal <= varData(3)) && ((mod(oldVal, 1) == 0) || ~isInteger)
    % If conditions are satisfied, use the original value sent in
    newVal = oldVal;
else
    % If the value is not finite, real, and between the bounds, throw out
    % the new value and replace it with the most recent valid value
    newVal = varData(1);
    
end
return

function neuronStruct = initializeNeuronStruct
% CHECKED
% Initialize all variables
resultsDataAllocated = 100; % Number of slots to preallocate for results data for this neuron
numParameters = 3; % Tracking 3 parameters (alpha, beta, and w)

% Use variables to initialize other parameters
resultsDataCurLoc = 1; % Current location along the resultsData matrix
resultsData = nan(resultsDataAllocated, 4); % TODO SILICON: For each stimulation, list which channel in the array is matched to this neuron % Preallocate data for resultsData, with 4 columns to track strength, duration, result, and time
parameters = []; % Parameters array to track the fitted alpha, beta, and w for the neuron
optS = []; % Structure that represents the metadata around optimization (i.e. the calculation time, stimulation indices used, and the thetas calculated)
laserEQ = []; % The parameters of the laser equation in TDT (of the form [A B C D], where power = polyval(laserEQ, laserControlStrength))
% characterizedDurs = []; % Durations that have been characterized for this neuron
voltageTraces = cell(resultsDataAllocated,1); % Preallocate a cell array to hold onto the filtered voltage trace data from each stimulation
voltageTracesRaw = cell(resultsDataAllocated,1); % Preallocate a cell array to hold onto the raw (unfiltered) voltage trace data from each stimulation
pcaSortCodes = cell(resultsDataAllocated,1); % TODO SILICON: May need new sorting parameters after tetrode sorting is introduced* Preallocate a cell array to hold the sort codes for this neuron for each stimulation
defaultTriggerSettings = [.5, 0, 30]; % The default settings for the triggerVoltage, and the triggerTimespan lower and upper bounds
% SDCurvePoints = []; % The strength-duration curve data that will be plot (strengths are the best fit values for a 50% probability of firing)
triggerVoltage = nan(resultsDataAllocated, 1); % The trigger voltage for each stimulation to count as a spike
triggerTimespan = nan(resultsDataAllocated,2); % The lower and upper bounds for the timespan in which a spike must occur
sequenceInfo = []; % This matrix holds all of the meta data associated with a running a sequence.  Column 1 tracks the index of the first stimulation in the sequence (the first slot), columns 2-[numSlots + 1] track if the neuron was supposed to fire (1) or not (0) during each slot, and the last two columns record the strength and duration (in that order) that was used to index this neuron during the given sequence
sequenceInfoCurLoc = 1; % The next index for new sequence info
trialData = struct('numRunsPerTrial', [], 'numStimsPerSequence', [], 'numSequencesPerRun', [], 'expectedResults', [], 'interStimPeriod', [], 'runOrders', [], 'stimTimes', [], 'stimInds', [], 'opt', [], 'blockStartInfo', [], 'TDTPenetNum', [], 'trialCompleted', false); % Keeps track of information needed for the main experiment
% adaptiveInfo = struct('epochStartSeqs', [], 'searchPoints', [], 'verifyPoints', [], 'escapePoints', []); % Keeps track of the search- and escapePoints for each each epoch of all experiments, the verification points for each experiment, the id's of each of the search-/escape-/verify- points in that epoch's meshPointLocations, the starting sequence number for each epoch, and the mesh point ID's in each epoch
% adaptiveInfo = struct('meshPointLocations', [], 'epochStartSeqs', [], 'searchPoints', [], 'verifyPoints', [], 'escapePoints', [], 'searchCodes', [], 'verifyCodes', [], 'escapeCodes', []); % Keeps track of the search- and escapePoints for each each epoch of all experiments, the verification points for each experiment, the id's of each of the search-/escape-/verify- points in that epoch's meshPointLocations, the starting sequence number for each epoch, and the mesh point ID's in each epoch

% Create a structure initialized with the correct fields and starting values
neuronStruct.resultsDataAllocated = resultsDataAllocated;
neuronStruct.numParameters = numParameters;
neuronStruct.resultsDataCurLoc = resultsDataCurLoc;
neuronStruct.resultsData = resultsData;
% neuronStruct.SDCurvePoints = SDCurvePoints;
neuronStruct.parameters = parameters;
neuronStruct.optS = optS;
neuronStruct.laserEQ = laserEQ;
% neuronStruct.characterizedDurs = characterizedDurs;
neuronStruct.voltageTraces = voltageTraces;
neuronStruct.voltageTracesRaw = voltageTracesRaw;
neuronStruct.pcaSortCodes = pcaSortCodes;
neuronStruct.defaultTriggerSettings = defaultTriggerSettings;
neuronStruct.triggerVoltage = triggerVoltage;
neuronStruct.triggerTimespan = triggerTimespan;
neuronStruct.sequenceInfo = sequenceInfo;
neuronStruct.sequenceInfoCurLoc = sequenceInfoCurLoc;
neuronStruct.trialData = trialData;
neuronStruct.useIrradiance = true; % The script now uses irradiance control
% neuronStruct.adaptiveInfo = adaptiveInfo;

function allNeuronStructs = initializeAllNeuronStructs(numNeurons)
% CHECKED %
% This function is used to initialize the allNeuronStructs structure, which
% holds results information on each neuron.  It is assumed that the index
% of each neuron corresponds to the same numbered physical electrode
% connected to the TDT system (i.e. allNeuronStructs{1} holds information
% on the neuron being recorded by Electrode_1).

% Create and populate final cell array
allNeuronStructs = cell(numNeurons, 1); % Initialize the cell array to hold the structs

% For each one of the neurons in this program, create a newly initialized
% neuronStruct
for neuronInd = 1:numNeurons
    allNeuronStructs{neuronInd} = initializeNeuronStruct;
end

function handles = checkNewVectorInput(hObject, handles, curValsName, boundsName)
% CHECKED %

% Function called by the durations edit boxes upon editing.  Ensures that
% the input values are properly formatted and within range.

newVals = [];
try
    eval(['newVals = ' get(hObject, 'String') ';']); % Attempt to assign the new values into a varaible, newVals
catch e
    % If assignment didn't work, replace the object's string with its last
    % known value
    disp(e)
    newVals = handles.(curValsName);
end

% Check every value in the string
for i = 1:length(newVals)
    % Make sure that each new value is within bounds.  If not, replace it
    % with a NaN
    newVals(i) = returnValidValue(newVals(i), [NaN handles.(boundsName)]);
end

% Set all NaN values to be empty
newVals(isnan(newVals)) = [];
if isempty(newVals)
    % If newVals is now empty, because all values were out of bounds (set
    % to NaN's), replace the entire vector with the last valid vector
    newVals = handles.(curValsName);
end

% Sort the newVals vector, in ascending order
newVals = sort(newVals);

% Set newVals to be the new durations
handles.(curValsName) = newVals; % Update the last known valid value
set(hObject, 'String', vec2String(newVals)); % Update the GUI with valid values
guidata(hObject, handles);

%% Additional Hardware Functions
function handles = startOptServeComm(handles)
% This function will attempt to connect to the optimization server

% If the optimization server is already connected, then don't bother trying
% to reconnect
handles = checkOptServeConnectivity(handles);
if handles.optServeConnected
    return;
end

% Confirm that the user has turned on the server (otherwise this will block
% for 10 seconds or so, leading to much annoyance)
choice = questdlg('Has the optimization server been started?','Server ready?','Yes','Cancel','Cancel');
if strcmp(choice, 'Cancel')
   return; 
end

% Open optimization server connection
try
   fopen(handles.optServeT);
end

% Check the connection and update the GUI
handles = checkOptServeConnectivity(handles);

% If the communication channel is open, send the irradiation ratio
if handles.optServeConnected
    tcpMatWrite(handles.optServeT, handles.irrRatio);
end

function handles = startTDTComm(handles)
% CHECKED %
% This function is used to start communication with the TDT system. It will
% create an active x control object, and initialize it to the default
% values needed to communicate with the TDT system.  The connection can be
% closed by calling endTDTComm on the DA object created by this function.

%% TDT Communication Notes
% Read variable from TDT
% matlabTargetVar = DA.GetTargetVal('RZ2.TDT_Var_To_Read');
%
% Read buffer from TDT
% matlabTargetVar = DA.ReadTargetVEX('RZ2.TDT_Buffer_To_Read',offsetFromBeginning,numWordsToReceive,sourceTypeAsString,destinationTypeAsString);
% Ex: Electrode1_Spikes = DA.ReadTargetVEX('RZ2.Electrode1_Spikes',0,Index,'I32','F32');
%
% Write variable to TDT
% DA.SetTargetVal('RZ2.TDT_Var_To_Write',matlabSourceVar);
%
% Put Into idle mode
% DA.SetSysMode(0)

%% Start communication with TDT
% 07/28/2015
% MP
try
    % If this is a Windows computer
    DA = actxcontrol('TDevAcc.X');
    %     startTime = clock;
    %     while etime(clock, startTime) < 5
    DA.ConnectServer('Local');
    %         handles = checkTDTConnectivity(handles);
    %         if handles.TDTConnected
    %            break;
    %         end
    %     end
    
    % Turn TDT into Matlab State
    DA.SetTargetVal('RZ2.Matlab State On',1);
    
    % Check if SpikePAC is being used on the TDT server
    handles.TDTUsePCA = DA.GetTargetVal('RZ2.usePCA');
    
    % TODO SILICON: Check also for tetrode sorting, possibly check for silicon
    % probe vs carbostar %
    
    handles.DA = DA; % Store DA in handles
end

% Check the TDT connection
handles = checkTDTConnectivity(handles);

% Try to extract the laser equation parameters from TDT
handles = setMaxStrengthFromIrradiance(handles);
handles = setLaserParameters(handles);
handles = setControlVoltages(handles);

function handles = checkTDTConnectivity(handles)
% CHECKED %

% This function will check the TDT connection, and update the GUI to
% represent the current TDT connection state

% Check TDT connection
connected = false;
try
    if handles.DA.CheckServerConnection
        % If DA valid (no error), and it is connected to TDT
        % (sending/receiving information), then list the program as
        % connected.
        connected = true;
    end
catch e
    disp(e);
end

% Update handles' version of the connection status
handles.TDTConnected = connected;

handles = updateGUI(handles); % Update the GUI with respect to the new connected status
checkDefaults(handles); % Check if the default parameters are currently in place


function handles = checkOptServeConnectivity(handles)
% This function will check the optimization server connection, and update
% the GUI to represent the current optimization server connection state

handles.optServeConnected = strcmp(handles.optServeT.Status, 'open');

handles = updateGUI(handles); % Update the GUI with respect to the new connected status
checkDefaults(handles); % Check if the default parameters are currently in place

function isRecording = checkTDTRecording(handles)
% CHECKED %

% This function will check if TDT is currently set to a recording state
if ~isempty(handles.DA)
    isRecording = handles.DA.GetSysMode == 3;
else
    isRecording = false;
end

function handles = updateGUI(handles)
% TODO: Any new GUI elements should probably be updated in this routine %
% Update the GUI (colors, enabled functions, etc)

% Toggle the enable status of the "stimulation" (either real or fake)
% buttons
if  ~(handles.TDTConnected || get(handles.testing, 'Value'))
    % If the program is not connected to TDT, and the user is not testing
    
    % Change enabled status of all stimulation buttons
    handles = enableAllStim(handles, false);
else
    % If the program is in testing mode or is connected to TDT
    
    % Change enabled status of all stimulation buttons
    handles = enableAllStim(handles, true);
end

% Determine if sequence start button should be enabled
% if (handles.TDTConnected || get(handles.testing, 'Value')) && (size(handles.allNeuronsStruct{1}.sequenceInfo,1) == size(handles.allNeuronsStruct{2}.sequenceInfo,1))
%     set(handles.controlSequenceStart, 'Enable', 'on');
% else
%     set(handles.controlSequenceStart, 'Enable', 'off');
% end

% Toggle things related to the actual TDT connectivity
if handles.TDTConnected
    % If the system is currently connected to the TDT platform
    
    % Change Connect/Disconnect button
    set(handles.TDTConnectButton, 'String', 'Disconnect');
    
    % See if the TDT program is set to recording or not
    if checkTDTRecording(handles)
        % Change recording button
        set(handles.setToRecording, 'Enable', 'off');
        
        % Change connection flag
        set(handles.connectedFlag, 'ForegroundColor', handles.brightGreen);
        set(handles.connectedFlag, 'String' , 'Connected and Recording');
    else
        % Change recording button
        set(handles.setToRecording, 'Enable', 'on');
        
        % Change connection flag
        set(handles.connectedFlag, 'ForegroundColor', handles.yellow);
        set(handles.connectedFlag, 'String' , 'Connected, Not Recording');
    end
    
else
    % If the program is not connected to TDT
    
    % Change flag
    set(handles.connectedFlag, 'ForegroundColor', handles.red);
    set(handles.connectedFlag, 'String' , 'Not Connected');
    
    % Change Connect/Disconnect button
    set(handles.TDTConnectButton, 'String', 'Connect');
    
    % Change recording button
    set(handles.setToRecording, 'Enable', 'off');
end


if handles.optServeConnected
    % Change the connect/disconnect button text
    set(handles.optServeConnect, 'String', 'Disconnect');
    
    % Set the connection status banner
    set(handles.optServeStatus, 'String', 'Connected', 'ForegroundColor', handles.green);
else
    % Change the connect/disconnect button text
    set(handles.optServeConnect, 'String', 'Connect');
    
    % Set the connection status banner
    set(handles.optServeStatus, 'String', 'Not Connected', 'ForegroundColor', handles.red);
end

% Always reset the fitSDCurve buttons, listboxes, and connect button to on
buttonHandles = [handles.fitSDCurve handles.neuron1StimsListbox handles.neuron2StimsListbox handles.TDTConnectButton handles.optServeConnect];
set(buttonHandles, 'Enable', 'on');

% Always set the trial displayed to be off if running a trial
if handles.runningTrial
    str = 'on';
else
    str = 'off';
end
displayHandles = [handles.runNumberText, handles.runNumberOut, handles.sequenceNumberText, handles.sequenceNumberOut, handles.elapsedTimeText, handles.elapsedTimeOut, handles.runningEstTimeText, handles.runningEstTimeOut, handles.stopTrial, handles.pauseTrial];
set(displayHandles, 'visible', str);

drawnow;

function handles = enableAllStim(handles, doEnable)
% CHECKED %

% This is a shortcut function to setting the enabled status of all of the
% buttons that will stimulate

if doEnable
    str = 'on';
    handles.stimEnabled = true;
else
    str = 'off';
    handles.stimEnabled = false;
end

buttonHandles = [handles.characterizeNeurons handles.manStim handles.manStim1 handles.manStim2 handles.trialStart handles.swapGTs];

% Change enabled status of all stimulation buttons
set(buttonHandles, 'Enable', str);

% Warn the user if the current csf file has not been modified or if the
% program is not currently set to recording
if get(handles.testing, 'Value') || (handles.haveControllerFilePath(handles.modCSFIndex) && handles.controllerFileModified(handles.modCSFIndex) && checkTDTRecording(handles))
    set(handles.trialStart, 'BackgroundColor', handles.gray);
else
    set(handles.trialStart, 'BackgroundColor', handles.yellow);
end

% --- Executes on selection change in neuron2StimsListbox.
function neuron2StimsListbox_Callback(hObject, eventdata, handles)
% hObject    handle to neuron2StimsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns neuron2StimsListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from neuron2StimsListbox
handles = stimsListboxCallback(handles, get(hObject, 'Value'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function neuron2StimsListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neuron2StimsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in applyTrigWindowSelected2.
function applyTrigWindowSelected2_Callback(hObject, eventdata, handles)
% hObject    handle to applyTrigWindowSelected2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

neuronNum = get(hObject, 'UserData'); % The number of the selected neuron
applyTrigParamsCallback(hObject, handles, handles.selectedStims{neuronNum}) % Apply the parameters to only the selected neurons

% --- Executes on button press in setTriggerVoltage2.
function setTriggerVoltage2_Callback(hObject, eventdata, handles)
% hObject    handle to setTriggerVoltage2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = setTriggerVoltage(hObject, handles); % Call the routine to set the trigger voltage, for this neuron
guidata(hObject, handles); % Update handles

% --- Executes on button press in setTriggerTimespan2.
function setTriggerTimespan2_Callback(hObject, eventdata, handles)
% hObject    handle to setTriggerTimespan2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = setTriggerTimespan(hObject, handles);
guidata(hObject, handles);

% --- Executes on button press in applyTrigWindowAll2.
function applyTrigWindowAll2_Callback(hObject, eventdata, handles)
% hObject    handle to applyTrigWindowAll2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

neuronNum = get(hObject, 'UserData'); % The number of the selected neuron
numStims = handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc - 1; % The number of stimulations that are currently being tracked for this neuron
applyTrigParamsCallback(hObject, handles, 1:numStims) % Apply the parameters to only the selected neurons

function strength2Min_Callback(hObject, eventdata, handles)
% hObject    handle to strength2Min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strength2Min as text
%        str2double(get(hObject,'String')) returns contents of strength2Min as a double
handles = checkNewInput(hObject, handles, 'strength2MinVarData');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function strength2Min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strength2Min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function strength2Max_Callback(hObject, eventdata, handles)
% hObject    handle to strength2Max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strength2Max as text
%        str2double(get(hObject,'String')) returns contents of strength2Max as a double
handles = checkNewInput(hObject, handles, 'strength2MaxVarData');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function strength2Max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strength2Max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function durations2_Callback(hObject, eventdata, handles)
% hObject    handle to durations2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of durations2 as text
%        str2double(get(hObject,'String')) returns contents of durations2 as a double

handles = checkNewVectorInput(hObject, handles, 'durations2CurVals', 'durations2VarDataBounds');

% Record the new values
neuronNum = get(hObject, 'UserData');
handles.characterizationDurations{neuronNum} = handles.durations2CurVals;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function durations2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to durations2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha2_Callback(hObject, eventdata, handles)
% hObject    handle to alpha2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha2 as text
%        str2double(get(hObject,'String')) returns contents of alpha2 as a double

% --- Executes during object creation, after setting all properties.
function alpha2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta2_Callback(hObject, eventdata, handles)
% hObject    handle to beta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta2 as text
%        str2double(get(hObject,'String')) returns contents of beta2 as a double

% --- Executes during object creation, after setting all properties.
function beta2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w2_Callback(hObject, eventdata, handles)
% hObject    handle to w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w2 as text
%        str2double(get(hObject,'String')) returns contents of w2 as a double

% --- Executes during object creation, after setting all properties.
function w2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = loadNeuronStruct(handles, neuronStructIn)
% CHECKED %
% This function will load in the inputted allNeuronsStruct, by putting it
% into the handles structure as well as modifying the GUI as needed.  It
% assumes that the inputted allNeuronsStruct is an initialized (not
% necessarily filled) allNeuronsStruct object

% First, add it to the handles struct
handles.allNeuronsStruct = neuronStructIn; % Where the struct will be referenced later
handles.lastSavedNeuronsStruct = neuronStructIn; % Indicating that this struct does not need to be saved as is (because it is either loaded from a file, or just initialized)

% Load the appropriate GUI information for this neuron

% Load all of the stimulations, initialize the currently selected
% voltage trace, and redraw the plot
handles = loadStimulationsIntoGUI(handles);

for neuronNum = 1:handles.numNeurons
    % Load this struct's parameters into the GUI, if applicable
    
    % Get the most recently calculated parameters
    if ~isempty(neuronStructIn{neuronNum}.parameters)
        aNew = neuronStructIn{neuronNum}.parameters(end, 1);
        bNew = neuronStructIn{neuronNum}.parameters(end, 2);
        wNew = neuronStructIn{neuronNum}.parameters(end, 3);
    else
        aNew = '';
        bNew = '';
        wNew = '';
    end
    set(handles.(['alphaOut' num2str(neuronNum)]), 'String', aNew);
    set(handles.(['betaOut' num2str(neuronNum)]), 'String', bNew);
    set(handles.(['wOut' num2str(neuronNum)]), 'String', wNew);
    
    %     % Clear the example spike data
    %     handles.exampleSpikeCollection{neuronNum} = [];
end

% Ensure that the GUI does not try to load a sequence from this neuron
% handles.hasFiredSequence = false;
% updateGUISequenceIndicators(handles);
handles = updateGUI(handles); % Ensure that a new sequence will only be run if the two neurons' sequenceInfo are synchronized

% Clear any previous autosave data for this electrode in preparation for
% the new neuron
handles = clearAutoSave(handles);

% If this neuron has been characterized, plot it (the check occurs
% within the function)
handles = drawSDPlot(handles);

% Record that the struct has been loaded for this neuron
handles.structsLoaded = true(1, handles.numNeurons);

function handles = clearAutoSave(handles)
% CHECKED %

% Clears the auto save data
handles.autoSaveFiles = '';
set(handles.autoSaveStatus, 'BackgroundColor', 'red', 'String', 'No Save');

function handles = loadStimulationsIntoGUI(handles)
% TODO SILICON: Include selection for tetrode channel? %
% This function loads all of the stimulations associated with neuron number
% neuronNum into the GUI (specifically, into neuronXStimsListbox).

value = zeros(1,handles.numNeurons);
listBoxStrings = cell(1, handles.numNeurons);
for neuronNum = 1:handles.numNeurons
    % For each neuron, load new stimulations into the listbox, and from
    % there into the voltage trace plot
    if handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc == 1
        % If the current location is 1, meaning that no stimulations have been
        % recorded, display the associated message in the listbox
        listBoxStrings{neuronNum} = 'No Stims';
        handles.stimsExist(neuronNum) = false; % No stimulations exist
        value(neuronNum) = 1; % The first (and only) item must be selected
    else
        % If there are stimulations, load the information for each one into
        % the box
        numStims = handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc - 1; % The number of stimulations is equal to the current location minus 1
        listBoxStrings{neuronNum} = cell(numStims,1);
        for i = 1:numStims
            listBoxStrings{neuronNum}{i} = sprintf('s=%0.1f, d=%0.2f', handles.allNeuronsStruct{neuronNum}.resultsData(i, 1), handles.allNeuronsStruct{neuronNum}.resultsData(i, 2));
        end
        handles.stimsExist(neuronNum) = true; % Stimulations exist
        value(neuronNum) = numStims; % Make sure the last stimulation is selected so that it is displayed
    end
end

for neuronNum = 1:handles.numNeurons
    % For each neuron, load new stimulations into the listbox, and from
    % there into the voltage trace plot
    listboxHandle = handles.(['neuron' num2str(neuronNum) 'StimsListbox']); % Get the handle associated with this neuron's listBox
    
    set(listboxHandle, 'String', listBoxStrings{neuronNum}, 'Value', value(neuronNum)); % Set the strings and choice for the list box
end

% Call the listbox's callback function, to load in the trigger
% voltage/timespan, as well as draw the voltage trace
handles = stimsListboxCallback(handles); % As if the new stimulation was selected

% --- Executes on button press in saveNeuron.
function saveNeuron_Callback(hObject, eventdata, handles)
% hObject    handle to saveNeuron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = saveNeuronStruct(handles);
guidata(hObject, handles);

function handles = saveNeuronStruct(handles, varargin)
% CHECKED %
% Save the incoming neuron struct to both a .mat file and as a variable in
% the workspace

createFileName = true;
if nargin > 2
    createFileName = false;
    varFileName = varargin{1};
end

neuronStruct = handles.allNeuronsStruct; % Get the struct that will be saved

% Check if this struct to be saved is the same as the last version of the
% struct that was saved. If it is the same, do not bother saving it.
if (~isEqualNeuronStruct(handles.lastSavedNeuronsStruct, neuronStruct))
    % It is different, so save this struct
    % Update the status
    set(handles.autoSaveStatus, 'BackgroundColor', 'yellow', 'String', 'Saving...');
    drawnow;
    
    % Save the inputted neuron struct to the lastSavedNeuronStruct array,
    % so that it can be compared to later
    handles.lastSavedNeuronsStruct = neuronStruct;
    
    if createFileName
        % If there was no filename that was inputted
        if get(handles.autoSave, 'Value') && ~isempty(handles.autoSaveFiles)
            % Save using the existing autosave file name
            varFileName = handles.autoSaveFiles;
        else
            % Check for a name in the current directory is that is not currently used
            baseName = handles.baseFileName; % The base name for the variable in the base workspace
            curClock = clock;
            %             dateTimeString = sprintf('_%0.2d-%0.2d-%0.2d_%0.2d-%0.2d', curClock(2), curClock(3), curClock(1) - 2000, curClock(4), curClock(5));
            dateTimeString = sprintf('%d%.2d%.2d_%.2d%.2d', curClock(1), curClock(2), curClock(3), curClock(4), curClock(5));
            varFileName = [pwd filesep baseName dateTimeString '.mat']; % The file name that this structure will be stored in
            varNameModifier = 2; % A counter to add to the base name to differentiate different neurons
            % Check if the date-time version already exists as a file
            while exist(varFileName, 'file') % Check the name
                % If so, iterate through names by adding a variable modifier to
                % the end
                varFileName = [pwd filesep baseName dateTimeString sprintf('_%0.2d', varNameModifier) '.mat']; % Modify the file name
                varNameModifier = varNameModifier + 1; % Iterate the modifer
            end
        end
    end
    
    % Make sure that the variable name is the same as handles.structVarName
    eval([handles.structVarName ' = neuronStruct;']);
    
    % Save data to file
    % Data is saved in case of a Matlab error
    save(varFileName, handles.structVarName);
    
    % Save the name of the file
    handles.autoSaveFiles = varFileName;
    
    % Update the status
    set(handles.autoSaveStatus, 'BackgroundColor', 'green', 'String', 'Saved');
    drawnow;
end

function equal = isEqualNeuronStruct(A, B)
% CHEKED %
% This function will compare two neuronStructs, A and B, and output if they
% are equal to each other.
equal = false; % Are the two structs equal to each other?

% Check that both structs are in fact neuronStructs
if ~(isNeuronStruct(A) && isNeuronStruct(B))
    return; % If they are not both neuron structs, then they cannot be equal (or are both not neuron structs, which is outside of this function's scope)
end

% Go through each field, and make sure that the contents are equal
allFields = fieldnames(A); % Both should have the same fieldnames
for i = 1:length(allFields)
    if ~isequalwithequalnans(A.(allFields{i}),  B.(allFields{i}))
        return; % If one of the matching fields between the two structs containts content that is not equal to the other, return from the function with equal = false
    end
end

% If the function has not ended before here, then the two must be equal
% neuronstructs
equal = true;
return



% --- Executes on button press in clearData.
function clearData_Callback(hObject, eventdata, handles)
% hObject    handle to clearData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clearDataCallback(hObject, handles);

% --- Executes on button press in clearData2.

function clearDataCallback(hObject, handles)
% CHECKED %

% Check that the user doesn't delete any information that they wanted,
% then clear the current data from the given neuron

isBlank = false(1,handles.numNeurons);
for neuronNum = 1:handles.numNeurons
    % If all of the neuronStructs are blank, then don't try to clear (will
    % just be a waste of time)
    isBlank(neuronNum) = handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc == 1;
end

if all(isBlank)
    % All structs are blank, so the gui doesn't need to be cleared
    return;
end

% Confirm clearing
choice = questdlg('Are you sure you want to clear all data?','Clear Data?','Yes','No','Yes');
if strcmp(choice, 'No')
    % If the user does not want to clear the data, exit the function
    return;
end

% Check if the structs have been saved
for neuronNum = 1:handles.numNeurons
    notSaved = false;
    % Check if this data has been saved
    if ~(isEqualNeuronStruct(handles.lastSavedNeuronsStruct{neuronNum}, handles.allNeuronsStruct{neuronNum}))
        notSaved = true;
        break;
    end
end

% If they have not, offer to save them
if notSaved
    % Check if the user wants to save the data before clearing it
    choice = questdlg('This data has not been saved.  Save data?  ("Cancel" to cancel data deletion.)','Unsaved Changes','Yes','No','Cancel','Yes');
    switch choice
        case 'Yes'
            % Save the struct first
            handles = saveNeuronStruct(handles); % Handles is returned, so that the lastSavedNeuronStruct will be maintained (it wouldn't be otherwise, because it is overwritten later in this function)
        case 'Cancel'
            % Do not clear the neuron data
            return;
    end
end

% Clear the neuronStruct by loading in a blank structure, carrying over the
% trigger information
newNeuronStruct = cell(1, handles.numNeurons);
for neuronNum = 1:handles.numNeurons
    newNeuronStruct{neuronNum} = initializeNeuronStruct; % Create a new, blank neuron struct
    newNeuronStruct{neuronNum}.defaultTriggerSettings = handles.GUITriggerSettings(neuronNum,:); % Load in the GUI's current trigger settings
end

handles = loadNeuronStruct(handles, newNeuronStruct);

% Reset the current trial number
handles.curTrialNum = 1;

% Save the data
guidata(hObject, handles);

function voltageTrace1_CreateFcn(hObject, eventdata, handles)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% CHECKED %

% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

if ~isempty(handles)
    % Check if neuron data has been saved
    ready = false;
    while ~ready
        % Check if all of the data has been saved
        notSaved = false;
        for neuronNum = 1:handles.numNeurons
            if ~(isEqualNeuronStruct(handles.lastSavedNeuronsStruct{neuronNum}, handles.allNeuronsStruct{neuronNum})) && (handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc ~= 1)
                % The current data has not been saved
                notSaved = true;
                break;
            end
        end
        
        if notSaved
            % If the data has not been saved, then ask if the user wants to save it
            % Check if the user wants to save the data before clearing it
            choice = questdlg('Neuron data has not been saved.  Save data?  ("Cancel" to cancel exiting.)','Unsaved Changes','Yes','No','Cancel','Yes');
            switch choice
                case 'Yes'
                    % Save the struct first
                    handles = saveNeuronStruct(handles); % Handles is returned, so that the lastSavedNeuronStruct will be maintained (it wouldn't be otherwise, because it is overwritten later in this function)
                    
                    % Save changes to handles
                    guidata(hObject, handles);
                    ready = true;
                case 'Cancel'
                    % Do not exit the code
                    return;
                case 'No'
                    % Delete the struct, but confirm first
                    choice = questdlg('Are you sure? ("Yes" to NOT save data)','Unsaved Changes','Yes','No','No');
                    switch choice
                        case 'Yes'
                            % Delete the struct
                            ready = true;
                    end
                    
                    % If the user answers no here, they will be asked again
                    % (hence the while-loop)
            end
        else
            % If the data has been saved, then the figure is ready for
            % closing
            ready = true;
        end
    end
    
    % Disconnect from TDT
    if (handles.TDTConnected)
        handles = endTDTComm(handles); % End communication with TDT
    end
    
    % Disconnect from the optimization server
    if (handles.optServeConnected)
       handles = endOptServeComm(handles); % End communication with the optimization server 
    end
end

% Update handles
guidata(hObject, handles);

% Close the figure
delete(hObject);

function handles = endTDTComm(handles)
% CHECKED %
% This function will end communication with the TDT software

if ispc
    % Reset TDT's Matlab state
    handles.DA.SetTargetVal('RZ2.Mat Reset',1);
    handles.DA.SetTargetVal('RZ2.Mat Reset',0);
    
    % Set TDT's Matlab state to off
    handles.DA.SetTargetVal('RZ2.Matlab State On',0);
    
    % End TDT communication
    handles.DA.CloseConnection;
end

% Make sure that this program knows that TDT is not recording anymore
handles.TDTRecording = false;

function handles = endOptServeComm(handles)
% This function will end the connection to the optimization server

% Check the connection
handles = checkOptServeConnectivity(handles);

if handles.optServeConnected
    % First, send a "-1" to tell the server that we are disconnecting
    % (because otherwise, the damn thing wouldn't even know!)
    tcpMatWrite(handles.optServeT, -1);
    
    % Close the connection
    fclose(handles.optServeT);
end

% Check the connection and update the GUI
handles = checkOptServeConnectivity(handles);

% --- Executes on button press in clearTriggerParametersSelection1.
function clearTriggerParametersSelection1_Callback(hObject, eventdata, handles)
% hObject    handle to clearTriggerParametersSelection1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clearParametersSelection(hObject, handles);


function clearParametersSelection(hObject, handles)
% CHECKED %
% This function clears the selected trigger voltage and trigger timespan
% from handles, as well as clearing it from the GUI

neuronNum = get(hObject, 'UserData');

% Clear the handles object
handles = clearSelectedTriggerParams(handles, neuronNum);

% Update the handles object
guidata(hObject, handles);

% Redraw the voltage plot (now without the selected trigger voltage and
% timespan)
handles = drawVoltageTraces(handles, neuronNum);
guidata(hObject, handles);


% --- Executes on button press in clearTriggerParametersSelection2.
function clearTriggerParametersSelection2_Callback(hObject, eventdata, handles)
% hObject    handle to clearTriggerParametersSelection2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clearParametersSelection(hObject, handles);


% --- Executes on button press in testing.
function testing_Callback(hObject, eventdata, handles)
% hObject    handle to testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of testing

testing = get(hObject, 'Value');

if testing
    set(handles.fitAllCurve, 'Visible', 'On');
else
    set(handles.fitAllCurve, 'Visible', 'Off');
end
handles = updateGUI(handles);

% Update handles
guidata(hObject, handles);


% --- Executes on button press in calcControlParams.
function calcControlParams_Callback(hObject, eventdata, handles)
% hObject    handle to calcControlParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = updateGUI(handles);

% Update handles
guidata(hObject, handles);


function controlStrength2_Callback(hObject, eventdata, handles)
% hObject    handle to controlStrength2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of controlStrength2 as text
%        str2double(get(hObject,'String')) returns contents of controlStrength2 as a double
handles = checkNewInput(hObject, handles, 'manStrength2VarData');
handles = updateUnitCategorization(handles); % Update which unit is unit A and unit B

handles = setControlVoltages(handles);
set(handles.controlV2, 'String', num2str(handles.curControlVoltages(2)));

% Redraw the SD curves
handles = drawSDPlot(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function controlStrength2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to controlStrength2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function controlDuration2_Callback(hObject, eventdata, handles)
% hObject    handle to controlDuration2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of controlDuration2 as text
%        str2double(get(hObject,'String')) returns contents of controlDuration2 as a double
handles = checkNewInput(hObject, handles, 'manDuration2VarData');
handles = updateUnitCategorization(handles); % Update which unit is unit A and unit B

% Redraw the SD curves
handles = drawSDPlot(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function controlDuration2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to controlDuration2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function controlStrength1_Callback(hObject, eventdata, handles)
% hObject    handle to controlStrength1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of controlStrength1 as text
%        str2double(get(hObject,'String')) returns contents of controlStrength1 as a double
handles = checkNewInput(hObject, handles, 'manStrength1VarData');
handles = updateUnitCategorization(handles); % Update which unit is unit A and unit B

handles = setControlVoltages(handles); % Calculate the voltage needed to generate this power

% Redraw the SD curves
handles = drawSDPlot(handles);

guidata(hObject, handles);

function handles = updateUnitCategorization(handles)
% Update the categorization of each unit, showing which is unit A and which
% is unit B.  By convention, unit A prefers short/strong pulses, while unit
% B prefers long/weak pulses
g1 = str2double(get(handles.controlStrength1, 'String'));
t1 = str2double(get(handles.controlDuration1, 'String'));
g2 = str2double(get(handles.controlStrength2, 'String'));
t2 = str2double(get(handles.controlDuration2, 'String'));

if g1 > g2 && t1 < t2
    % g1 is higher and t1 is shorter
    handles.unit1IsUnitA = true;
else
    % g1 is lower and t1 is longer
    handles.unit1IsUnitA = false;
end

% --- Executes during object creation, after setting all properties.
function controlStrength1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to controlStrength1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function controlDuration1_Callback(hObject, eventdata, handles)
% hObject    handle to controlDuration1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of controlDuration1 as text
%        str2double(get(hObject,'String')) returns contents of controlDuration1 as a double
handles = checkNewInput(hObject, handles, 'manDuration1VarData');
handles = updateUnitCategorization(handles); % Update which unit is unit A and unit B

% Redraw the SD curves
handles = drawSDPlot(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function controlDuration1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to controlDuration1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SDSelect2.
function SDSelect2_Callback(hObject, eventdata, handles)
% hObject    handle to SDSelect2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select the new stimulations
handles = SDSelectCallback(hObject, handles);

% Update handles
guidata(hObject, handles);

% --- Executes on button press in SDSelect1.
function SDSelect1_Callback(hObject, eventdata, handles)
% hObject    handle to SDSelect1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select the new stimulations
handles = SDSelectCallback(hObject, handles);

% Update handles
guidata(hObject, handles);

function handles = SDSelectCallback(hObject, handles)
% CHECKED %

% Allow the user to define a rectangle in the SD plot inside which all
% stimulations will be selected

neuronNum = get(hObject, 'UserData'); % The neuron currently being worked on
axisHandle = handles.(['SDCurve' num2str(neuronNum)]); % Handle to the SD curve plot

% Get the strength and duration ranges from the user
[x1, y1, button] = ginput(1); % Get the first x and y point
if button ~= 1
    % If the user didn't click the mouse, and instead hit a key, assume
    % that they want to cancel the selection
    return;
end

set(gcf, 'WindowButtonMotionFcn', @(object, eventdata) drawSelectionRect(handles, axisHandle, x1, y1));
% pause();
% [x2, y2] = ginput(1); % Get the second x and y point
set(gcf, 'WindowButtonDownFcn', @(object, eventData) getCurLoc(object, handles, axisHandle));

% Wait for the second [x y] pair to be returned
if waitforbuttonpress == 1
    % If the user didn't click the mouse, and instead hit a key, assume
    % that they want to cancel the selection
    set(gcf, 'WindowButtonMotionFcn', '');
    set(gcf, 'WindowButtonDownFcn', '');
    handles = drawSDPlot(handles);
    return;
end
handles = guidata(hObject); % Constantly refresh handles, to see if getCurLoc is done yet)
x2 = handles.lastMousePos(1);
y2 = handles.lastMousePos(2);
% disp([x1 x2 y1 y2]);

% Return the button functions to normal (do nothing)
set(gcf, 'WindowButtonMotionFcn', '');
set(gcf, 'WindowButtonDownFcn', '');

% Sort the x and y values
xl = min([x1 x2]); % x lower
xu = max([x1 x2]); % x upper
yl = min([y1 y2]); % y lower
yu = max([y1 y2]); % y upper

% Use the strength and duration ranges to select neurons
strengthData = handles.allNeuronsStruct{neuronNum}.resultsData(:,1);
durationData = handles.allNeuronsStruct{neuronNum}.resultsData(:,2);
selectedStims = find((yl < strengthData) & (strengthData < yu) & (xl < durationData) & (durationData < xu));

% Set these stimulation values to the listbox
listBoxHandle = handles.(['neuron' num2str(neuronNum) 'StimsListbox']);
set(listBoxHandle, 'Value', selectedStims); % Set the values to the listbox (selecting the stims that were found here)
handles = stimsListboxCallback(handles); % Call the callback function for this object, as if the user just entered new values

function getCurLoc(hObject, handles, axisHandle)
% CHECKED %
% This function sets the current cursor position to handles.lastMousePos
C = get(axisHandle, 'CurrentPoint'); % Get the mouse position from the axis
mousePos = [C(1,1) C(1,2)]; % Extract the cursor position
handles.lastMousePos = mousePos; % Store the mouse position in handles

% Update handles
guidata(hObject, handles);


function drawSelectionRect(handles, axisHandle, x1, y1)
% CHECKED %

% This is the function that will be called each time the mouse moves while
% the user is selecting a rectangle on the SD plot

% Draw all underlying data on the plot
handles = drawSDPlot(handles);

C = get(axisHandle, 'CurrentPoint'); % Get the current mouse position
xCur = C(1,1);
yCur = C(1,2);

hold(axisHandle, 'on'); % Turn hold on, so the new rectangle can be drawn
axes(axisHandle);
drawRect(axisHandle, [x1 xCur], [y1 yCur], 'k--');
hold(axisHandle, 'off'); % Turn hold off

function drawRect(axisHandle, x, y, ls)
% CHECKED %
% This function will draw a rectangle on the given axis.  The rectangle is
% defined by two opposite corners, (x1, y1) and (x2, y2).  The input is
% defined as: x = [x1 x2], y = [y1 y2]. ls is a linespec string that will
% be passed to the plot function.

x1 = x(1);
x2 = x(2);
y1 = y(1);
y2 = y(2);

plot(axisHandle, [x1 x1], [y1 y2], ls, [x1 x2], [y1 y1], ls, [x1 x2], [y2 y2], ls, [x2 x2], [y1 y2], ls);

% function selectFromAxis_Callback(hObject, eventdata, handles)
% handles = selectFromAxisCallback(hObject, eventdata, handles, handles.manStimStrength, handles.manStimDuration);
% guidata(hObject, handles);
% 
% % --- Executes on button press in selectFromAxis1.
% function selectFromAxis1_Callback(hObject, eventdata, handles)
% % hObject    handle to selectFromAxis1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% handles = selectFromAxisCallback(hObject, eventdata, handles, handles.controlStrength1, handles.controlDuration1);
% guidata(hObject, handles);
% 
% % --- Executes on button press in selectFromAxis2.
% function selectFromAxis2_Callback(hObject, eventdata, handles)
% % hObject    handle to selectFromAxis2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% handles = selectFromAxisCallback(hObject, eventdata, handles, handles.controlStrength2, handles.controlDuration2);
% guidata(hObject, handles);
% 
% function handles = selectFromAxisCallback(hObject, eventdata, handles, strHandle, durHandle)
% % CHECKED %
% 
% % Get strength/duration pair from one of the SD plots
% [newDur, newStr, button] = ginput(1);
% 
% if button ~= 1
%     % If the user didn't click the mouse, and instead hit a key, assume
%     % that they want to cancel the selection
%     return;
% end
% 
% % Write this strength/duration pair into the manual stimulation edit boxes
% set(strHandle, 'String', num2str(newStr));
% set(durHandle, 'String', num2str(newDur));
% 
% % Call each callback function so that it will check the new strength and
% % duration
% manStimStrength_Callback(strHandle, eventdata, handles);
% handles = guidata(hObject); % Get the new handles structure that the previous callback just saved
% manStimDuration_Callback(durHandle, eventdata, handles);
% 
% % Redraw the SD curves
% handles = drawSDPlot(handles);


% --- Executes on key press with focus on neuron1StimsListbox and none of its controls.
function neuron1StimsListbox_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to neuron1StimsListbox (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% handles = listBoxKeyPressCallback(hObject, eventdata, handles);

% Update handles
guidata(hObject, handles);

% --- Executes on key press with focus on neuron2StimsListbox and none of its controls.
function neuron2StimsListbox_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to neuron2StimsListbox (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% handles = listBoxKeyPressCallback(hObject, eventdata, handles);

% Update handles
guidata(hObject, handles);

function playWarn(handles)
if handles.doPlayWarn
    try
        play(handles.warnAudio);
    catch
        disp('Could not play audio');
    end
end

function handles = listBoxKeyPressCallback(hObject, eventdata, handles)
% This callback will allow the user to delete some stimulations using the
% backspace key.  It will confirm with the user before deleting.

neuronNum = get(hObject, 'UserData');

% if strcmp(eventdata.Key, 'backspace')
%     % User hit the backspace key. Make sure that they want to delete this stimulations
%     numStims = length(handles.selectedStims{neuronNum}); % The number of stimulations selected
%     playWarn(handles);
%     choice = questdlg(sprintf('Are you sure that you want to delete %d stimulations?  This cannot be undone. ("Yes" to delete stimulations, "No" to not delete)', numStims), 'Warning: Delete Stimulations?', 'Yes', 'No', 'No');
%     if strcmp(choice, 'Yes')
%         % User confirmed wanting to delete the stimulations.  Remove them from the resultsData for this neuron
%         handles = deleteStims(neuronNum, handles.selectedStims{neuronNum}, handles);
%     end
% end

drawVoltageTraces(handles, neuronNum);
handles = stimsListboxCallback(handles);


% --- Executes on button press in alignSpikes1.
function alignSpikes1_Callback(hObject, eventdata, handles)
% hObject    handle to alignSpikes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alignSpikes1
neuronNum = get(hObject, 'UserData');
handles = drawVoltageTraces(handles, neuronNum);
guidata(hObject, handles);

% --- Executes on button press in alignSpikes2.
function alignSpikes2_Callback(hObject, eventdata, handles)
% hObject    handle to alignSpikes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alignSpikes2
neuronNum = get(hObject, 'UserData');
handles = drawVoltageTraces(handles, neuronNum);
guidata(hObject, handles);


% --- Executes on button press in unfreeze.
function unfreeze_Callback(hObject, eventdata, handles)
% hObject    handle to unfreeze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Re-enable the GUI if frozen (if there was an exception during one of the
% main routines, and the entire GUI is grayed out)
handles = updateGUI(handles);

% Update handles
guidata(hObject, handles);


% --- Executes on button press in viewSequenceTimeseries.
function viewSequenceTimeseries_Callback(hObject, eventdata, handles)
% CHECKED %
% hObject    handle to viewSequenceTimeseries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Check if the figure is currently open
if ishandle(handles.sequenceTimeseriesFig)
    % If so, just bring it to the front
    figure(handles.sequenceTimeseriesFig);
else
    % If not, create a new figure and assign it to the sequenceTimeseriesFig
    % variable
    handles.sequenceTimeseriesFig = figure;
end

% Now that we have a figure, redraw it
redrawSequenceTimeseries(handles);

guidata(hObject, handles);

function redrawSequenceTimeseries(handles)
if isempty(handles.sequenceTimeseriesFig) || ~ishandle(handles.sequenceTimeseriesFig)
    % If there is no handle, don't draw it
    return;
end

% If the figure already exists, draw the information on the sequence
% results figure
plotFontSize = 18;
linewidth = 2;
legendFontSize = 18;
markerSize = 15;
% doAdaptive = get(handles.useAdaptive, 'Value');

% Draw a new axis
clf(handles.sequenceTimeseriesFig); % Clear the figure
ax = axes('Parent', handles.sequenceTimeseriesFig);
hold(ax, 'on');

% Calculate the needed information
[TP, FA, stimP, timeAxis, numSequences, errStr] = calculateSequenceProbabilities(handles);

% Check if there was an error
if ~isempty(errStr)
    % Do not plot the figure, but display the error message
    title(ax, sprintf('Error: %s', errStr), 'fontsize', plotFontSize);
    return;
end

% Plot the probabilities
% sampInds = 1:10:length(timeAxis);
legendEntries = {'Neuron 1 True Positive', 'Neuron 2 True Positive', 'Neuron 1 False Positive', 'Neuron 2 False Positive'};
% h = plot(ax, timeAxis, TP(1,:), 'b-', timeAxis, TP(2,:), 'g-', timeAxis, FA(1,:), 'b--', timeAxis, FA(2,:), 'g--', timeAxis(sampInds), TP(1,sampInds), 'bx', timeAxis(sampInds), TP(2,sampInds), 'gx', timeAxis(sampInds), FA(1,sampInds), 'bo', timeAxis(sampInds), FA(2,sampInds), 'go', 'linewidth', linewidth, 'markersize', markerSize);
h = plot(ax, timeAxis, TP(1,:), 'bx-', timeAxis, TP(2,:), 'ro-', timeAxis, FA(1,:), 'bo--', timeAxis, FA(2,:), 'rx--', 'linewidth', linewidth, 'markersize', markerSize);
% h = h(5:8); % Only keep references to some of the plots (for the legend)
% if ~any(isnan(adaptationBlockChange(:))) && logical(length(timeAxis) >= max(adaptationBlockChange(:)))
%     numAdaptations = length(adaptationBlockChange);
%     
%     % Add lines to show when adaptation occurs
%     j = plot(ax, repmat(numAdaptations', 2, 1), repmat([0;1], 1, numAdaptations), 'm:', 'linewidth', 2);
%     h = [h;j];
%     legendEntries = [legendEntries, 'Adaptation Block Borders'];
%     % Create a patch that will highlight when adaptation is on
%     numAdaptationBlocks = size(adaptationBlockChange, 1);
%     
%     % If there is at least one adaptation block to indicate
%     if numAdaptationBlocks > 0
%         allX = [timeAxis(adaptationBlockChange(:, 1))' timeAxis(adaptationBlockChange(:, 2))']'; % Arrange the times of onset and offset of adaptation
%         shapedX = zeros(2, numAdaptationBlocks);
%         shapedX(:) = allX;
%         X = repelem(shapedX, 2, 1);
%         Y = repmat([1 0 0 1]', 1, numAdaptationBlocks);
%         
%         hold(ax, 'on');
%         %     a = plot(ax, repmat(adaptationBlockChange,2,1), repmat([0;1], 1, size(GTChangeLocs,2)), 'm--', 'linewidth', linewidth);
%         patch(ax, X, Y, [0 0 1], 'facealpha', .1, 'edgealpha', 0);
%         hold(ax, 'off');
%         %    h(end + 1) = a;
%         %    legendEntries{end + 1} = 'New GT Pairs';
%     end
% end
% if ~isnan(adaptationEvents)
%     x = repelem(timeAxis(adaptationEvents), 2, 1);
%     y = repmat([0;1], 1, length(adaptationEvents));
%     
%     hold(ax, 'on');
%     a = plot(ax, x, y, 'm:');
%     hold(ax, 'off');
%     
%     h(end + 1) = a(1);
%     legendEntries{end + 1} = 'Adaptation Update Event';
% end

title(ax, sprintf('Success over time (n = %d seqs, P[Sa] = %0.3f, P[Sb] = %0.3f)', numSequences, stimP(1), stimP(2)), 'fontsize', plotFontSize);
xlabel(ax, 'Sequence Start Time (minutes)', 'fontsize', plotFontSize);
ylabel(ax, 'Probability', 'fontsize', plotFontSize)
set(ax, 'fontsize', plotFontSize);
legend(ax, h, legendEntries, 'fontsize', legendFontSize, 'location', 'best');

% function redrawSequencePSpace(handles)
% % CHECKED %
% % If the figure already exists, draw the information on the sequence
% % results figure
% plotFontSize = 18;
% linewidth = 2;
% legendFontSize = 18;
% markerSize = 15;
% 
% if (isempty(handles.sequencePSpaceFig1)|| isempty(handles.sequencePSpaceFig2)) || (~ishandle(handles.sequencePSpaceFig1) || ~ishandle(handles.sequencePSpaceFig2))
%     % If there is no handle, don't draw it
%     return;
% end
% 
% % Draw a new axis
% clf(handles.sequencePSpaceFig1); % Clear the figure
% clf(handles.sequencePSpaceFig2); % Clear the figure
% a1 = axes('Parent', handles.sequencePSpaceFig1);
% axes(a1);
% colorbar;
% a2 = axes('Parent', handles.sequencePSpaceFig2);
% axes(a2);
% colorbar;
% 
% % Calculate the needed information
% [TP, FA, stimP, timeAxis, numSequences, errStr] = calculateSequenceProbabilities(handles);
% 
% % Check if there was an error
% if ~isempty(errStr)
%     % Do not plot the figure, but display the error message
%     title(a1, sprintf('Error: %s', errStr), 'fontsize', plotFontSize);
%     title(a2, sprintf('Error: %s', errStr), 'fontsize', plotFontSize);
%     return;
% end
% 
% % Plot the probabilities
% tEnd = timeAxis(end);
% 
% % Plot the trajectory through probablity space
% z = zeros(size(timeAxis));
% 
% axes(a1)
% surface([TP(1,:);TP(1,:)],[FA(2,:);FA(2,:)],[z;z],[timeAxis;timeAxis],'facecol','no','edgecol','interp','linew',2);
% title(sprintf('Probability Trajectory (window = %d, n = %d, P[Sa] = %0.3f, P[Sb] = %0.3f)', smoothingSize, numSequences, stimP(1), stimP(2)), 'fontsize', plotFontSize);
% xlabel('True Positives', 'fontsize', plotFontSize);
% ylabel('False Alarms', 'fontsize', plotFontSize)
% set(a1, 'XColor', handles.blue, 'YColor', handles.green, 'fontsize', plotFontSize);
% ylim([0 1]);
% xlim([0 1]);
% caxis([0 tEnd])
% 
% axes(a2);
% surface([TP(2,:);TP(2,:)],[FA(1,:);FA(1,:)],[z;z],[timeAxis;timeAxis],'facecol','no','edgecol','interp','linew',2);
% title(sprintf('Probability Trajectory (window = %d, n = %d, P[Sa] = %0.3f, P[Sb] = %0.3f)', smoothingSize, numSequences, stimP(1), stimP(2)), 'fontsize', plotFontSize);
% xlabel('True Positives', 'fontsize', plotFontSize);
% ylabel('A False Alarms', 'fontsize', plotFontSize)
% set(a2, 'XColor', handles.green, 'YColor', handles.blue, 'fontsize', plotFontSize);
% ylim([0 1]);
% xlim([0 1]);
% caxis([0 tEnd])

function [TP, FA, stimP, timeAxis, numSequences, errStr] = calculateSequenceProbabilities(handles)
% CHECKED %

% Preallocate values for each output, in case of error
TP = nan;
FA = nan;
stimP = nan;
timeAxis = nan;
numSequences = nan;
adaptationBlockChange = nan;
adaptationEvents = nan;
errStr = '';
% smoothingSize = str2double(get(handles.smoothingSize, 'String')); % Number of sequences to smooth across (will calculate probabilities for this number of sequences in a boxcar smoothing manner)
% trueSmoothingRadius = floor(smoothingSize/2); % Calculate size of smoothing window
% trueSmoothingDiameter = 2*trueSmoothingRadius;

struct1 = handles.allNeuronsStruct{1};
struct2 = handles.allNeuronsStruct{2};

% Check that sequence results can be drawn for this pair of neurons

% Check that sequenceInfo exists for both neurons
if isempty(struct1.sequenceInfo) || isempty(struct2.sequenceInfo)
    errStr = 'At least one unit has had no sequences yet.';
    return;
end

% Extract sequence data (to compare and make sure that these two
% neuronStructs are compatible)
numStimsPerSeq = handles.numSlots;
% Find where the last trialData left off

% sIOffset1 = fix((struct1.sequenceInfoCurLoc - 2)/numSeqsPerExp)*numSeqsPerExp + 1;
% sIOffset2 = fix((struct2.sequenceInfoCurLoc - 2)/numSeqsPerExp)*numSeqsPerExp + 1;
if isempty(struct1.trialData) || length(struct1.trialData) == 1
    sIOffset1 = 1;
    sIOffset2 = 1;
else
    % If there have been multiple trials started
    lastStimOfLastTrial1 = max(struct1.trialData(end - 1).stimInds(:));
    lastStimOfLastTrial2 = max(struct2.trialData(end - 1).stimInds(:));
    sIOffset1 = lastStimOfLastTrial1/numStimsPerSeq + 1;
    sIOffset2 = lastStimOfLastTrial2/numStimsPerSeq + 1;
end

% sIOffset1 = fix((struct1.sequenceInfoCurLoc - 2)/numSeqsPerExp)*numSeqsPerExp + 1;
% sIOffset2 = fix((struct2.sequenceInfoCurLoc - 2)/numSeqsPerExp)*numSeqsPerExp + 1;
seqInf1 = struct1.sequenceInfo(sIOffset1:(struct1.sequenceInfoCurLoc - 1), :);
seqInf2 = struct2.sequenceInfo(sIOffset2:(struct2.sequenceInfoCurLoc - 1), :);
aTimes = seqInf1(:, 1);
bTimes = seqInf2(:, 1);

% Check that the same number of sequences occurred with both neurons
if length(aTimes) ~= length(bTimes) || all(bTimes ~= aTimes)
    errStr = 'Sequence information is not the same between neurons.';
    return;
end

% Extract all needed data
% Extract the correct form of aExpected/bExpected
sequenceLength = handles.numSlots; % The number of stimulations in each sequence
eR1T = seqInf1(:,2:(1 + sequenceLength))'; % The expected results that the neuron will have for each sequence
eR2T = seqInf2(:,2:(1 + sequenceLength))'; % --This should be the same as ~aExpected--
doAdapting = seqInf1(:, 4 + sequenceLength); % Which type of block (adapting or static) is this sequence in?
aControlGTs = seqInf1(:, (1:2) + (sequenceLength + 1)); % The control GT used to activate this neuron
bControlGTs = seqInf2(:, (1:2) + (sequenceLength + 1)); % The control GT used to activate this neuron
% aResultsTimes = struct1.resultsData(:,4); % When was this stimulation done?
% bResultsTimes = struct2.resultsData(:,4);

sequenceTimes = aTimes; % The times that a sequence occurs
numSequences = length(aTimes); % It is assumed that both a and b have the same number of sequences, and that they happen at identical times
% newFileTime = max(aResultsTimes(1), bResultsTimes(1));
numAllStims = numSequences*sequenceLength;

eR1Col = eR1T(:);
eR2Col = eR2T(:);

sizeAdaptationBlock = handles.adaptiveBufferReplacementSize; % Number of sequences per adaptation block
allBlockBorders = 1:sizeAdaptationBlock*sequenceLength:(numAllStims + 1);
numAdaptationBlocks = length(allBlockBorders) - 1;
allBlocks = 1:numAdaptationBlocks;
TP1_block = ones(numAdaptationBlocks, 1);
TP2_block = ones(numAdaptationBlocks, 1);
FA1_block = ones(numAdaptationBlocks, 1);
FA2_block = ones(numAdaptationBlocks, 1);
% timeAxis = ones(numAdaptationBlocks, 1);

%% Preprocess sequence information
% Get indices of first of sequence
% [~, aFirstInds] = ismember(aTimes, aResultsTimes); % The indicies in resultsData that correspond to the first stimulation of each sequence
% [~, bFirstInds] = ismember(bTimes, bResultsTimes); % The indicies in resultsData that correspond to the first stimulation of each sequence

% Make an index matrix for each neuron, where each value represents the
% index in the resultsData matrix (e.g. 5 represents the 5th
% stimulation in resultsData). Rows are sequences, columns are indices
% in that sequence.
% aSequenceInds = bsxfun(@plus, aFirstInds, 0:(sequenceLength - 1));
% bSequenceInds = bsxfun(@plus, bFirstInds, 0:(sequenceLength - 1));

% Make a stimulation results matrix for each neuron
% trialStimIndsT1 = struct1.trialData(end).stimInds';
% trialStimIndsT2 = struct2.trialData(end).stimInds';
% trialStimInds1 = trialStimIndsT1(:);
% trialStimInds2 = trialStimIndsT2(:);
stimInds1 = find(struct1.resultsData(:, 4) == sequenceTimes(1), 1):(find(struct1.resultsData(:, 4) == sequenceTimes(end), 1, 'last') + sequenceLength - 1);
stimInds2 = find(struct2.resultsData(:, 4) == sequenceTimes(1), 1):(find(struct2.resultsData(:, 4) == sequenceTimes(end), 1, 'last') + sequenceLength - 1);
sR1Col = struct1.resultsData(stimInds1, 3); % Did a spike occur on this stimulation?
sR2Col = struct2.resultsData(stimInds2, 3);

% Calculate number of timesteps after smoothing
% thisTrueSmoothingDiameter = min(trueSmoothingDiameter, numSequences);
% trueSmoothingSize = thisTrueSmoothingDiameter + 1;
% thisTrueSmoothingRadius = ceil(thisTrueSmoothingDiameter/2);
% numTimeSteps = numSequences - smoothingSize + 1;
% if numTimeSteps < 1
%     numTimeSteps = 1;
%     smoothingSize = 1;
% end

% Prepare for calculation of when adaptation events occur
% allGTs = [aControlGTs bControlGTs];

% Calculate a time axis
% sequenceTimes = (sequenceTimes' - newFileTime)/60; % Adjust sequence times to minutes after the most recent newly created file was made
% timeAxis = sequenceTimes((floor(smoothingSize/2) + 1):(floor(smoothingSize/2) + numTimeSteps));
% timeAxis = sequenceTimes((end - numTimeSteps + 1):end);

% TP = zeros(2,numTimeSteps); % P{a|a}
% FA = zeros(2,numTimeSteps); % P{a|b}
% smoothedDoAdapting = zeros(numTimeSteps, 1);
% smoothedDoAdapting = doAdapting;

% Calculate smoothing window
% if smoothingSize > 1
%     n = 0:(smoothingSize - 1);
%     w = cos((pi*n)/(smoothingSize - 1) - pi/2);
%     %         cosW = repmat(cosW', numSlots, 1); % Smoothing window is applied to each stim individually, but it is applied to the columnization of TPLogicalX./eRX, each of which are numSequences x numSlots, which means the column contains [seq1slot1, seq2slot1, ..., seqNslot1, seq1slot2, ...] (lists all sequences' first slot, then all sequences' second slot, etc.), so the smoothing window is repeated, because it must be applied according to each sequence, not the slot number
%     w = w'/sum(w); % Normalize the window so it does not scale the data
% else
%     w = 1;
% end

%% Calculate the smoothed probabilities over time
for blockInd = allBlocks
    %% Extract needed information
    % Sequences that will be used in this calculation
%     sequenceInds = blockNum:(blockNum + smoothingSize - 1);
    thisAllInds = allBlockBorders(blockInd):(allBlockBorders(blockInd + 1) - 1);
    thisERLogicals1 = logical(eR1Col(thisAllInds));
    thisERLogicals2 = logical(eR2Col(thisAllInds));
    thisERInds1 = thisAllInds(thisERLogicals1);
    thisERInds2 = thisAllInds(thisERLogicals2);
    
%     thisInds1 = thisAllInds(thisERInds1);
%     thisInds2 = thisAllInds(thisERInds2);
    
    % Neuron 1 stimulations
    thisSR11 = sR1Col(thisERInds1);
    thisSR21 = sR2Col(thisERInds1);
    TP1_block(blockInd) = mean(thisSR11);
    FA2_block(blockInd) = mean(thisSR21);
    
    % Neuron 2 stimulations
    thisSR12 = sR1Col(thisERInds2);
    thisSR22 = sR2Col(thisERInds2);
    TP2_block(blockInd) = mean(thisSR22);
    FA1_block(blockInd) = mean(thisSR12);
    
%     timeAxis(blockInd) = sequenceTimes((blockInd - 1)*sizeAdaptationBlock + 1);
    
%     % Determine if this segment is considered adapting or not
%     smoothedDoAdapting(blockNum) = logical(round(mean(doAdapting(sequenceInds))));
%     
%     %% Analyze sequence information (determine results for stimulation in sequence)
%     % Find true positives
%     aTPLogical = thisAStimResults & thisAExpected; % A spike when it was expected
%     bTPLogical = thisBStimResults & thisBExpected;
%     
%     % Find false positives
%     aFPLogical = thisAStimResults & ~thisAExpected; % A spike when it was not expected
%     bFPLogical = thisBStimResults & ~thisBExpected;
%     
%     %% Calculate results (get probabilities of each condition)
%     % Total true positives probabilities
%     TP(1, blockNum) = sum(w.*(sum(aTPLogical, 2)./sum(thisAExpected, 2))); % P{a|a}
%     TP(2, blockNum) = sum(w.*(sum(bTPLogical, 2)./sum(thisBExpected, 2))); % P{b|b}
%     
%     % Total false positives probabilities
%     FA(1, blockNum) = sum(w.*(sum(aFPLogical, 2)./sum(~thisAExpected, 2))); % P{a|b}
%     FA(2, blockNum) = sum(w.*(sum(bFPLogical, 2)./sum(~thisBExpected, 2))); % P{b|a}
end

timeAxis = (sequenceTimes((allBlocks - 1)*sizeAdaptationBlock + 1) - sequenceTimes(1))/60;
FA = [FA1_block'; FA2_block'];
TP = [TP1_block'; TP2_block'];



% smoothedDoAdapting = [0; smoothedDoAdapting;0]; % Zero pad the doAdapting vector, which allows us to calculate where to draw a patch to highlight in the right place
% doAdapting = [0;doAdapting;0];

%% Calculate the places where the adaptation block changed
% adaptationBlockChange = [(find(diff(doAdapting) > 0)) (find(diff(doAdapting) < 0)) - 1]; % The first column indicates sequences where the block type moved from static -> adaptation.  The second column indicates sequences where the block type moved from adaptation -> static.
adaptationBlockChange = allBlockBorders(1:(end - 1));
% adaptationEvents = 1 + find(any(diff(allGTs, 1) ~= 0, 2));

% Calculate probabilities of stimulating neuron 1 or 2
% Define helper functions
fractionOfTotal = @(mat) sum(logical(mat(:)))/numel(mat); % This function gets the percentage of logical true's in a 2D matrix
stimP = [fractionOfTotal(eR1Col) fractionOfTotal(eR2Col)];

function sequenceReps_Callback(hObject, eventdata, handles)
% hObject    handle to sequenceReps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sequenceReps as text
%        str2double(get(hObject,'String')) returns contents of sequenceReps as a double
handles = checkNewInput(hObject, handles, 'numSequencesVarData', true);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sequenceReps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sequenceReps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stopSequence.
function stopSequence_Callback(hObject, eventdata, handles)
% hObject    handle to stopSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% The user is attempting to stop the sequence from running
handles.stoppingSequence = true;
guidata(hObject, handles);

% Make this button invisible again
set(hObject, 'Visible', 'off');

% Change the GUI so that it shows the user that the sequences are stopping
% set(handles.controlSequenceStart, 'String', 'Stopping...');
% set(handles.controlSequenceStart, 'BackgroundColor', [1 1 0]);


% --- Executes on slider movement.
function stimStatusSlider_Callback(hObject, eventdata, handles)
% hObject    handle to stimStatusSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function stimStatusSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimStatusSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in manStim2.
function manStim2_Callback(hObject, eventdata, handles)
% hObject    handle to manStim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get strength and duration from the GUI
strength = str2double(get(handles.controlStrength2, 'String'));
duration = str2double(get(handles.controlDuration2, 'String'));

manStim(strength, duration, hObject, handles);

% --- Executes on button press in manStim1.
function manStim1_Callback(hObject, eventdata, handles)
% hObject    handle to manStim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get strength and duration from the GUI
strength = str2double(get(handles.controlStrength1, 'String'));
duration = str2double(get(handles.controlDuration1, 'String'));

manStim(strength, duration, hObject, handles);



function smoothingSize_Callback(hObject, eventdata, handles)
% hObject    handle to smoothingSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smoothingSize as text
%        str2double(get(hObject,'String')) returns contents of smoothingSize as a double
handles = checkNewInput(hObject, handles, 'smoothingSizeVarData', true);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function smoothingSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothingSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in playWarnCheck.
function playWarnCheck_Callback(hObject, eventdata, handles)
% hObject    handle to playWarnCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of playWarnCheck

handles.doPlayWarn = get(hObject, 'Value'); % Set whether or not the program should play warnings
guidata(hObject, handles); % Update handles


function warnVolumeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to warnVolumeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of warnVolumeEdit as text
%        str2double(get(hObject,'String')) returns contents of warnVolumeEdit as a double
handles = checkNewInput(hObject, handles, 'warnVolumeVarData');

% Check to see if the value has changed
newVal = str2double(get(hObject, 'String'))/100;
if handles.warnVolume ~= newVal
    % Divide the input volume by 100 (so it's between 0 and 1), and save it
    handles.warnVolume = newVal;
    
    % Apply this new volume to the audio object and save it
    handles.warnAudio = audioplayer(handles.warnVolume*handles.rawWarn, 22050);
end
guidata(hObject, handles);
playWarn(handles); % Test the new volume

% --- Executes during object creation, after setting all properties.
function warnVolumeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to warnVolumeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function alignSpikes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alignSpikes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in usePCACheck1.
function usePCACheck1_Callback(hObject, eventdata, handles)
% hObject    handle to usePCACheck1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usePCACheck1

% Reevaluate the result status of all the stimulations
neuronNum = get(hObject, 'UserData');
numStims = handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc - 1;
handles = recheckResults(handles, neuronNum, 1:numStims); % Recheck ALL stimulations
checkDefaults(handles);
guidata(hObject, handles);

% --- Executes on button press in viewRaw2.
function viewRaw2_Callback(hObject, eventdata, handles)
% hObject    handle to viewRaw2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of viewRaw2

% Redraw the voltage traces
neuronNum = get(hObject, 'UserData');
handles = drawVoltageTraces(handles, neuronNum);
guidata(hObject, handles);

% --- Executes on button press in viewRaw1.
function viewRaw1_Callback(hObject, eventdata, handles)
% hObject    handle to viewRaw1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of viewRaw1

% Redraw the voltage traces
neuronNum = get(hObject, 'UserData');
handles = drawVoltageTraces(handles, neuronNum);
guidata(hObject, handles);

% --- Executes on button press in usePCACheck2.
function usePCACheck2_Callback(hObject, eventdata, handles)
% hObject    handle to usePCACheck2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usePCACheck2

% Reevaluate the result status of all the stimulations
neuronNum = get(hObject, 'UserData');
numStims = handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc - 1;
handles = recheckResults(handles, neuronNum, 1:numStims); % Recheck ALL stimulations
checkDefaults(handles);
guidata(hObject, handles);

% --- Executes on button press in viewShape1.
function viewShape1_Callback(hObject, eventdata, handles)
% hObject    handle to viewShape1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function h = createFig(figNum)
% CHECKED %
% Creates a figure, identified by the specified figure number, and returns
% an axis.  Ensures only one figure with the specified ID exists.
h = figure(figNum);
clf(h);

% --- Executes on selection change in codeSelect1.
function codeSelect1_Callback(hObject, eventdata, handles)
% hObject    handle to codeSelect1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function codeSelect1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to codeSelect1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in viewShape2.
function viewShape2_Callback(hObject, eventdata, handles)
% hObject    handle to viewShape2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% neuronNum = get(hObject, 'UserData');
% createShapeFig(handles, neuronNum);

% --- Executes on selection change in codeSelect2.
function codeSelect2_Callback(hObject, eventdata, handles)
% hObject    handle to codeSelect2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function codeSelect2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to codeSelect2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autoSave2.
function autoSave2_Callback(hObject, eventdata, handles)
% hObject    handle to autoSave2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoSave2


% --- Executes during object creation, after setting all properties.
function saveNeuron_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveNeuron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in autosave.
function autoSave_Callback(hObject, eventdata, handles)
% hObject    handle to autosave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autosave


% --- Executes on button press in setToRecording.
function setToRecording_Callback(hObject, eventdata, handles)
% hObject    handle to setToRecording (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set the TDT system to recording mode
if handles.TDTConnected
    % If TDT is connected (it should be connected before this button is
    % available to press, but redundancy isn't a bad thing...)
    if handles.DA.GetSysMode ~= 3
        % If the system mode is NOT recording
        handles.DA.SetSysMode(3); % Set the mode to recording
    end
    
    handles.TDTRecording = true;
end

handles = updateGUI(handles); % Update the GUI

guidata(hObject, handles); % Update handles

function handles = updateTrialMatrices(handles)
% This function creates a new expectedResults matrix,
% which will eventually be sent to TDT

% Gather the meta-data
numStimsPerSequence = str2double(get(handles.stimsPerSequence, 'String'));
numRuns = str2double(get(handles.numRuns, 'String'));

% Calculate the expected results of the stimulations and add it to handles
[runOrders, expectedResultsAll] = createSequenceOrders(numStimsPerSequence, numRuns);
handles.curExpectedResults = expectedResultsAll;
handles.curRunOrders = runOrders;
% handles.curInterStimPeriods = interStimPeriodsAll;

% Calculate how long this trial would take
handles = updateTrialTimes(handles);

function handles = updateTrialTimes(handles)
% CHECKED %
% This function recalculates the amount of time that the current expected
% results and inter-stimulation periods matrices would take to complete,
% and then updates the GUI to reflect this calculation

% Get the inter-stim periods and expected results
expectedResults = handles.curExpectedResults;

if isempty(expectedResults)
    % If expectedResults is empty, exit this function
    return;
end

% Gather the other important information
numRuns = str2double(get(handles.numRuns, 'String'));
interSeqPeriodLow = str2double(get(handles.seqPeriodLow, 'String'));
interSeqPeriodHigh = str2double(get(handles.seqPeriodHigh, 'String'));

% Calculate the times, from both the inter-stimulation periods and the
% inter-sequence periods
baseTime = (str2double(get(handles.stimPeriods, 'String'))/1000)*(size(expectedResults, 1)*(size(expectedResults, 2) - 1));
totalTimeLow = (baseTime + (interSeqPeriodLow*size(expectedResults, 1))*numRuns - 1); % Total trial time, including inter-sequence periods and all runs, in seconds
totalTimeHigh = (baseTime + (interSeqPeriodHigh*size(expectedResults, 1))*numRuns - 1); % Total trial time, including inter-sequence periods and all runs, in seconds

% Update the GUI
set(handles.trialTimeEstLow, 'String', sprintf('%0.2f', totalTimeLow/60));
set(handles.trialTimeEstHigh, 'String', sprintf('%0.2f', totalTimeHigh/60));

% --- Executes on button press in trialStart.
function trialStart_Callback(hObject, eventdata, handles)
% hObject    handle to trialStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% First, check if the trial should be run at all
% First, make sure that we're not using an explore circuit
if ~get(handles.testing, 'Value') && handles.TDTConnected
    isExploreCircuit = handles.DA.GetTargetVal('RZ2.silExplore');
    if isExploreCircuit
        warndlg('Cannot perform experiment while TDT is in silicon explore mode.', 'Cannot run experiment');
        return;
    end
end

% Check if the values in the GUI are default values.  If they are not,
% confirm with the user
if ~checkDefaults(handles)
    choice = questdlg('Some trial parameters are not default values.  Continue?','Non-Default Parameters','Yes','No','No');
    if strcmp(choice, 'No')
        % If the user does not want to continue, then exit now
        return;
    end
end

%% Prepare trial progress parameters
success = false;
handles.runningTrial = true;
% handles.hasFiredSequence = true; % List that the unit has fired a sequence, so that the sequence monitoring will work

%% Prepare GUI
% Block out the GUI except for the stop button
% Disable the GUI components that would interfere
disableGUI(handles);

% Cahnge the display as needed
handles = updateGUI(handles);

% Get the GUI ready
displayHandles = [handles.runNumberText, handles.runNumberOut, handles.sequenceNumberText, handles.sequenceNumberOut, handles.elapsedTimeText, handles.elapsedTimeOut, handles.runningEstTimeText, handles.runningEstTimeOut, handles.stopTrial, handles.pauseTrial];
set(displayHandles, 'visible', 'on');

%% Gather the meta-data
doAdaptive = logical(get(handles.useAdaptive, 'Value')); % Should the adaptive algorithm be used
numStimsPerSequence = str2double(get(handles.stimsPerSequence, 'String'));
numRuns = str2double(get(handles.numRuns, 'String'));
interStimPeriod = str2double(get(handles.stimPeriods, 'String'));
interSeqPeriodLow = str2double(get(handles.seqPeriodLow, 'String'));
interSeqPeriodHigh = str2double(get(handles.seqPeriodHigh, 'String'));
expectedResultsRaw = handles.curExpectedResults; % Every value that is "true" in this matrix corresponds to a spot when unit A will be stimulated
runOrders = handles.curRunOrders; % The order of trials for each run
runOrdersT = runOrders';
runOrdersCol = runOrdersT(:); % Columnified run orders for all sequences

numSequencesPerRun = size(expectedResultsRaw, 1);
numSequencesTotal = numRuns*numSequencesPerRun; % Number of sequences, total
numStimsPerRun = numSequencesPerRun*numStimsPerSequence;
trialNum = handles.curTrialNum;

allG = [str2double(get(handles.controlStrength1, 'String')) str2double(get(handles.controlStrength2, 'String'))];
allT = [str2double(get(handles.controlDuration1, 'String')) str2double(get(handles.controlDuration2, 'String'))];
handles.curExpGT = [allG' allT']; % Num_neurons x [G T]

% Get the randomized trial order for this run
expectedResults = expectedResultsRaw(runOrdersCol,:); % Get the sequence order from the runOrders matrix

% Prepare TDT for the first run/sequence
% thisRunSeqInds = 1:numSequencesPerRun; % The numbers of the sequences in the first run
% thisRunExpectedResults = expectedResults(thisRunSeqInds,:);

% Send this data to TDT
sendExperimentDataToTDT(handles, numRuns, numSequencesPerRun, numStimsPerSequence, interSeqPeriodLow, interSeqPeriodHigh, interStimPeriod);
sendSequenceOrderToTDT(handles, expectedResults); % thisRunExpectedResults); % Send the entire expected results matrix to TDT
handles = sendGTDataToTDT(handles, handles.curControlVoltages, allT);
resetStoredSequences(handles); % Tell TDT to reset the number of stored sequences
TDTPenetNum = getPenetrationNumFromTDT(handles);

% Record all available meta-data in the trailData structure of the
% neuronStruct
trialData1 = struct('numRunsPerTrial', numRuns, 'numStimsPerSequence', numStimsPerSequence, 'numSequencesPerRun', numSequencesPerRun, 'expectedResults', expectedResultsRaw, 'interStimPeriod', interStimPeriod, 'runOrders', runOrders, 'stimTimes', [], 'stimInds', [], 'opt', [], 'blockStartInfo', [], 'TDTPenetNum', TDTPenetNum, 'trialCompleted', false);
trialData2 = trialData1; % Copy the trial data from unit A, with one modification:
trialData2.expectedResults = ~expectedResultsRaw; % Get the negation of the expected results
handles.allNeuronsStruct{1}.trialData(trialNum) = trialData1; % Add the trialData structure, in position "trialNum"
handles.allNeuronsStruct{2}.trialData(trialNum) = trialData2; % Add the trialData structure, in position "trialNum"

% Update handles
guidata(hObject, handles);

% Prepare for the trial loop
quitting = false; % Is the user currently quitting out of the trial?
stimTimes = nan(numSequencesTotal, numStimsPerSequence); % A matrix to store the start times for each of the stimulations
stimInds = nan(numSequencesTotal, numStimsPerSequence, handles.numNeurons); % Store the indices in each neuronStruct's resultsData that match the stimulations in the trial

handles.seqNum = 0; % Sequence number of the sequence that is being performed
handles.runNum = 0;
handles.blockNum = 0;

% TODO (later): Store information on interpolation function and/or underlying data
% (range of alphas/sigmas/strengths/durations that were calculated, etc.

% Set up adaptation blocking structure
if doAdaptive
    adaptingUseEpochs = logical(get(handles.useAdaptingEpochs, 'Value'));
    adaptingRunsPerEpoch = str2double(get(handles.adaptingRunsPerEpoch, 'String'));
%     lastPrepareAdaptiveSeq = 0;
%     numSequencesSinceUpdate = 0; % Number of sequences that have been administered since the last adaptation update
    
    if handles.doParallelStimAndCompute
        % A GT pair is be needed to send on the first update
       handles.lastCalculatedG = handles.curControlVoltages;
       handles.lastCalculatedT = allT;
    end
    
    % Precalculate the behaviors for each block (when each GT pair gets
    % calculated and when they get sent)
    numBlocksPerRun = (numSequencesPerRun/handles.adaptiveBufferReplacementSize); % Number of blocks per run
    if adaptingUseEpochs
        numBlocksPerEpoch = numBlocksPerRun*adaptingRunsPerEpoch;
    else
        numBlocksPerEpoch = numBlocksPerRun*numRuns; % Do not use any epochs
    end
    numBlocks = numRuns*numBlocksPerRun; % Calculate the total number of blocks being dealt with
    
    % Determine which blocks require a new GT to be sent before they start
    allDoSendGTPre = repmat([true(numBlocksPerEpoch + 1, 1);false(numBlocksPerEpoch - 1, 1)], ceil(numBlocks/(2*numBlocksPerEpoch)), 1); % GT pairs should be sent before every block except those in the no-adaptation epochs (however, the first block in a no-adaptation epoch will get a new GT pair).  I know it's complicated...this makes much more sense if you diagram it out, I promise
    allDoSendGT = allDoSendGTPre(1:numBlocks);
%     allDoSendGT(1) = false; % No GT pairs should be sent on the first block, we'll just use the ones we calculated after characterization
    
    % Determine which blocks require a new GT pair to be calculating while
    % TDT runs the stimulations (to be sent to TDT right before the next
    % block)
    if handles.doParallelStimAndCompute
        allDoCalcGT = circshift(allDoSendGT, [-1 0]); % If a new GT is going to be sent before the NEXT block, then calculate a new one THIS block
        allDoCalcGT(end) = false; % No need to calculate GT during the last block
    else
        allDoCalcGTPre = repmat([true(numBlocksPerEpoch, 1);false(numBlocksPerEpoch, 1)], ceil(numBlocks/(2*numBlocksPerEpoch)), 1);
        allDoCalcGT = allDoCalcGTPre(1:numBlocks);
    end
    
    allDoCalcGT(1) = false; % No GT pair should be calculated on the first block, there is no data to calculate a new GT based on
    
    % Calculate which runs are adaptive or not
    if adaptingUseEpochs
        allRunsAdaptive = repelem(repmat([true;false], ceil(numRuns/adaptingRunsPerEpoch), 1), adaptingRunsPerEpoch);
        allRunsAdaptive = allRunsAdaptive(1:numRuns);
    else
        allRunsAdaptive = ones(numRuns, 1);
    end
    
    allSequencesAdaptive = repelem(allRunsAdaptive, numSequencesPerRun);
    thisRunAdapting = allRunsAdaptive(1);
    blockStartInfo = [1 thisRunAdapting]; % Initialize this array, which keeps track of which stimulation is the first of a new adaptive block (of either adapting or not adapting)
    
    % To indicate the first run of the control loop
    firstLoop = true;
else
    thisRunAdapting = false;
    blockStartInfo = [];
end

% Let TDT know that we're using or not using adaptation
setTDTAdapting(handles, doAdaptive);

% Start the sequence
startTrial(handles); % Send a signal to begin the trial in the TDT system
% startTime = clock; % Get the start time of the trial (only used for updating the user)

%% Start the experiment
% Start the first sequence
handles.trialStartTime = startNextRun(handles);
% TODO TDT: Get a time-stamp from TDT for the start of the first stimulation.
% Save this timestamp, then use it to calculate the stimTimes for the
% experiment

% Prepare for main loop
handles.lastSequenceAnalyzed = 0; % The sequence number of the last sequence that was analyzed
% numRecordedSequences = 0; % The number of sequences recorded
lastNumStoredSequences = 0;
lastBlockEnded = 0;
% Execute control loop
while true
    % Continue to loop through this process until the experiment ends, as
    % indicated by TDT
    
    %% Extract meta-data from TDT
    handles.seqNum = getTDTSeqNum(handles) + 1; % Sequence number of the sequence that is being performed
    handles.runNum = getTDTRunNum(handles) + 1; % Run number of the run that is being performed (or just was, if endOfRun is true)
    handles.blockNum = getTDTBlockNum(handles) + 1; % Block number of the block that is being performed (or just was, if endOfBlock is true)
    endOfRun = isTDTRunDone(handles); % Check if a TDT run has finished since the last time we checked
    endOfExperiment = isTDTExperimentDone(handles);
    endOfBlock = isTDTBlockDone(handles); % Check if a TDT block has finished since the last time we checked
    
    %     if (endOfRun || endOfExperiment || endOfBlock)
    %         pause(2);
    %     end
    
    numStoredSequences = tdtNumStoredSequences(handles);
    %     newData = numStoredSequences > 0; % Has TDT finished producing some data since the last time some was extracted?
    
%     numRecordedSequences = numRecordedSequences + numStoredSequences;
    %     Breakpoint if: ((endOfRun || endOfExperiment) && mod(numRecordedSequences, numSequencesPerRun) ~= 0)
    %% Perform data processing
    % Continuing to adapt (or not adapt) while waiting for the end of
    % the current run
    if lastNumStoredSequences ~= numStoredSequences
        %% Update the GUI
        seqNumInRun = mod(handles.seqNum - 1, numSequencesPerRun) + 1;
        elapsedTimeInSeconds = etime(clock, handles.trialStartTime);
        totalRunSequences = numSequencesPerRun*numRuns;
        curSequenceTotal = (handles.runNum - 1)*numSequencesPerRun + handles.seqNum - 1;
        
        set(handles.runNumberOut, 'String', sprintf('%d out of %d', handles.runNum, numRuns));
        set(handles.sequenceNumberOut, 'String', sprintf('%d out of %d', seqNumInRun, numSequencesPerRun));
        set(handles.elapsedTimeOut, 'String', sprintf('%0.2f minutes', elapsedTimeInSeconds/60));
        set(handles.runningEstTimeOut, 'String', sprintf('%0.2f minutes', (elapsedTimeInSeconds/60)*(totalRunSequences/curSequenceTotal)));
        
        drawnow;
    end
    
    lastNumStoredSequences = numStoredSequences;
    
    if (endOfBlock || endOfRun || endOfExperiment) && ~firstLoop
        %% Extract data from TDT
        handles.testingCurSequence = runOrdersCol(handles.seqNum); % Store the current sequence number (used in getTDTSequenceData for testing conditions)
        
        % Retreive data from TDT
        [electrodeData, electrodeDataRaw, sortCodes, interSequenceIntervals, allControlVoltages, allControlTs] = getTDTSequenceData(handles, numStimsPerSequence); % TODO SILICON: may need to refactor for silicon probes
        numNewSequences = size(electrodeData,1); % Number of new sequences that have been retreived
        
        %% Process and store data from TDT
        % Determine which sequences need to be analyzed
        analyzeSequences = handles.lastSequenceAnalyzed + (1:numNewSequences); % Analyze all sequences from the one after the last one analyzed, up until to the sequence before the current one (which is not finished yet)
        numSeqs = length(analyzeSequences);
        
        % Loop through all of the sequences that were extracted
        for seqInd = 1:numSeqs
            % TODO: Extract the electrodeData and sort codes
            thisElectrodeData = electrodeData(seqInd, :);
            thisElectrodeDataRaw = electrodeDataRaw(seqInd, :);
            thisSortCodes = sortCodes(seqInd, :);
            thisSeqStartTime = nextSeqTime;
            
            % Get the strengths and durations used for this sequence
%             controlGs = allControlVoltages(seqInd, :);
            controlVoltages = allControlVoltages(seqInd, :); % Get the control voltages for this sequence
            controlGs = zeros(size(controlVoltages));
            for gInd = 1:size(controlVoltages, 2)
               controlGs(gInd) = handles.fIrrFromControlV(controlVoltages(gInd), handles.curLaserParameters); % Calculate the laser power at each of the control voltages
            end
            controlTs = allControlTs(seqInd, :);
            
            % Prepare indexing of expectedResults and other variables
            thisSeqNum = analyzeSequences(seqInd); % The current sequence being analyzed
            thisER = expectedResults(thisSeqNum, :); % The expected results of this sequence (in binary form, i.e. true=N1, false=N2)
            thisERN = expectedResults2NeuronNums(thisER); % The expected results of this sequence (in integer form, i.e. 1=N1, 2=N2)
            if doAdaptive
                thisSequenceAdaptive = allSequencesAdaptive(thisSeqNum); % Is this sequence an adaptive sequence?
            else
                thisSequenceAdaptive = false;
            end
            
            % Set up metadata
            thisStimTimes = thisSeqStartTime + interStimPeriod*(0:(handles.numSlots - 1))/1000; % Get the stimulation times, given the interstim periods and the (approximate) start time of the sequence
            stimTimes(thisSeqNum,:) = thisStimTimes;
            thisGT = [controlGs(thisERN); controlTs(thisERN)]'; % Get the series of strength/durations of stimulation applied (e.g. [G1 G2 G2 G1 G2; T1 T2 T2 T1 T2]')
            
            % Store the newly acquired data from TDT
            [handles, ~, indsInResultsData] = addSequenceToNeuronStruct(handles, thisElectrodeData, thisElectrodeDataRaw, thisSortCodes, thisStimTimes, thisGT, numStimsPerSequence);
            
            % Acquire remaining metadata
            stimInds(thisSeqNum, :, :) = permute(indsInResultsData, [3 1 2]);
            
            for neuronNum = 1:handles.numNeurons
                % Add the sequence information to the sequenceInfo matrix
                handles.allNeuronsStruct{neuronNum}.sequenceInfo(handles.allNeuronsStruct{neuronNum}.sequenceInfoCurLoc, :) = [stimTimes(thisSeqNum,1) (thisER == mod(neuronNum, handles.numNeurons)) controlGs(neuronNum) controlTs(neuronNum) thisSequenceAdaptive]; % Amount of time since the beginning of the session, whether or not this unit was SUPPOSED to have fired, the strength and duration used to activate it, and whether or not adaptation is on for this run
                handles.allNeuronsStruct{neuronNum}.sequenceInfoCurLoc = handles.allNeuronsStruct{neuronNum}.sequenceInfoCurLoc + 1; % Increment the sequence info location
            end
            
            % Calculate the next sequence's start time, given by the
            % intersequence interval
            nextSeqTime = thisSeqStartTime + interSequenceIntervals(seqInd);
        end
        
        % TEST
        %         fprintf('Next resultsData location is %d.\n', handles.allNeuronsStruct{1}.resultsDataCurLoc);
        
        %% Finish up
        try
            handles.lastSequenceAnalyzed = analyzeSequences(end); % Get the index of the last analyzed sequence
        catch
            handles.lastSequenceAnalyzed = handles.lastSequenceAnalyzed + numSeqs;
        end
        
        % Load the stimulations from this last run into the GUI
        handles = loadStimulationsIntoGUI(handles);
        
%         redrawSequencePSpace(handles);
        drawnow;
    end
    
    %% End of the last run that TDT performed. Perform inter-Run tidying
    if endOfBlock || endOfRun || endOfExperiment
        %% Update the GUI
        % Load the stimulations from this last run into the GUI
        handles = loadStimulationsIntoGUI(handles);
        
        %% Check pausing/stopping conditions
        if endOfExperiment
            % Now that data processing has been done, if it is the end of the experiment, then break out of this loop
            break;
        end
        
        % Auto save the data between each run
        handles = autosave(handles);
        
        newHandles = guidata(hObject); % Get the most recent version of handles
        if newHandles.pausingTrial
            % First, see if the user wants to pause before the next run begins
            playWarn(handles); % Let the user know that the experiment has been paused
            set(handles.trialSettingsPanel, 'BackgroundColor', handles.yellow);
            tdtSetPaused(handles, true); % Pause TDT
            while newHandles.pausingTrial && ~newHandles.stoppingTrial
                % While the trial is still paused (and the user hasn't decided to just stop the trial), keep updating handles
                pause(.1);
                newHandles = guidata(hObject);
            end
            
            % User decided to continue
            set(handles.trialSettingsPanel, 'BackgroundColor', handles.gray);
            
            % Unpause TDT
            tdtSetPaused(handles, false);
        end
        
        if newHandles.stoppingTrial
            % If the user is trying to stop, break out of this loop
            quitting = true;
            break;
        end
        
        % Redraw the Sequence probabilities graph with the new information
        handles.sequenceTimeseriesFig = newHandles.sequenceTimeseriesFig; % Get the figure axis to the current results window, if it exists
        redrawSequenceTimeseries(handles);
        
        guidata(hObject, handles); % Update handles
        
        %% Send new data to TDT
%         thisRunSeqInds = (handles.runNum - 1)*numSequencesPerRun + (1:numSequencesPerRun); % The numbers of the sequence in this run
%         thisRunExpectedResults = expectedResults(thisRunSeqInds,:);
        
        %% Update Adapting procedure
        if doAdaptive && adaptingUseEpochs
            %             thisRunAdapting = xor(mod(handles.runNum - 1, adaptingRunsPerEpoch) == 0, thisRunAdapting); % Toggle thisRunUpdating if adaptingRunsPerEpoch number of Runs have gone by
            thisRunAdapting = allRunsAdaptive(handles.runNum);
            blockStartInfo = [blockStartInfo; [handles.runNum*numStimsPerRun + 1 thisRunAdapting]];
            
            if thisRunAdapting
                set(handles.adaptationBlock, 'Visible', 'off'); % Let the user know that we are not adapting now
            else
                set(handles.adaptationBlock, 'Visible', 'on'); % Let the user know that we are adapting now
            end
        end
        
        %         %% Start the next run in TDT
        %         startNextRun(handles);
    end
    
    
    %% Check if the experiment should be adapting now
    if doAdaptive && ((endOfBlock && lastBlockEnded ~= handles.blockNum) || firstLoop)
        % "If doing adaptive, and this is either a) the first loop, or b)
        % the end of a block that has not been processed before"
        
        % If Matlab needs to update TDT
        thisDoCalcGT = allDoCalcGT(handles.blockNum); % See if a new GT must be calculated during the next block
        thisDoSendGT = allDoSendGT(handles.blockNum); % See if the previously calculated GT should be sent before the next block
        
        fprintf('Preparing for block %d now (doCalc = %d)\n', handles.blockNum, thisDoCalcGT);
        
        if handles.doParallelStimAndCompute
            if thisDoSendGT
                % If doing asynchronous updates and this block requires a new GT, send the last calculated GT
                handles = sendGTDataToTDT(handles, handles.lastCalculatedG, handles.lastCalculatedT);
            end
            % Tell TDT to start the next adaptation block
            nextSeqTime = tdtNextAdaptationBlock(handles);
        end
        
        if thisDoCalcGT
            % If this block needs adaptation, perform calculations for adaptation blocks
            handles = prepareAdaptive(handles);
            
            if ~handles.doParallelStimAndCompute
                % If doing synchronous updates after an adaptation calculation was performed, send the last calculated GT
                handles = sendGTDataToTDT(handles, handles.lastCalculatedG, handles.lastCalculatedT);
            end
        end
        
        if ~handles.doParallelStimAndCompute
            % Tell TDT to start the next adaptation block
            nextSeqTime = tdtNextAdaptationBlock(handles);
        end
        
%         guidata(hObject, handles); % Update handles
        
        % Prepare for next round
        firstLoop = false;
    end
    
    if endOfBlock
       lastBlockEnded = handles.blockNum; 
    end
    
    %     if doAdaptive
    %         numSequencesSinceUpdate = numSequencesSinceUpdate + numStoredSequences;
%         doUpdate = numSequencesSinceUpdate >= handles.adaptiveBufferReplacementSize; % If there are handles.adaptiveBufferReplacementSize number of new sequences in TDT, then perform an update step
%         
%         %         if (handles.runNum > 1 || handles.seqNum > handles.gracePeriod)
%         if thisRunAdapting && doUpdate && (lastPrepareAdaptiveSeq ~= handles.seqNum) && ~(get(handles.testing, 'Value') && (mod(handles.seqNum, 5) ~= 0))
%             % If using adaptive (and this is an adapting block, we are past
%             % the grace period), and at least
%             % handles.adaptiveBufferReplacementSize sequences have been
%             % administered since the last time prepareAdaptive() was run
%             % (and if in testing mode while we are on one of every 5
%             % sequences), then perform adaptive optimizer duties
%             
%             set(handles.adaptationBlock, 'Visible', 'off'); % Let the user know that we are adapting now
%             
%             % Let TDT know that we're doing adaptation
%             setTDTAdapting(handles, true);
%             
%             % Perform calculations for adaptation blocks
%             handles = prepareAdaptive(handles);
%             
%             % Prepare for next round
%             numSequencesSinceUpdate = 0;
%             lastPrepareAdaptiveSeq = handles.seqNum; % Record the run number that this adaptive calculation was performed on
%         else
%             if ~(doAdaptive && thisRunAdapting)
%                 % If not using adaptive, wait until TDT has finished this sequence
%                 set(handles.adaptationBlock, 'Visible', 'on'); % Let the user know that we are not adapting now
%                 
%                 % Let TDT know that we're NOT doing adaptation
%                 setTDTAdapting(handles, false);
%             end
%             
%             %             while ~isTDTSequenceDone(handles)
%             %                 % Wait for the sequence to be done
%             %             end
%         end
%         %         else
%         %             set(handles.adaptationBlock, 'Visible', 'on'); % Let the user know that we are not adapting now
%         %         end
%     end
    
end

if quitting
    % If the user quits out of the trial, tell TDT that the trial is over
    tdtSendTrialDone(handles);
else
    success = true; % If the user did not quit, the trial ended successfully
    handles = updateTrialMatrices(handles); % If the trial ended successfully, scramble the sequence order after
end

% fprintf('\nStored %d sequences.\n\n', numRecordedSequences);

handles.stoppingTrial = false;
handles.pausingTrial = false;
handles.curTrialNum = handles.curTrialNum + 1; % Interate the trial number

% Record if the trial succeeded or failed, and add the final data
handles.allNeuronsStruct{1}.trialData(trialNum).trialCompleted = success; % Record whether or not the trial was successfully completed
handles.allNeuronsStruct{1}.trialData(trialNum).stimTimes = stimTimes; % Record the times at which the stimulations occurred
handles.allNeuronsStruct{1}.trialData(trialNum).stimInds = stimInds(:, :, 1); % Record the indices in resultsData which represent the stimulations in the trial
handles.allNeuronsStruct{1}.trialData(trialNum).blockStartInfo = blockStartInfo; % Record when each adaptation block starts, and what type of block (adapting or static) it is
handles.allNeuronsStruct{2}.trialData(trialNum).trialCompleted = success; % Record whether or not the trial was successfully completed
handles.allNeuronsStruct{2}.trialData(trialNum).stimTimes = stimTimes; % Record the times at which the stimulations occurred
handles.allNeuronsStruct{2}.trialData(trialNum).stimInds = stimInds(:, :, 2); % Record the indices in resultsData which represent the stimulations in the trial
handles.allNeuronsStruct{2}.trialData(trialNum).blockStartInfo = blockStartInfo; % Record when each adaptation block starts, and what type of block (adapting or static) it is

% % Prepare adaptive algorithm for next time
% handles = resetAdaptive(handles);

handles.runningTrial = false;

% Auto save
handles = autosave(handles);

guidata(hObject, handles); % Update handles

% Notify the user that the sequence has ended
playWarn(handles);

% Re-enable the GUI
updateGUI(handles);
set(handles.pauseTrial, 'String', 'Pause');
set(handles.adaptationBlock, 'Visible', 'off');

function neuronNums = expectedResults2NeuronNums(expectedResults)
% This function converts any sized matrix of data in the form of
% expectedResults (binary data that represents when a stimulation is given
% attempting to control a neuron, where true=N1 and false=N2) to an equal
% sized matrix of integers (where 1=N1 and 2=N2).
neuronNums = ~expectedResults + 1;


function handles = prepareAdaptive(handles)
% This function is run before each sequence in an experiment, and prepares
% the TDT system for the next sequence when using the adaptive filter.

% if handles.doParallelStimAndCompute
%     % If the GT updates are asynchronous...
%     % Send last GT pair, then begin computing the next GT pair
%     sendGTDataToTDT(handles, handles.lastCalculatedG, handles.lastCalculatedT);
%     tdtNextAdaptationBlock(handles); % Tell TDT to start the next adaptation block
% end

% Get the indicies in both neurons' resultsData of the last
% numSequencesToFit sequences of data
loc1 = handles.allNeuronsStruct{1}.resultsDataCurLoc;
loc2 = handles.allNeuronsStruct{2}.resultsDataCurLoc;

doNeuronOpt = loc1 > 1 && loc2 > 1; % Do optimization if more than 1 stimulation has been given to both neurons

%% Generate a theta for both neurons with this data (this part takes a WHILE)
if doNeuronOpt
    tic;
    handles = calculateOptimalThetas(handles, true, etime(handles.trialStartTime, handles.startTime), true);
end

fprintf('New strengths/durations chosen (%0.2f seconds):\n\nNeuron 1:\nG = %0.2f\nT = %0.2f\n\nNeuron 2:\nG = %0.2f\nT = %0.2f\n', toc, handles.curExpGT(1, 1), handles.curExpGT(1, 2), handles.curExpGT(2, 1), handles.curExpGT(2, 2));

%% Pass the GT's to TDT
% Send the next sequence of mesh points to TDT
% if ~get(handles.testing, 'Value')
%     if handles.TDTConnected
% This is a real experiment, and we are connected to TDT

% Prepare the GT values
handles.lastCalculatedG = [handles.fControlVFromIrr(handles.curExpGT(1, 1), handles.curLaserParameters) handles.fControlVFromIrr(handles.curExpGT(2, 1), handles.curLaserParameters)];
handles.lastCalculatedT = handles.curExpGT(:, 2);

%         if handles.doParallelStimAndCompute
%             % If the TDT udpates are asynchronous, then store the most
%             % recently calculated GT pairs until the TDT system wants them
%             handles.lastCalculatedG = allG;
%             handles.lastCalculatedT = allT;
%         else
%             % If the updates are synchronous, the experiment has been waiting
%             % until this finished to continue, so send the GT data
%             sendGTDataToTDT(handles, allG, allT);
%             tdtNextAdaptationBlock(handles); % Tell TDT to start the next adaptation block
%         end
%     end
% end

% Return the updated adaptive structure back to handles
% handles.adaptiveStruct = handles;

function handles = calculateOptimalThetas(handles, isAdaptive, varargin)
% This function will calculate the optimal models bestTheta1 and bestTheta2
% to fit to results data H1 and H2, given the stimulations parameterized by
% G and T. Aside from handles, all input is columnar with the same length
    
constraintNeuronBParameters = false; % Should the parameter space for neuron B be constrained
doGlobalSearch = true; % Should the optimal GT pairs be found using a global search?
doCustomGlobal = true; % Should the custom global implementation be used (because the Matlab Global Optimization Toolbox uses so much goddamn overhead)
useOldLinInterp = false; % Testing for the old linear interpolation function

% Let the user know that we're adapting now...
set(handles.adaptiveDisplay, 'Visible', 'on');
drawnow;
fprintf('Neuron parameter optimization started...\n');

% Set up the parallel pool if needed
if handles.useParallel && isempty(gcp)
   parpool('IdleTimeout', 5*60);
end

%% Extract G, T, and H data over last numSequencesToFit sequences
numStims = handles.numSlots*handles.numSequencesToFit; % The number of stims to include, if the constraint is to be used
% Get the indicies in both neurons' resultsData of the last
% numSequencesToFit sequences of data
loc1 = handles.allNeuronsStruct{1}.resultsDataCurLoc;
loc2 = handles.allNeuronsStruct{2}.resultsDataCurLoc;

if nargin > 2
    % If there is a start time, then this is not an adaptive optimization,
    % it was initiated by the user simply to characterize the neuron
    charTime = varargin{1}; % The time of the most recent characterization.  Only data after this time will be used in the fit
    useNumStimsLimit = false; % Should the number of stimulations being optimized be limited by handles.numSequencesToFit?
    if nargin > 3
        useNumStimsLimit = varargin{2};
    end
    
    resultsAfterTime1 = handles.allNeuronsStruct{1}.resultsData(:,4) > charTime; % Indices in resultsData to use (only data from stimulations after the most recent characterization, if applicable)
    resultsAfterTime2 = handles.allNeuronsStruct{2}.resultsData(:,4) > charTime; % Indices in resultsData to use (only data from stimulations after the most recent characterization, if applicable)
    
    firstInd1 = find(resultsAfterTime1, 1);
    firstInd2 = find(resultsAfterTime2, 1);
    
    if useNumStimsLimit
       firstInd1 =  max(firstInd1, loc1 - numStims);
       firstInd2 =  max(firstInd2, loc2 - numStims);
    end
    
    lastInd1 = find(resultsAfterTime1, 1, 'last');
    lastInd2 = find(resultsAfterTime2, 1, 'last');
else
    % If there is no start time, then this is an adaptive optimization,
    % performed during the experiment
    
    % Find the first index ('max' needed in case very few stimulations have
    % been performed)
    firstInd1 = max(1, loc1 - numStims);
    firstInd2 = max(1, loc2 - numStims);
    
    lastInd1 = loc1 - 1;
    lastInd2 = loc2 - 1;
end

resultsDataIndsToUse1 = firstInd1:lastInd1;
resultsDataIndsToUse2 = firstInd2:lastInd2;

% Extract the results data
trainingResultsData1 = handles.allNeuronsStruct{1}.resultsData(resultsDataIndsToUse1,:); % The results data associated with this neuron
trainingResultsData2 = handles.allNeuronsStruct{2}.resultsData(resultsDataIndsToUse2,:); % The results data associated with this neuron

% Extract G, T, H1, and H2
G = trainingResultsData1(:, 1);
T = trainingResultsData1(:, 2);
H1 = trainingResultsData1(:, 3);
H2 = trainingResultsData2(:, 3);

%% Determine which neuron should be neuron A and which should be B
% Use the user (or previous optimization) supplied control GT's to make the
% decision.

t0 = [.1 .01 0.01]; % Initial t to test, as failsafe

if ~isempty(handles.allNeuronsStruct{1}.parameters)
    bestTheta1 = handles.allNeuronsStruct{1}.parameters(end, :);
else
    bestTheta1 = t0;
end

if ~isempty(handles.allNeuronsStruct{2}.parameters)
    bestTheta2 = handles.allNeuronsStruct{2}.parameters(end, :);
else
    bestTheta2 = t0;
end

if handles.unit1IsUnitA
    HA = H1;
    HB = H2;
    lastThetaA = bestTheta1;
    lastThetaB = bestTheta2;
    lastGTA = handles.curExpGT(1,:);
    lastGTB = handles.curExpGT(2,:);
else
    HA = H2;
    HB = H1;
    lastThetaA = bestTheta2;
    lastThetaB = bestTheta1;
    lastGTA = handles.curExpGT(2,:);
    lastGTB = handles.curExpGT(1,:);
end

% Initialize the optimal thetas and gts to the previous ones
bestThetaA = lastThetaA;
bestThetaB = lastThetaB;
bestGTA = lastGTA;
bestGTB = lastGTB;

if handles.useOptServe && handles.optServeConnected 
    % Calculate the optimization using the server
    
    %% Send the data
    % Package the data to be sent
    GTH = [G T HA HB];
    lastThetas = [lastThetaA;lastThetaB];
    lastGTs = [lastGTA;lastGTB];
    
    if ~isempty(GTH)
        % Clear the buffer from the server of all other data (just in case)
        if handles.optServeT.BytesAvailable > 0
            fread(handles.optServeT, handles.optServeT.BytesAvailable);
        end
        
        % Send all of the data to the server (the order is very important!)
        tcpMatWrite(handles.optServeT, GTH);
        tcpMatWrite(handles.optServeT, lastThetas);
        tcpMatWrite(handles.optServeT, lastGTs);
        
        %% Wait while the server calculates
        while strcmp(handles.optServeT.Status, 'open') && handles.optServeT.BytesAvailable == 0
        end
        
        % Check for error
        if strcmp(handles.optServeT.Status, 'closed')
            % If the connection closed for some reason
            error('Connection to optimization server was interrupted.');
        end
        
        %% Retreive the data
        % Get the raw data
        newThetas = tcpMatRead(handles.optServeT, true);
        newGTs = tcpMatRead(handles.optServeT, true);
        calcTime = tcpMatRead(handles.optServeT, true);
        
        % Rearrange the data for storage
        bestThetaA = newThetas(1, :);
        bestThetaB = newThetas(2, :);
        
        bestGTA = newGTs(1, :);
        bestGTB = newGTs(2, :);
        
    else
        calcTime = [0 0];
    end
    
    t_optTheta = calcTime(1);
    t_optGT = calcTime(2);
else
    
    %% Prepare the optimizations
    % Convert GT to a cell (speeds up the optimization, because interpolation
    % can only be done on individual points rather than over vectors, so
    % cellfun must be used JK NOT ANYMORE CUZ SAM FIXED IT)
    
    if useOldLinInterp
        GT = mat2cell([G T], ones(size(G, 1), 1), 2);
    else
        GT = [G T];
    end
    
    % Get the initial conditions and bounds for the cost optimizing
    % strength/duration for each neuron
    startGTA = [30 2];
    startGTB = [10 15];
    GTLB = [0 0];
    GTUB = [max(handles.allStrs) max(handles.allDurs)];
    
    % t0 = [.1 .125/handles.irrRatio .01]; % Initial t to test
    LB = [0 0 .002]; % Lower bounds for t in the optimization
    UB = [max(handles.allAlphas) min(handles.b0*(max(handles.allStrs)/max(G)), handles.maxAllowableIrradiance) max(handles.allSigmas)]; % TODO: Beta restriction may be too severe TODO_POW: Multiply by strength scaling factor?  Ensure that all optimization fits with new power-based model
    
    
    % Get the error and approximation functions to be used from handles
    fErr = @FPErrorCalc;
    allStrs = handles.allStrs;
    allDurs = handles.allDurs;
    allAlphas = handles.allAlphas;
    allSigmas = handles.allSigmas;
    PAll = handles.PAll;
    b0 = handles.b0;
    approxFun = @(G, T, t) lininterpn(allStrs, allDurs, allAlphas, allSigmas, PAll, G*(t(2)/b0), T, t(1), t(3));
    
    if useOldLinInterp
        fApprox = @(GT, t) cellfun(@(GT) lininterpn(allStrs, allDurs, allAlphas, allSigmas, PAll, GT(1)*(t(2)/b0), GT(2), t(1), t(3)), GT);
    else
        fApprox = @(GT, t) approxFun(GT(:, 1), GT(:, 2), t); % lininterpn(handles.allStrs, handles.allDurs, handles.allAlphas, handles.allSigmas, handles.PAll, GT(1)*(t(2)/handles.b0), GT(2), t(1), t(3));
    end
    
    % Prepare the cost optimization function
    costFunLambda = 1e-5; %.01;
    % handles.costFun = @(P1, P2, GT) -(P1 * (1 - P2)) + costFunLambda*(GT(1).^2)*GT(2); % Cost function to optimize using control
    handles.costFun = @(P1, P2, GT) -(P1 * (1 - P2)) + costFunLambda*(GT(1).^2); % Cost function to optimize using control
    optCostFun = @(GT, t1, t2) handles.costFun(approxFun(GT(1), GT(2), t1), approxFun(GT(1), GT(2), t2), GT);
    
    algorithm = 'sqp'; % 'interior-point; % Changed to sqp on 5/16/2019
    thetaOpt  = optimoptions('fmincon', 'Algorithm', algorithm, 'Display', 'none'); % Optimization options;
    GTOpt  = optimoptions('fmincon', 'Algorithm', algorithm, 'Display', 'none'); % Optimization options;
    
    %% Calculate the optimal theta for both neurons
    tic;
    if doGlobalSearch
        if doCustomGlobal
            % Create a starting point set for the theta calculation
            %         n = 3;
            %         alphaPts = linspace(LB(1), UB(1), n); % [.05 .2];
            %         betaPts = linspace(LB(2), UB(2), n); % [.0005 .01];
            %         sigmaPts = linspace(LB(3), UB(3), n); % [.01 .1];
            
            alphaPts = linspace(LB(1), UB(1), 4); % [.05 .2];
            betaPts = linspace(LB(2), UB(2), 4); % [.0005 .01];
            sigmaPts = linspace(LB(3), UB(3), 3); % [.0005 .01];
            alphaPts = alphaPts(2:3);
            betaPts = betaPts(2:3);
            sigmaPts = sigmaPts(round(end/2));
            
            [aX, bX, sX] = ndgrid(alphaPts, betaPts, sigmaPts);
            thetaX = [aX(:) bX(:) sX(:)];
            %         thetaX = t0;
            %         t0StartPoints = CustomStartPointSet(thetaX);
            
            % Create a starting point set for the GT calculation
            %         n = 2;
            %         gPts = linspace(5, 50, n); %[5 20 50];
            %         tPts = linspace(2, 15, n); %[2 7 15];
            
            gPts = linspace(0, GTUB(1), 4); %[5 20 50];
            tPts = linspace(0, GTUB(2), 4); %[2 7 15];
            gPts = gPts(2:3);
            tPts = tPts(2:3);
            [gX, tX] = ndgrid(gPts, tPts);
            gtX = [gX(:) tX(:)];
            %         GTStartPoints = CustomStartPointSet(gtX);
        end
        
        if handles.useParallel
            gs = MultiStart('UseParallel', handles.useParallel);
        else
            gs = GlobalSearch;
        end
    end
    
    % Calculate the optimal t for the first neuron
%     bestThetaA = [];
%     bestThetaB = [];
    bestGT1 = handles.curExpGT(1, :);
    bestGT2 = handles.curExpGT(2, :);
    % try
    if doGlobalSearch
        if doCustomGlobal
            thisThetaX = [thetaX;lastThetaA]; % Include the last theta of this neuron as a starting point
            bestThetaA = multipleInitialPoint(@(x) fmincon(@(t) fErr(fApprox(GT, t), HA), x, [], [], [], [], LB, UB, [], thetaOpt), thisThetaX, handles.useParallel);
        else
            costOptimProblem = createOptimProblem('fmincon', 'objective', @(t) fErr(fApprox(GT, t), HA), 'x0', t0, 'lb', LB, 'ub', UB);
            if handles.useParallel
                bestThetaA = run(gs, costOptimProblem, t0StartPoints);
            else
                bestThetaA = run(gs, costOptimProblem);
            end
        end
    else
        bestThetaA = fmincon(@(t) fErr(fApprox(GT, t), HA), t0, [], [], [], [], LB, UB, [], thetaOpt);
    end
    % catch e
    %     throw(e);
    % end
    
    % Calculate the optimal t for the first neuron
    try
        if constraintNeuronBParameters
            % Limit the non-target alpha/beta pair to satisfy the sufficient and
            % necessary conditions (must be changed if changing the order of
            % calculation
            thisUB = [bestThetaA(1) bestThetaA(2) UB(3)];
            nonlcon = @(x) deal(x(1)/x(2) - (bestThetaA(1)/bestThetaA(2)), 0);
            if doGlobalSearch
                % Global, nonlinear constraint
                if doCustomGlobal
                    thisThetaX = [thetaX;lastThetaB]; % Include the last theta of this neuron as a starting point
                    bestThetaB = multipleInitialPoint(@(x) fmincon(@(t) fErr(fApprox(GT, t), HB), x, [], [], [], [], LB, thisUB, nonlcon, thetaOpt), thisThetaX, handles.useParallel);
                else
                    costOptimProblem = createOptimProblem('fmincon', 'objective', @(t) fErr(fApprox(GT, t), HB), 'x0', t0, 'lb', LB, 'ub', UB, 'nonlcon', nonlcon);
                    if handles.useParallel
                        bestThetaB = run(gs, costOptimProblem, t0StartPoints);
                    else
                        bestThetaB = run(gs, costOptimProblem);
                    end
                end
            else
                % Local, nonlinear constraint
                bestThetaB = fmincon(@(t) fErr(fApprox(GT, t), HB), t0, [], [], [], [], LB, thisUB, nonlcon, thetaOpt);
            end
        else
            if doGlobalSearch
                % Global, no nonlinear constraint
                if doCustomGlobal
                    thisThetaX = [thetaX;lastThetaB]; % Include the last theta of this neuron as a starting point
                    bestThetaB = multipleInitialPoint(@(x) fmincon(@(t) fErr(fApprox(GT, t), HB), x, [], [], [], [], LB, UB, [], thetaOpt), thisThetaX, handles.useParallel);
                else
                    costOptimProblem = createOptimProblem('fmincon', 'objective', @(t) fErr(fApprox(GT, t), HB), 'x0', t0, 'lb', LB, 'ub', UB);
                    if handles.useParallel
                        bestThetaB = run(gs, costOptimProblem, t0StartPoints);
                    else
                        bestThetaB = run(gs, costOptimProblem);
                    end
                end
            else
                % Local, no nonlinear constraint
                bestThetaB = fmincon(@(t) fErr(fApprox(GT, t), HB), t0, [], [], [], [], LB, UB, [], thetaOpt);
            end
        end
    end
    
    t_optTheta = toc;
    
    %% Use the cost function to generate the optimal GT to stimulate both neurons
    % Generate optimal GT pairs for both neurons
    
    % If optimization should be done, then do so.  Otherwise, keep the
    % previous value in the handles.curExpGT matrix
    fprintf('Neuron 1 Cost Optimization started...\n');
    bestGTA = [];
    bestGTB = [];
    % try
    % Try to perform the optimization.  If the optimization fails, then
    % simply do not change the experimental stimulus
    if doGlobalSearch
        %     gtX = startGTA;
        if doCustomGlobal
            thisGTX = [gtX;lastGTA]; % Include the last GT of this neuron as a starting point
            bestGTA = multipleInitialPoint(@(x) fmincon(@(T) optCostFun(T, bestThetaA, bestThetaB), x, [], [], [], [], GTLB, GTUB, [], GTOpt), thisGTX, false); % handles.useParallel);
        else
            costOptimProblem = createOptimProblem('fmincon', 'objective', @(T) optCostFun(T, bestThetaA, bestThetaB), 'x0', startGTA, 'lb', GTLB, 'ub', GTUB);
            if handles.useParallel
                bestGTA = run(gs, costOptimProblem, GTStartPoints);
            else
                bestGTA = run(gs, costOptimProblem);
            end
        end
    else
        bestGTA = fmincon(@(T) optCostFun(T, bestThetaA, bestThetaB), startGTA, [], [], [], [], GTLB, GTUB, [], GTOpt);
    end
    % catch e
    %     throw(e)
    % end
    
    fprintf('Neuron 2 Cost Optimization started...\n');
    try
        % Try to perform the optimization.  If the optimization fails, then
        % simply do not change the experimental stimulus
        if doGlobalSearch
            %         gtX = startGTB;
            if doCustomGlobal
                thisGTX = [gtX;lastGTB]; % Include the last GT of this neuron as a starting point
                bestGTB = multipleInitialPoint(@(x) fmincon(@(T) optCostFun(T, bestThetaB, bestThetaA), x, [], [], [], [], GTLB, GTUB, [], GTOpt), thisGTX, false); % handles.useParallel);
            else
                costOptimProblem = createOptimProblem('fmincon', 'objective', @(T) optCostFun(T, bestThetaB, bestThetaA), 'x0', startGTB, 'lb', GTLB, 'ub', GTUB);
                if handles.useParallel
                    bestGTB = run(gs, costOptimProblem, GTStartPoints);
                else
                    bestGTB = run(gs, costOptimProblem);
                end
            end
        else
            bestGTB = fmincon(@(T) optCostFun(T, bestThetaB, bestThetaA), startGTB, [], [], [], [], GTLB, GTUB, [], GTOpt);
        end
    catch e
        throw(e)
    end
    
    % Finalize meta-data
    t_optGT = toc - t_optTheta;
    calcTime = [t_optTheta t_optGT];
end

%% Assign calculated values to proper neurons
if handles.unit1IsUnitA
    if ~isempty(bestThetaA)
        bestTheta1 = bestThetaA;
    end
    if ~isempty(bestGTA)
        bestGT1 = bestGTA;
    end
    if ~isempty(bestThetaB)
        bestTheta2 = bestThetaB;
    end
    if ~isempty(bestGTB)
        bestGT2 = bestGTB;
    end
else
    if ~isempty(bestThetaB)
        bestTheta1 = bestThetaB;
    end
    if ~isempty(bestGTB)
        bestGT1 = bestGTB;
    end
    if ~isempty(bestThetaA)
        bestTheta2 = bestThetaA;
    end
    if ~isempty(bestGTA)
        bestGT2 = bestGTA;
    end
end

% Update which unit is considered "unit A"
handles = updateUnitCategorization(handles);

%% Save metadata from the optimization
optS.calcTime = calcTime; % Store the time to calculate theta, and each optimal stimulation GT
optS1 = optS;
optS2 = optS;

optS1.stimInds = [firstInd1 lastInd1]; % The indices used to make the calculation
optS2.stimInds = [firstInd2 lastInd2];
optS1.thetas = bestTheta1; % The theta calculated
optS2.thetas = bestTheta2;
optS1.isAdaptive = isAdaptive; % This was done outside of an adaptive trial
optS2.isAdaptive = isAdaptive;

handles.allNeuronsStruct{1}.parameters = [handles.allNeuronsStruct{1}.parameters;bestTheta1];
handles.allNeuronsStruct{2}.parameters = [handles.allNeuronsStruct{2}.parameters;bestTheta2];
handles.allNeuronsStruct{1}.optS = [handles.allNeuronsStruct{1}.optS;optS1];
handles.allNeuronsStruct{2}.optS = [handles.allNeuronsStruct{2}.optS;optS2];

% Store the GT pairs to use for control
handles.curExpGT = [bestGT1; bestGT2];

if isAdaptive
    % Store references to this specific optimization in the trialData of
    % this experiment
    optNum1 = size(handles.allNeuronsStruct{1}.optS, 1);
    optNum2 = size(handles.allNeuronsStruct{2}.optS, 1);
    trialNum = handles.curTrialNum;
    
    handles.allNeuronsStruct{1}.trialData(trialNum).opt = [handles.allNeuronsStruct{1}.trialData(trialNum).opt;optNum1];
    handles.allNeuronsStruct{2}.trialData(trialNum).opt = [handles.allNeuronsStruct{2}.trialData(trialNum).opt;optNum2];
end

% Output the alpha/beta/w values to the GUI
for neuronNum = 1:handles.numNeurons
    set(handles.(['alphaOut' num2str(neuronNum)]), 'String', handles.allNeuronsStruct{neuronNum}.parameters(end, 1));
    set(handles.(['betaOut' num2str(neuronNum)]), 'String', handles.allNeuronsStruct{neuronNum}.parameters(end, 2));
    set(handles.(['wOut' num2str(neuronNum)]), 'String', handles.allNeuronsStruct{neuronNum}.parameters(end, 3));
end

% Let the user know that we're done
set(handles.adaptiveDisplay, 'Visible', 'off');
fprintf('Theta opt: %0.2f seconds, GT opt: %0.2f seconds.\n', t_optTheta, t_optGT);
drawnow;

function TDTTrialNum = getPenetrationNumFromTDT(handles)
% CHECKED %
TDTTrialNum = 0;
if ~get(handles.testing, 'Value') && handles.TDTConnected
    % If the program is not in testing mode
    TDTTrialNum = handles.DA.GetTargetVal('RZ2.PenetrationNumber'); % Get the Trial number from TDT
end

function endOfExperiment = isTDTExperimentDone(handles)
% CHECKED %
endOfExperiment = false;
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % If the experiment has finished
        endOfExperiment = handles.DA.GetTargetVal('RZ2.TrialFullStopped'); % Find out if TDT is done with the experiment
    else
        TDTComError;
    end
else
    endOfExperiment = getTDTSeqNum(handles) == 1 && getTDTRunNum(handles) == (str2double(get(handles.numRuns, 'String')) + 1);
end

function endOfRun = isTDTRunDone(handles)
% CHECKED %
endOfRun = false;
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % If the most recent run has finished
        endOfRun = handles.DA.GetTargetVal('RZ2.RunHasCompleted'); % Find out if TDT is done with the last run
        
        if endOfRun
            handles.DA.SetTargetVal('RZ2.RunCompletionChecked', 1);
            handles.DA.ZeroTarget('RZ2.RunCompletionChecked');
        end
    else
        TDTComError;
    end
else
    endOfRun = handles.seqNum == 1 && handles.runNum ~= 1;
end

function endOfBlock = isTDTBlockDone(handles)
% CHECKED %
endOfBlock = false;
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % If the most recent block has finished
        endOfBlock = handles.DA.GetTargetVal('RZ2.BlockHasCompleted'); % Find out if TDT is done with the last block
        
        if endOfBlock
            handles.DA.SetTargetVal('RZ2.BlockCompletionChecked', 1);
            handles.DA.ZeroTarget('RZ2.BlockCompletionChecked');
        end
    else
        TDTComError;
    end
else
    endOfBlock = (handles.seqNum == 1 || handles.seqNum == 11);
end

function setTDTAdapting(handles, isAdapting)
% Let the TDT system know that we are doing adaptation now, so it should
% pause after reaching the adaptation buffer refill size
if ~get(handles.testing, 'Value') && handles.TDTConnected
    % If the program is not in testing mode
    handles.DA.SetTargetVal('RZ2.MatlabDoingAdaptive', isAdapting); % Tell TDT that we are in adaptive mode now
end

function tdtSetPaused(handles, isPaused)
% Tell TDT that we should pause the experiment for now
if ~get(handles.testing, 'Value') && handles.TDTConnected
    % If the program is not in testing mode
    handles.DA.SetTargetVal('RZ2.MatlabPaused', isPaused); % Tell TDT that we are in adaptive mode now
end

function startTime = tdtNextAdaptationBlock(handles)
% This function tells TDT that Matlab is ready for TDT to continue onto the
% next adaptation block, so it does not need to wait for Matlab anymore
if ~get(handles.testing, 'Value') && handles.TDTConnected
    % If the program is not in testing mode
    handles.DA.SetTargetVal('RZ2.NextAdaptationBlock', 1); % Tell TDT that Matlab is ready for the next adaptation block to begin
    handles.DA.ZeroTarget('RZ2.NextAdaptationBlock');
    %     startTime = getTDTTimestamp(handles);
    % else
    %     startTime = 0;
end
startTime = etime(clock, handles.startTime);

function tdtSendTrialDone(handles)
% CHECKED %
if ~get(handles.testing, 'Value') && handles.TDTConnected
    % If the program is not in testing mode
    handles.DA.SetTargetVal('RZ2.MatlabEndTrial', 1); % Tell TDT to end the trial
    handles.DA.ZeroTarget('RZ2.MatlabEndTrial');
end

function numStoredSeqs = tdtNumStoredSequences(handles)
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        numStoredSeqs = handles.DA.GetTargetVal('RZ2.NumSequencesStored');
    else
        TDTComError;
    end
else
    numStoredSeqs = (handles.runNum - 1)*size(handles.curExpectedResults, 1) + handles.seqNum - handles.lastSequenceAnalyzed - 1; % Get the number of sequences that have gone by since the last analysis
%     numStoredSeqs = double(~isTDTExperimentDone(handles));
end

function done = isTDTSequenceDone(handles)
% CHECKED %
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % If the program is not in testing mode
        done = false;
        if handles.DA.GetTargetVal('RZ2.NewSequence') && ~logical(handles.DA.GetTargetVal('RZ2.SequenceRunning'))
            % TODO TDT: Update this and the TDT circuit to see if there is new
            % data that Matlab can pull.  Matlab may need to directly reset
            % a sequence counter in TDT
            %
            % If it is a "new" sequence (this data has not already been
            % read by Matlab) and the sequence is not running, then the
            % newest sequence is done.  Send a signal to TDT to reset the
            % NewSequence value.
            done = true;
            handles.DA.SetTargetVal('RZ2.MatlabReadSequenceDone', 1); % TODO TDT: Does this need to be modified in the circuit?
            handles.DA.ZeroTarget('RZ2.MatlabReadSequenceDone');
            handles.DA.SetTargetVal('RZ2.MatResultsRead', 1);
            handles.DA.ZeroTarget('RZ2.MatResultsRead');
        end
    end
else
    done = true;
end

function blockNum = getTDTBlockNum(handles)
a = ((handles.runNum - 1)*size(handles.curExpectedResults, 1) + handles.seqNum);
b = a/handles.adaptiveBufferReplacementSize;
blockNum = floor(b); % For testing case, gets overwritten during actual experiment.

if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % Extract the block number from TDT
        blockNum = handles.DA.GetTargetVal('RZ2.blockNum'); % Find out the current run number in the TDT software (the one that is being run now [or if isTDTRunDone is true, the one that will be run soon]
    else
        TDTComError;
    end
end

function runNum = getTDTRunNum(handles)
runNum = handles.runNum - double(handles.seqNum ~= 1); % For testing case, gets overwritten during actual experiment.  Weird logic, but trust me.
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % Extract the run number from TDT
        runNum = handles.DA.GetTargetVal('RZ2.runNum'); % Find out the current run number in the TDT software (the one that is being run now [or if isTDTRunDone is true, the one that will be run soon]
    else
        TDTComError;
    end
end

function seqNum = getTDTSeqNum(handles)
seqNum = mod(handles.seqNum, size(handles.curExpectedResults, 1)); % It will be iterated outside of the function
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % Extract the sequence number from TDT
        seqNum = handles.DA.GetTargetVal('RZ2.sequenceNum'); % Find out the current sequence number in the TDT software (the one that is being run now [or if isTDTRunDone is true, the one that will be run soon]
    else
        TDTComError;
    end
end

function sendExperimentDataToTDT(handles, numRuns, numSequencesPerRun, numStimsPerSequence, interSeqPeriodLow, interSeqPeriodHigh, interStimPeriod)
% Here is where new information can be sent to TDT, used once before experiment: TODO TDT: add interStimPeriod %
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % If the program is not in testing mode
        % Send trial data
        handles.DA.SetTargetVal('RZ2.NumRunsPerTrial', numRuns); % Set the number of runs
        handles.DA.SetTargetVal('RZ2.NumSequencesPerRun', numSequencesPerRun); % Set the number of sequences per run
        handles.DA.SetTargetVal('RZ2.NumStimsPerSequence', numStimsPerSequence); % Set the number of stims per sequence
        handles.DA.SetTargetVal('RZ2.InterStimPeriod', interStimPeriod); % Set the number of stims per sequence
        sendInterSequence(handles, interSeqPeriodLow, interSeqPeriodHigh);
        
        if logical(get(handles.useAdaptive, 'Value'))
            % If adaptation is being used this trial, then send the
            % information required for the buffer
            handles.DA.SetTargetVal('RZ2.AdaptiveBufferReplacementSize', handles.adaptiveBufferReplacementSize); % Set the size of the adaptive buffer
        end
    end
end

function handles = sendGTDataToTDT(handles, allG, allT)
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % Send the GT data to TDT (must be in column vector form, one G and
        % one T for each neuron for each stimulation, to facilitate the
        % adaptive filter which may use a different GT for each
        % stimulation)
        
        % Separate out the G's for each neuron
        G1 = allG(1);
        G2 = allG(2);
        T1 = allT(1);
        T2 = allT(2);
        
        % Send the separated vectors to TDT
        handles.DA.SetTargetVal('RZ2.G1', G1);
        handles.DA.SetTargetVal('RZ2.G2', G2);
        handles.DA.SetTargetVal('RZ2.T1', T1);
        handles.DA.SetTargetVal('RZ2.T2', T2);
    end
end

% Set the GUI information
set(handles.controlStrength1, 'String', num2str(handles.fIrrFromControlV(allG(1), handles.curLaserParameters)));
set(handles.controlStrength2, 'String', num2str(handles.fIrrFromControlV(allG(2), handles.curLaserParameters)));
set(handles.controlDuration1, 'String', num2str(allT(1)));
set(handles.controlDuration2, 'String', num2str(allT(2)));
handles = setControlVoltages(handles);


function sendSequenceOrderToTDT(handles, expectedResults)
% This function sends the expectedResults to TDT
% Deal with inputs

if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % Send expected results and inter-stim periods
        % Convert to column vector, ordered by each row first (row 1, row 2,
        % etc.)
        eRTrans = expectedResults';
        handles.DA.WriteTargetVEX('RZ2.ExpectedResults', 0, 'F32', eRTrans(:)');
    end
end

function sendInterSequence(handles, interSeqPeriodLow, interSeqPeriodHigh)
% CHECKED %
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % Send the intersequence period data
        handles.DA.SetTargetVal('RZ2.InterSequencePeriodLow', interSeqPeriodLow); % Set the low value for the inter-sequence period
        handles.DA.SetTargetVal('RZ2.InterSequencePeriodHigh', interSeqPeriodHigh); % Set the high value for the inter-sequence period
    end
end

function penetrationNumberDecrement(handles)
% CHECKED %
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        % Decrement the penetration number
        handles.DA.SetTargetVal('RZ2.penNumDecrement', 1);
        handles.DA.ZeroTarget('RZ2.penNumDecrement');
    end
end

function startTrial(handles)
% CHECKED %
if ~get(handles.testing, 'Value') && handles.TDTConnected
    % If the program is not in testing mode
    handles.DA.SetTargetVal('RZ2.MatlabStartTrial', 1); % Start the sequence
    handles.DA.ZeroTarget('RZ2.MatlabStartTrial');
end

function startTime = startNextRun(handles)
% CHECKED %
if ~get(handles.testing, 'Value') && handles.TDTConnected
    % If the program is not in testing mode
    % TODO TDT: Add a MatlabStartRun hook, and start an entire run rather
    % than just a sequence
    handles.DA.SetTargetVal('RZ2.MatlabStartRun', 1); % Start the sequence
    handles.DA.ZeroTarget('RZ2.MatlabStartRun');
%     startTime = getTDTTimestamp(handles);
else
%     startTime = 0;
end
startTime = clock;
% startTime = etime(clock, handles.startTime); % The start time of the most recent sequence, in seconds

function startTime = getTDTTimestamp(handles)
% This function gets the timestamp from TDT, by retrieving the number of
% samples elapsed since starting, and converting it to time, by using the
% TDT sampling rate
startTime = 0;
if ~get(handles.testing, 'Value') && handles.TDTConnected
    startTime = tdtSamps2Sec(handles, handles.DA.GetTargetVal('RZ2.zTime'));
end

function time = tdtSamps2Sec(handles, numSamps)
% This function converts number of TDT samples to seconds
time = numSamps/handles.TDTSampRate;

function [electrodeData, electrodeDataRaw, sortCodes, interSequenceIntervals, controlGs, controlTs] = getTDTSequenceData(handles, numStimsPerSequence)
testing = get(handles.testing, 'Value');

if ~testing
    % Get TDT to pause while we read the data
    setTDTReadingData(handles, true);
    
    % Wait for TDT to be done with this sequence
    while tdtSequenceRunning(handles)
        % Do nothing...
    end
    
    % Start extracting the data...
    
    
    % First, determine how many sequences are being extracted
    numSequences = tdtNumStoredSequences(handles);
else
    % Get all run orders
    runOrders = handles.curRunOrders; % The order of trials for each run
    runOrdersT = runOrders';
    runOrdersCol = runOrdersT(:);
    
    % Get the sequences in expectedResults that will be simulated now
    numSequences = tdtNumStoredSequences(handles);
    sequenceInds = handles.lastSequenceAnalyzed + (1:numSequences);
    
    expectedResultsAll = handles.curExpectedResults(runOrdersCol, :); % The expected results for all sequences
    expectedResults = expectedResultsAll(sequenceInds, :); % The expected results for the sequences of interest
%     eRT = expectedResults';
%     eRCol = eRT(:);
    
    % Get testing data from the model
    % Gather needed information
    
    %     seqNum = handles.testingCurSequence; % Store the current sequence number (used in getTDTSequenceData for testing conditions)
    
    allG = [str2double(get(handles.controlStrength1, 'String')) str2double(get(handles.controlStrength2, 'String'))]; % Get the stimulation powers
    allT = [str2double(get(handles.controlDuration1, 'String')) str2double(get(handles.controlDuration2, 'String'))];
%     G = allG(~eRCol(sequenceInds) + 1);
%     T = allT(~eRCol(sequenceInds) + 1);
    G = allG(~expectedResults + 1);
    T = allT(~expectedResults + 1);
end

% Preallocate cells to hold all of the incoming data
electrodeData = cell(numSequences, handles.numNeurons);
electrodeDataRaw = cell(numSequences, handles.numNeurons);
sortCodes = cell(numSequences, handles.numNeurons);
interSequenceIntervals = zeros(1, numSequences);
controlGs = zeros(numSequences, 2);
controlTs = zeros(numSequences, 2);

% Store traces, raw traces, and PCA codes properly
traces1 = zeros(handles.sizeBuffer, numStimsPerSequence);
tracesRaw1 = zeros(handles.sizeBuffer, numStimsPerSequence);
traces2 = zeros(handles.sizeBuffer, numStimsPerSequence);
tracesRaw2 = zeros(handles.sizeBuffer, numStimsPerSequence);
pcaCodes1 = zeros(handles.sizeBuffer, numStimsPerSequence);
pcaCodes2 = zeros(handles.sizeBuffer, numStimsPerSequence);


% If the program is not in testing mode

% TODO TDT: Move some components to time slices
% TODO TDT: Change the RecordingStimData calculation (above the iterate
% loop)?

bufSize = handles.sizeBuffer;

isSilCircuit = handles.DA.GetTargetVal('RZ2.isSil'); % Determine if the circuit is a silicon probe circuit
isTetCircuit = isSilCircuit && handles.DA.GetTargetVal('RZ2.isTet');

if ~testing
    for seqNum = 1:numSequences
        for stimNum = 1:numStimsPerSequence
            traces1(:,stimNum) = handles.DA.ReadTargetVEX(sprintf('RZ2.VoltageTraces1_%d', stimNum), (seqNum - 1)*bufSize + 1, bufSize, 'F32', 'F32');
            tracesRaw1(:,stimNum) = handles.DA.ReadTargetVEX(sprintf('RZ2.VoltageTracesRaw1_%d', stimNum), (seqNum - 1)*bufSize + 1, bufSize, 'F32', 'F32');
            traces2(:,stimNum) = handles.DA.ReadTargetVEX(sprintf('RZ2.VoltageTraces2_%d', stimNum), (seqNum - 1)*bufSize + 1, bufSize, 'F32', 'F32');
            tracesRaw2(:,stimNum) = handles.DA.ReadTargetVEX(sprintf('RZ2.VoltageTracesRaw2_%d', stimNum), (seqNum - 1)*bufSize + 1, bufSize, 'F32', 'F32');
            %                 if isTetCircuit
            %                     pcaCodes1(:,stimNum) = handles.DA.ReadTargetVEX(sprintf('RZ2.SortCodes1_%d', stimNum), (seqNum - 1)*bufSize + 1, bufSize, 'I16', 'I16'); % Retrieve the values as integers, as they will be interpretted as a 16-bit mask
            %                     pcaCodes2(:,stimNum) = handles.DA.ReadTargetVEX(sprintf('RZ2.SortCodes2_%d', stimNum), (seqNum - 1)*bufSize + 1, bufSize, 'I16', 'I16');
            %                 else
            pcaCodes1(:,stimNum) = handles.DA.ReadTargetVEX(sprintf('RZ2.SortCodes1_%d', stimNum), (seqNum - 1)*bufSize + 1, bufSize, 'F32', 'F32');
            pcaCodes2(:,stimNum) = handles.DA.ReadTargetVEX(sprintf('RZ2.SortCodes2_%d', stimNum), (seqNum - 1)*bufSize + 1, bufSize, 'F32', 'F32');
            %                 end
            %             timeStamps(stimNum) = tdtSamps2Sec(handles, handles.DA.ReadTargetVEX(sprintf('RZ2.TimeStamp_%d', stimNum), (seqNum - 1), 1, 'I16', 'F32')); % TODO TDT: Check that this is correct
        end
        electrodeData(seqNum, :) = {traces1', traces2'};
        electrodeDataRaw(seqNum, :) = {tracesRaw1', tracesRaw2'};
        
        % If this is a tetrode probe circuit, perform some final processing
        if isTetCircuit
            % For explanation of each step, see similar processing in the
            % stimulate() function
            allCodes1Row = uint16(pcaCodes1(:))'; % The data is received as a float, so interpret the value as a uint16 so that the mask makes sense
            allCodes2Row = uint16(pcaCodes2(:))';
            
            N = 16;
            allSortBits1 = bitget(repmat(allCodes1Row, N, 1), repmat((1:N)', 1, length(allCodes1Row)));
            allSortBits2 = bitget(repmat(allCodes2Row, N, 1), repmat((1:N)', 1, length(allCodes2Row)));
            
            % Find the least significant bit of each code (the
            % "first" sort code that each waveform is sorted into).
            %  We are throwing out some sorting information, but we
            %  don't really need to know EVERY sort code that each
            %  waveform falls into
            [isValid1, finalCodes1] = max(allSortBits1, [], 1);
            [isValid2, finalCodes2] = max(allSortBits2, [], 1);
            % Make sure that we don't count any zero values
            isValid1 = logical(isValid1);
            isValid2 = logical(isValid2);
            finalCodes1(~isValid1) = 0;
            finalCodes2(~isValid2) = 0;
            
            % Store the final sort codes (feed the new results column-wise
            % back into the original matrix)
            pcaCodes1(:) = double(finalCodes1);
            pcaCodes2(:) = double(finalCodes2);
        end
        
        sortCodes(seqNum, :) = {pcaCodes1', pcaCodes2'};
        
        %     if ~testing
        interSequenceIntervals(seqNum) = handles.DA.ReadTargetVEX('RZ2.InterSequenceIntervals', (seqNum - 1), 1, 'F32', 'F32')/1000; % Get this sequence's intersequence interval, in seconds
        controlGs(seqNum, :) = [handles.DA.ReadTargetVEX('RZ2.SequenceG1s', (seqNum - 1), 1, 'F32', 'F32') handles.DA.ReadTargetVEX('RZ2.SequenceG2s', (seqNum - 1), 1, 'F32', 'F32')];
        controlTs(seqNum, :) = [handles.DA.ReadTargetVEX('RZ2.SequenceT1s', (seqNum - 1), 1, 'F32', 'F32') handles.DA.ReadTargetVEX('RZ2.SequenceT2s', (seqNum - 1), 1, 'F32', 'F32')];
        %     else
        %         interSequenceIntervals(stimNum) = seqNum - 1; % Sue me.
        %         controlGs(seqNum, :) = [handles.fControlVFromIrr(handles.curExpGT(1, 1), handles.curLaserParameters) handles.fControlVFromIrr(handles.curExpGT(2, 1), handles.curLaserParameters)]';
        %         controlTs(seqNum, :) = handles.curExpGT(:, 2)';
        %     end
    end
    
    % Let TDT know that we've retreived all of the sequences
    resetStoredSequences(handles)
    
    % Unpause the run
    setTDTReadingData(handles, false);
else
    for seqNum = 1:numSequences
        for stimNum = 1:numStimsPerSequence
            [electrodeDataStim, electrodeDataRawStim, sortCodesStim] = simulatedNeuron(G(seqNum, stimNum), T(seqNum, stimNum), handles);
            traces1(:,stimNum) = electrodeDataStim{1};
            tracesRaw1(:,stimNum) = electrodeDataRawStim{1};
            traces2(:,stimNum) = electrodeDataStim{2};
            tracesRaw2(:,stimNum) = electrodeDataRawStim{2};
            pcaCodes1(:,stimNum) = sortCodesStim{1};
            pcaCodes2(:,stimNum) = sortCodesStim{2};
            
            interSequenceIntervals(stimNum) = seqNum - 1; % Sue me.
        end
        controlGs(seqNum, :) = [handles.fControlVFromIrr(handles.curExpGT(1, 1), handles.curLaserParameters) handles.fControlVFromIrr(handles.curExpGT(2, 1), handles.curLaserParameters)]';
        controlTs(seqNum, :) = handles.curExpGT(:, 2)';
        
        electrodeData(seqNum, :) = {traces1', traces2'};
        electrodeDataRaw(seqNum, :) = {tracesRaw1', tracesRaw2'};
        sortCodes(seqNum, :) = {pcaCodes1', pcaCodes2'};
    end
end

function resetStoredSequences(handles)
% This function tells TDT to reset the counter on stored sequences
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        handles.DA.SetTargetVal('RZ2.MatlabResetStoredSequences', 1);
        handles.DA.ZeroTarget('RZ2.MatlabResetStoredSequences');
    end
end

function setTDTReadingData(handles, newVal)
% Let TDT know that it needs to pause while we read some data
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        handles.DA.SetTargetVal('RZ2.MatlabReadingData', newVal);
    end
end

function isRunning = tdtSequenceRunning(handles)
% Determine if TDT is waiting after a sequence
if ~get(handles.testing, 'Value')
    if handles.TDTConnected
        isRunning = handles.DA.GetTargetVal('RZ2.SequenceRunning');
    end
end

function [runOrders, expectedResultsAll] = createSequenceOrders(numStims, numRuns)
% This function will create one run's worth of expected results data and
% the inter-stim periods that go between each stimulation.  The expected
% results will be NxS, where N is the number of sequences per trial, and S
% is the number of stims per sequence.  The interStimPeriods will be
% Nx(S-1), and its values will be in ms.
% NumStims is the number of stimulations per sequence, and stimPeriods is
% a vector containing each inter-stim period to be tested.
% tic;
% numISPeriods = length(stimPeriods); % The number of inter-stim periods
numPositives = ceil(numStims/2); % The number of stims for the given unit in this sequence

% First, create the base for the sequences
prototypeSequence = zeros(1,numStims);
prototypeSequence(1:numPositives) = 1; % Make the prototype for the sequence (some number of positives, and some number of negatives
nonUniqueSequences = prototypeSequence(perms(1:numStims)); % Create all possible sequence orders (this implies inidividuality between the different "1"s and "0s", the next line removes this)
expectedResultsAll = unique(nonUniqueSequences, 'rows'); % Get all sequences which are functionally identical
if mod(numStims, 2) ~= 0
    % If the number of stimulations is odd
    expectedResultsAll = [expectedResultsAll;~expectedResultsAll]; % If odd, add the sequences in which there are numPositives-1 positives as well (e.x. if numStims = 5, get sequences where there are 3 neuronA + 2 neuronB, and also where there are 2 neuronA + 3 neuronB)
end
numTotalSequencesPerRun = size(expectedResultsAll, 1); % The number of sequences per run (in the uniqueSequences matrix)

% % Create the interstim periods
% numInterStimSlots = numStims - 1; % Number of inter-stimulation periods in a given sequence
% numPeriodsOrders = numISPeriods^numInterStimSlots; % How many unique orderings are there of inter-stimulation periods?
% uniqueInterStimPeriods = zeros(numPeriodsOrders, numInterStimSlots); % Make room for all IS periods
% allISPeriodInds = 1:numISPeriods; % A vector that represents each index in the inter-stim periods matrix for a given sequence
% for i = 1:numInterStimSlots
%     allInds = repmat(allISPeriodInds,numISPeriods^(numInterStimSlots - i), numISPeriods^(i - 1));
%     uniqueInterStimPeriods(:, i) = stimPeriods(allInds(:));
% end
possibleDivisors = divisors(numRuns);
if length(possibleDivisors) > 1
    possibleDivisors = possibleDivisors(1:(end - 1)); % Remove self, if possible (i.e. if numRuns ~= 1)
end
numGroups = max(gcd(possibleDivisors, numTotalSequencesPerRun));
numTrialsPerGroup = numTotalSequencesPerRun/numGroups;

% Generate random trials using Latin Squares
% Determine how many squares are needed
numSquares = ceil(numRuns/numGroups);
trialsInEachGroup = zeros(numSquares*numGroups,numTrialsPerGroup); % A matrix to hold the trials represented by each group number (column)
groupOrders = zeros(numSquares*numGroups, numGroups); % The matrix to hold each group order (row represents run number, column represents the group to be used)
for squareNum = 1:numSquares
    % Prepare the trial order for each group
    thisSquareTrials = zeros(numGroups, numTrialsPerGroup);
    thisSquareTrials(:) = randperm(numTotalSequencesPerRun);
    trialsInEachGroup(((squareNum - 1)*numGroups + 1):squareNum*numGroups, :) = thisSquareTrials;
    
    % Prepare for square construction
    groupAdd = (squareNum - 1)*numGroups; % Amount to add to this loop's square, so that each number in the square represents a unique group
    
    % Create a simple Latin Square
    thisS = mod(bsxfun(@plus, 0:(numGroups - 1), (0:(numGroups - 1))'), numGroups) + 1 + groupAdd;
    
    % Permute the simple square's rows and columns
    thisS = thisS(:, randperm(numGroups));
    thisS = thisS(randperm(numGroups), :);
    
    % Add the Latin Square into the groupOrders matrix
    groupOrders(((squareNum - 1)*numGroups + 1):squareNum*numGroups, :) = thisS;
end

% Truncate the groupOrders to only use as many runs as are needed (it may
% have been overfilled if numRuns/numGroups is not an integer)
groupOrders = groupOrders(1:numRuns,:);

% Randomize the group orders (otherwise, blocks of numGroups number of
% runs will all have the same trials grouped together, albeit in
% different orders).
groupOrders = groupOrders(randperm(numRuns),:);

% Convert the groupOrders matrix into a trialOrder matrix
% First, make a matrix with with a randomized order of each group
% (randomized per run, so every time it appears in the experiment, it
% the group'ss contents are in a different order (within-row
% randomization code by "the cyclist" on www.mathworks.com forum).
randomizedGroupedTrials = trialsInEachGroup(groupOrders',:);
numRandomizedGroups = size(randomizedGroupedTrials,1);
rowIndex = repmat((1:numRandomizedGroups)',[1 numTrialsPerGroup]);
[~,randomizedColIndex] = sort(rand(numRandomizedGroups,numTrialsPerGroup),2);
randomizedGroupedTrials = randomizedGroupedTrials(sub2ind([numRandomizedGroups,numTrialsPerGroup],rowIndex,randomizedColIndex));

% Next, feed these randomized trials into a trialOrder matrix
% (numRunsxnumTrialsPerRun)
runOrders = zeros(numTotalSequencesPerRun, numRuns);
runOrders(:) = randomizedGroupedTrials';
runOrders = runOrders';


function  numRuns_Callback(hObject, eventdata, handles)
% hObject    handle to numRuns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numRuns as text
%        str2double(get(hObject,'String')) returns contents of numRuns as a double
lastNumRuns = handles.numRunsVarData(1);
handles = checkNewInput(hObject, handles, 'numRunsVarData', true);
if lastNumRuns ~= handles.numRunsVarData(1)
    % If the number of runs changed, then update the trial matrices  (takes
    % a long time)
    handles = updateTrialMatrices(handles); % Update the trial time
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function numRuns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numRuns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stimPeriods_Callback(hObject, eventdata, handles)
% hObject    handle to stimPeriods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stimPeriods as text
%        str2double(get(hObject,'String')) returns contents of stimPeriods as a double

handles = checkNewInput(hObject, handles, 'stimPeriodsVarData');
handles = updateTrialMatrices(handles); % Update the matrices
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function stimPeriods_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimPeriods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stimsPerSequence_Callback(hObject, eventdata, handles)
% hObject    handle to stimsPerSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stimsPerSequence as text
%        str2double(get(hObject,'String')) returns contents of stimsPerSequence as a double

lastNumRuns = handles.stimsPerSequenceVarData(1);
handles = checkNewInput(hObject, handles, 'stimsPerSequenceVarData', true);
if lastNumRuns ~= handles.stimsPerSequenceVarData(1)
    % If the number of runs changed, then redo the matrices
    handles = updateTrialMatrices(handles); % Update the matrices
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function stimsPerSequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimsPerSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seqPeriodLow_Callback(hObject, eventdata, handles)
% hObject    handle to seqPeriodLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seqPeriodLow as text
%        str2double(get(hObject,'String')) returns contents of seqPeriodLow as a double

handles = checkNewInput(hObject, handles, 'seqPeriodLowVarData', false);
handles = seqPeriodCallback(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function seqPeriodLow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seqPeriodLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seqPeriodHigh_Callback(hObject, eventdata, handles)
% hObject    handle to seqPeriodHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seqPeriodHigh as text
%        str2double(get(hObject,'String')) returns contents of seqPeriodHigh as a double

handles = checkNewInput(hObject, handles, 'seqPeriodHighVarData', false);
handles = seqPeriodCallback(handles);
guidata(hObject, handles);

function handles = seqPeriodCallback(handles)
% When the inter-sequence period is updated by the user
handles = updateTrialTimes(handles); % Update the trial time
interSeqPeriodLow = str2double(get(handles.seqPeriodLow, 'String'));
interSeqPeriodHigh = str2double(get(handles.seqPeriodHigh, 'String'));
sendInterSequence(handles, interSeqPeriodLow, interSeqPeriodHigh)

% --- Executes during object creation, after setting all properties.
function seqPeriodHigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seqPeriodHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stopTrial.
function stopTrial_Callback(hObject, eventdata, handles)
% hObject    handle to stopTrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% The user is attempting to stop the sequence from running
handles.stoppingTrial = true;
guidata(hObject, handles);

% Make this button invisible again
set(hObject, 'Visible', 'off');



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to stimsPerSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stimsPerSequence as text
%        str2double(get(hObject,'String')) returns contents of stimsPerSequence as a double


% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimsPerSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton68.
function pushbutton68_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton69.
function pushbutton69_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton70.
function pushbutton70_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12


% --- Executes on button press in pauseTrial.
function pauseTrial_Callback(hObject, eventdata, handles)
% hObject    handle to pauseTrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.pausingTrial
    % If the trial is currently paused, the user wants to unpause the trial
    str = 'P';
    paused = false;
else
    % If the trial is currently running, the user wants to pause the trial
    str = 'Unp';
    paused = true;
end

handles.pausingTrial = paused;
set(hObject, 'String', sprintf('%sause', str));
guidata(hObject, handles); % Update handles
drawnow;

function checkSaveBannerState(handles)
% CHECKED %
% Check if the current neuronStruct is the same as the most recently saved
% neuronStruct.  If they are different, let the user know that a new save
% is required.

if ~isempty(handles.autoSaveFiles) && ~isEqualNeuronStruct(handles.lastSavedNeuronsStruct, handles.allNeuronsStruct)
    % If a structre has been saved, but the current struct and the last
    % saved struct are not equal, then change the banner
    
    set(handles.autoSaveStatus, 'BackgroundColor', 'yellow', 'String', 'Not Saved');
    drawnow;
end

function tdtBlockName_Callback(hObject, eventdata, handles)
% hObject    handle to tdtBlockName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tdtBlockName as text
%        str2double(get(hObject,'String')) returns contents of tdtBlockName as a double
set(hObject, 'String', handles.tdtBlockString);
guidata(hObject, handles);

function tdtCharBlockName_Callback(hObject, eventdata, handles)
% hObject    handle to tdtBlockName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tdtBlockName as text
%        str2double(get(hObject,'String')) returns contents of tdtBlockName as a double
set(hObject, 'String', handles.tdtCharBlockString);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tdtBlockName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tdtBlockName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function controllerFileName_Callback(hObject, eventdata, handles)
% hObject    handle to controllerFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of controllerFileName as text
%        str2double(get(hObject,'String')) returns contents of controllerFileName as a double
set(hObject, 'String', handles.controllerSaveName);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function controllerFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to controllerFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in viewSequencePSpace.
% function viewSequencePSpace_Callback(hObject, eventdata, handles)
% % CHECKED %
% % hObject    handle to viewSequencePSpace (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% if (~isempty(handles.sequencePSpaceFig1) && ~isempty(handles.sequencePSpaceFig2)) && (ishandle(handles.sequencePSpaceFig1) && ishandle(handles.sequencePSpaceFig2))
%     % If so, just bring it to the front
%     figure(handles.sequencePSpaceFig1);
%     figure(handles.sequencePSpaceFig2);
% else
%     % If not, create a new figure and assign it to the sequencePSpaceFig
%     % variable
%     handles.sequencePSpaceFig1 = figure;
%     handles.sequencePSpaceFig2 = figure;
% end
% 
% % Now that we have a figure, redraw it
% redrawSequencePSpace(handles);
% 
% guidata(hObject, handles);

% --- Executes on button press in swapGTs.
function swapGTs_Callback(hObject, eventdata, handles)
% hObject    handle to swapGTs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the old values from the GUI
n1G = str2double(get(handles.controlStrength1, 'String'));
n2G = str2double(get(handles.controlStrength2, 'String'));
n1T = str2double(get(handles.controlDuration1, 'String'));
n2T = str2double(get(handles.controlDuration2, 'String'));

% Modify the GUI
set(handles.controlStrength1, 'String', num2str(n2G));
set(handles.controlStrength2, 'String', num2str(n1G));
set(handles.controlDuration1, 'String', num2str(n2T));
set(handles.controlDuration2, 'String', num2str(n1T));

% Update the vardata
handles.manStrength1VarData(1) = n2G;
handles.manDuration1VarData(1) = n2T;
handles.manStrength2VarData(1) = n1G;
handles.manDuration2VarData(1) = n1T;

% Update which unit is unit A and unit B
handles = updateUnitCategorization(handles);

% Update the control voltages
handles = setControlVoltages(handles);

% Update the plot
handles = drawSDPlot(handles);

% Update handles
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function viewSequenceTimeseries_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viewSequenceTimeseries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function controllerExistingFilePath_Callback(hObject, eventdata, handles)
% hObject    handle to controllerExistingFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of controllerExistingFilePath as text
%        str2double(get(hObject,'String')) returns contents of controllerExistingFilePath as a double

handles = checkControllerFileStatus(handles);
guidata(hObject, handles);

function handles = checkControllerFileStatus(handles)
% TODO SILICON: Need to update this for tetrodes? %
% This function is used to check the status of a controller, specifically
% the dual SpikePac controller. It will first check if the user-provided
% file exists and is a .csf file, and then will check if the file has been
% modified or not.
newString = get(handles.controllerExistingFilePath, 'String');
handles.controllerFileModified(handles.modCSFIndex) = false;
handles.haveControllerFilePath(handles.modCSFIndex) = false;

% Check if the filepath currently in the edit box is a proper .csf file pathway
if checkControllerFilePath(handles, newString)
    handles.controllerFilePath{handles.modCSFIndex} = newString;
    set(handles.controllerExistingFilePath, 'BackgroundColor', handles.gray);
    handles.haveControllerFilePath(handles.modCSFIndex) = true;
    
    % If the filepath is a proper .csf file, check if it has been modified
    if checkControllerFileModified(handles)
        handles.controllerFileModified(handles.modCSFIndex) = true;
        
        set(handles.csfStatus, 'String', 'File Ready (Now load into controller)', 'BackgroundColor', handles.brightGreen);
        set(handles.modifyCSF, 'Enable', 'off');
    else
        set(handles.csfStatus, 'String', 'File not yet modified', 'BackgroundColor', handles.yellow);
        set(handles.modifyCSF, 'Enable', 'on');
    end
    
else
    set(handles.controllerExistingFilePath, 'BackgroundColor', handles.red);
    set(handles.csfStatus, 'String', '', 'BackgroundColor', handles.gray);
    set(handles.modifyCSF, 'Enable', 'off');
end

% Update the GUI
handles = updateGUI(handles);

% --- Executes on button press in modifyCSF.
function modifyCSF_Callback(hObject, eventdata, handles)
% hObject    handle to modifyCSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = modifyCSF(handles);
guidata(hObject, handles); % Update handles

function handles = modifyCSF(handles)
% CHECKED %
% This function modifies a CSF file so that both SpikePAC controllers'
% parameters are identical

if ~handles.haveControllerFilePath(handles.modCSFIndex) || (handles.haveControllerFilePath(handles.modCSFIndex) && handles.controllerFileModified(handles.modCSFIndex))
    % If there is no controller file, or it exists and is already modified,
    % simply return
    return;
end

% Extract the filecell from the csf file, and modify it
fileCell = loadTextFile(handles.controllerFilePath{handles.modCSFIndex});
if handles.modCSFIndex == 1
    sourceLinesOfInterest = extractTwoElectrodeCSFLinesOfInterest(fileCell, handles.sourceController{handles.modCSFIndex}, handles.targetLineStart);
    fileCell = insertTwoElectrodeLinesOfInterest(fileCell, handles.targetController, handles.targetLineStart, sourceLinesOfInterest); % Insert the source lines of interest into the target's lines of interest area
elseif handles.modCSFIndex == 2
    sourceLinesOfInterest = extractTetrodeCSFLinesOfInterest(fileCell, handles.sourceController{handles.modCSFIndex}, handles.tetFirstLine);
    fileCell = insertTetrodeLinesOfInterest(fileCell, handles.targetController{handles.modCSFIndex}, handles.tetFirstLine, sourceLinesOfInterest); % Insert the source lines of interest into the target's lines of interest area
end

% Write a second copy of the csf file, to save for later
writeSavedCSFFile(fileCell, sprintf('%s%s%s', fileparts(handles.controllerFilePath{handles.modCSFIndex}), filesep, handles.controllerSaveName{handles.modCSFIndex}));

% Write the file to disk
writeTextFile(fileCell, handles.controllerFilePath{handles.modCSFIndex});

% Check the files status now, and update the GUI
handles = checkControllerFileStatus(handles);
handles = updateGUI(handles);

function writeSavedCSFFile(fileCell, fileName)
% CHECKED %
% Write a text file with the given name (fileName should not include the
% file extension) and a number that does not already exist in the
% directory.
ext = '.csf';
fileInd = 1;

% Find a file name that has not been used yet
while true
    fullFileName = sprintf('%s_%.2d%s', fileName, fileInd, ext);
    if ~exist(fullFileName, 'file')
        break;
    end
    fileInd = fileInd + 1;
end

% Once a file name has been found, write the file contents to this file
writeTextFile(fileCell, fullFileName);

function writeTextFile(fileCell, filePath)
% CHECKED %
% Write a text file to the given filePath, using the contents of a fileCell
fid = fopen(filePath, 'w');
formatspec = '%s\r\n';
for ind = 1:length(fileCell)
    fprintf(fid, formatspec, fileCell{ind});
end
fclose(fid);

function fileCell = insertTwoElectrodeLinesOfInterest(fileCell, controllerName, targetLineStart, linesOfInterest)
% TODO SILICON: Need to update this for tetrodes? %
% Insert the chosen lines of interest into the target controller's lines of
% interest area
inSegment = false;
targetLinesInd = 1;

if ~iscell(targetLineStart)
    targetLineStart = {targetLineStart}; % Make a cell
end
numTargetLines = length(targetLineStart);

for ind = 1:length(fileCell)
    % Go through each line, check if it's in a controller's header
    if (length(fileCell{ind}) > 3) && strcmp(fileCell{ind}(1:3), '>>>')
        % A new PCA controller header was just hit
        if (length(fileCell{ind}) >= (3 + length(controllerName))) && strcmp(fileCell{ind}, sprintf('>>>%s', controllerName))
            % We are in the target header
            inSegment = true;
        else
            if inSegment
                % If we were previously in a segment, break out of the loop
                break;
            end
        end
    end
    
    % If inside the header segment for the target object
    if inSegment
        for tgtNum = 1:numTargetLines
            if (length(fileCell{ind}) > length(targetLineStart{tgtNum})) && strcmp(fileCell{ind}(1:length(targetLineStart{tgtNum})), targetLineStart{tgtNum})
                % If this line is longer than one of the target strings,
                % and starts with it (target strings being a string of
                % characters that indicates that this line should be
                % extracted
                if isempty(strfind(fileCell{ind}, controllerName))
                    % If the controller name is not in the line (otherwise,
                    % component renaming might occur)
                    if targetLinesInd <= length(linesOfInterest)
                        fileCell{ind} = linesOfInterest{targetLinesInd};
                        targetLinesInd = targetLinesInd + 1;
                    else
                        fileCell{ind} = [];
                    end
                end
            end
        end
    end
end

function fileCell = insertTetrodeLinesOfInterest(fileCell, controllerName, targetLineStart, linesOfInterest)
% TODO SILICON: Need to update this for tetrodes? %
% Insert the chosen lines of interest into the target controller's lines of
% interest area
targetLinesInd = 1;

for ind = 1:length(fileCell)
    % Go through each line, check if it's in a controller's header
    if (length(fileCell{ind}) > 3) && strcmp(fileCell{ind}(1:3), '>>>')
        % A new PCA controller header was just hit
        if (length(fileCell{ind}) >= (3 + length(controllerName))) && strcmp(fileCell{ind}, sprintf('>>>%s', controllerName))
            % We are in the source header
            startLine = ind;
            break;
        end
    end
end

inSegment = false;
endLineBeginning = '>>>';  % If the line starts with this string, exit the loop

doEndClean = true;
for ind = startLine:length(fileCell)
    % If inside the header segment for the source object
    if (length(fileCell{ind}) > length(targetLineStart)) && strcmp(fileCell{ind}(1:length(targetLineStart)), targetLineStart)
        % If this line is longer than one of the target strings,
        % and starts with it (target strings being a string of
        % characters that indicates that this line should be
        % extracted
        inSegment = true;
    end
    
    if inSegment
        % If we're entering the next controller's segment, then insert all
        % of the remaining lines, and end the loop
        if ((length(fileCell{ind}) > length(endLineBeginning)) && strcmp(fileCell{ind}(1:length(endLineBeginning)), endLineBeginning))
            fileCell = [fileCell(1:(ind - 1)) linesOfInterest(targetLinesInd:end) fileCell(ind:end)]; % Insert the rest of the target lines in here, then finish up
            return;
        end
        
        if targetLinesInd > length(linesOfInterest)
            % If we go past the linesOfInterest cell length, then make sure
            % there is nothing else between here and the next controller
            % block/end of file
            endLine = ind;
            break;
        end
        
        fileCell{ind} = linesOfInterest{targetLinesInd};
        targetLinesInd = targetLinesInd + 1;
    end
    
    if ind == length(fileCell)
        doEndClean = false;
    end
end

if doEndClean
    % Ensure that there is nothing after this segment ends, unless a new
    % controller segment starts
    lastLine = length(fileCell);
    for ind = endLine:length(fileCell)
        if ((length(fileCell{ind}) > length(endLineBeginning)) && strcmp(fileCell{ind}(1:length(endLineBeginning)), endLineBeginning))
            lastLine = ind - 1;
            break;
        end
    end
    fileCell(endLine:lastLine) = [];
else
    if targetLinesInd <= length(linesOfInterest)
        fileCell = [fileCell linesOfInterest(targetLinesInd:end)];
    end
end


function isCorrect = checkControllerFilePath(handles, newString)
% CHECKED %
% Check that the string represents a correct file path
isCorrect = false;
if (length(newString) > length(handles.controllerFileExtension)) && strcmp(newString((end - length(handles.controllerFileExtension) + 1):end), handles.controllerFileExtension)
    % If the inputted string has the correct extension
    if exist(newString, 'file')
        % If the string represents a file that exists
        isCorrect = true;
    end
end

function isModified = checkControllerFileModified(handles)
% CHECKED %

isModified = true;

fileCell = loadTextFile(handles.controllerFilePath{handles.modCSFIndex});
if handles.modCSFIndex == 1
    sourceLinesOfInterest = extractTwoElectrodeCSFLinesOfInterest(fileCell, handles.sourceController{handles.modCSFIndex}, handles.targetLineStart);
    targetLinesOfInterest = extractTwoElectrodeCSFLinesOfInterest(fileCell, handles.targetController{handles.modCSFIndex}, handles.targetLineStart);
elseif handles.modCSFIndex == 2
    sourceLinesOfInterest = extractTetrodeCSFLinesOfInterest(fileCell, handles.sourceController{handles.modCSFIndex}, handles.tetFirstLine);
    targetLinesOfInterest = extractTetrodeCSFLinesOfInterest(fileCell, handles.targetController{handles.modCSFIndex}, handles.tetFirstLine);
end
% Check if the source lines of interest and the target lines of interest
% are the same
if length(sourceLinesOfInterest) == length(targetLinesOfInterest)
    for ind = 1:length(sourceLinesOfInterest)
        if ~strcmp(sourceLinesOfInterest{ind}, targetLinesOfInterest{ind})
            isModified = false;
            break;
        end
    end
else
    isModified = false;
end

function linesOfInterest = extractTwoElectrodeCSFLinesOfInterest(fileCell, controllerName, targetLineStart)
% This function finds the lines of interest in an existing fileCell, and
% extracts them into a new cell for easy exploration
inSegment = false;
linesOfInterest{1} = [];
sourceLinesInd = 1;
if ~iscell(targetLineStart)
    targetLineStart = {targetLineStart}; % Make a cell
end
numTargetLines = length(targetLineStart);

for ind = 1:length(fileCell)
    % Go through each line, check if it's in a controller's header
    if (length(fileCell{ind}) > 3) && strcmp(fileCell{ind}(1:3), '>>>')
        % A new PCA controller header was just hit
        if (length(fileCell{ind}) >= (3 + length(controllerName))) && strcmp(fileCell{ind}, sprintf('>>>%s', controllerName))
            % We are in the source header
            inSegment = true;
        else
            if inSegment
                % If we were previously in a segment, break out of the loop
                break;
            end
        end
    end
    
    % If inside the header segment for the source object
    if inSegment
        for tgtNum = 1:numTargetLines
            if (length(fileCell{ind}) > length(targetLineStart{tgtNum})) && strcmp(fileCell{ind}(1:length(targetLineStart{tgtNum})), targetLineStart{tgtNum})
                % If this line is longer than one of the target strings,
                % and starts with it (target strings being a string of
                % characters that indicates that this line should be
                % extracted
                if isempty(strfind(fileCell{ind}, controllerName))
                    % If the controller name is not in the line (otherwise,
                    % component renaming might occur)
                    linesOfInterest{sourceLinesInd} = fileCell{ind};
                    sourceLinesInd = sourceLinesInd + 1;
                    break; % Break out of target line num loop
                end
            end
        end
    end
end

function linesOfInterest = extractTetrodeCSFLinesOfInterest(fileCell, controllerName, targetLineStart)
% This function finds the lines of interest in an existing fileCell, and
% extracts them into a new cell for easy exploration
inSegment = false;
linesOfInterest{1} = [];
sourceLinesInd = 1;

for ind = 1:length(fileCell)
    % Go through each line, check if it's in a controller's header
    if (length(fileCell{ind}) > 3) && strcmp(fileCell{ind}(1:3), '>>>')
        % A new PCA controller header was just hit
        if (length(fileCell{ind}) >= (3 + length(controllerName))) && strcmp(fileCell{ind}, sprintf('>>>%s', controllerName))
            % We are in the source header
            startLine = ind;
            break;
        end
    end
end

inSegment = false;
endLineBeginning = '>>>';  % If the line starts with this string, exit the loop
for ind = startLine:length(fileCell)
    % If inside the header segment for the source object
    if (length(fileCell{ind}) >= length(targetLineStart)) && strcmp(fileCell{ind}(1:length(targetLineStart)), targetLineStart)
        % If this line is longer than one of the target strings,
        % and starts with it (target strings being a string of
        % characters that indicates that this line should be
        % extracted
        inSegment = true;
    end
    
    if inSegment
        % If we're entering the next controller's segment, then end the loop
        if (length(fileCell{ind}) > length(endLineBeginning)) && strcmp(fileCell{ind}(1:length(endLineBeginning)), endLineBeginning)
            break;
        end
        
        linesOfInterest{sourceLinesInd} = fileCell{ind};
        sourceLinesInd = sourceLinesInd + 1;
    end
end

% Clear any blank lines from the end
for ind = length(linesOfInterest):-1:1
    if isempty(linesOfInterest{ind})
        linesOfInterest(ind) = [];
    else
        % As soon as we hit a non-blank line, exit
       break; 
    end
end

function fileCell = loadTextFile(fileName)
% CHECKED %
% This function loads a text file into a cell, one cell per line of text
% (No error check is done in this function
fileCell{1} = [];
fid = fopen(fileName);
line = fgetl(fid);
ind = 1;
while ischar(line)
    fileCell{ind} = line;
    ind = ind + 1;
    line = fgetl(fid);
end
fclose(fid);

% --- Executes during object creation, after setting all properties.
function controllerExistingFilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to controllerExistingFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browseCSF.
function browseCSF_Callback(hObject, eventdata, handles)
% hObject    handle to browseCSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.haveControllerFilePath(handles.modCSFIndex)
    startPath = handles.controllerFilePath{handles.modCSFIndex};
else
    [possiblePath, possiblePath2] = fileparts(handles.controllerFilePath);
    if exist(possiblePath, 'dir')
        startPath = possiblePath;
    elseif exist(possiblePath2, 'dir')
        startPath = possiblePath;
    else
        startPath = pwd;
    end
end


[fileName, pathName] = uigetfile({'*.csf', 'Controller File (.csf)'}, 'Choose Controller File', startPath);

if ischar(fileName) && ischar(pathName)
    set(handles.controllerExistingFilePath, 'String', [pathName fileName]); % Set the GUI editbox to be this new file name
    handles = checkControllerFileStatus(handles); % Check the new file
    guidata(hObject, handles); % Update handles
end


% --- Executes on key press with focus on figure1 and none of its controls.
function stimOnKeyPress(eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if ~handles.runningTrial
    val=double(get(handles.figure1,'CurrentCharacter')); % Get the current key being pressed
    
    if isempty(val)
        return;
    end
    
    if val == 28 && handles.stimEnabled
        % Stim 1
        manStim1_Callback(handles.manStim1, eventdata, handles);
    elseif val == 29 && handles.stimEnabled
        % Stim 2
        manStim2_Callback(handles.manStim2, eventdata, handles);
    elseif val == 30 || val == 31
        curLoc1 = get(handles.neuron1StimsListbox, 'Value');
        curLoc2 = get(handles.neuron2StimsListbox, 'Value');
        if val == 30
            % Move up (back in history)
            set(handles.neuron1StimsListbox, 'Value', max(curLoc1 - 1, 1));
            
            set(handles.neuron2StimsListbox, 'Value', max(curLoc2 - 1, 1));
        elseif val == 31
            % Move down (forward in history)
            boxLength1 = length(get(handles.neuron1StimsListbox, 'String'));
            set(handles.neuron1StimsListbox, 'Value', min(curLoc1 + 1, boxLength1));
            
            boxLength2 = length(get(handles.neuron2StimsListbox, 'String'));
            set(handles.neuron2StimsListbox, 'Value', min(curLoc2 + 1, boxLength2));
        end
        
        handles = stimsListboxCallback(handles);
        guidata(handles.neuron1StimsListbox, handles);
    end
end




% --- Executes on button press in fullTrace1.
function fullTrace1_Callback(hObject, eventdata, handles)
% hObject    handle to fullTrace1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fullTrace1
neuronNum = get(hObject, 'UserData');
handles = drawVoltageTraces(handles, neuronNum);
guidata(hObject, handles);

% --- Executes on button press in fullTrace2.
function fullTrace2_Callback(hObject, eventdata, handles)
% hObject    handle to fullTrace2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fullTrace2
neuronNum = get(hObject, 'UserData');
handles = drawVoltageTraces(handles, neuronNum);
guidata(hObject, handles);




function [handles, results, indsInResultsData] = addSequenceToNeuronStruct(handles, electrodeData, electrodeDataRaw, sortCodes, stimTimes, allGTs, numStimsPerSequence)
% TODO: Is this the best way to add data? %
indsInResultsData = zeros(numStimsPerSequence, handles.numNeurons);
results = zeros(numStimsPerSequence, handles.numNeurons);
allG = allGTs(:,1);
allT = allGTs(:,2);
for stimNum = 1:numStimsPerSequence
    % For each stimulation in this sequence
    for neuronNum = 1:handles.numNeurons
        % For each unit, get the results from the stimulation
        if get(handles.(sprintf('usePCACheck%d', neuronNum)), 'Value')
            % If this unit should be using PCA, then determine the results
            % using the sort codes from TDT
            results(stimNum, neuronNum) = detectSpikesFromTraces(electrodeData{neuronNum}(stimNum, :), handles.traceTime, handles.allNeuronsStruct{neuronNum}.defaultTriggerSettings, sortCodes{neuronNum}(stimNum, :), 1, get(handles.testing, 'Value'));
        else
            % If this unit should not be using PCA, then determine the results
            % using a voltage threshold
            results(stimNum, neuronNum) = detectSpikesFromTraces(electrodeData{neuronNum}(stimNum, :), handles.traceTime, handles.allNeuronsStruct{neuronNum}.defaultTriggerSettings, [], 1, get(handles.testing, 'Value'));
        end
        
        % Record the index in resultsData that each of the new results are
        % going into
        indsInResultsData(stimNum, neuronNum) = handles.allNeuronsStruct{neuronNum}.resultsDataCurLoc;
    end
    %             % For testing: Randomly toggle 70% of the results
    %             if get(handles.testing, 'Value')
    %                 toggleMask = rand(size(results(stimNum, :))) > .9;
    %                 results(stimNum, :) = xor(results(stimNum, :), toggleMask);
    %             end
    
    % Add the results to the neuronStruct
    handles.allNeuronsStruct = addResultToNeuronStruct(handles.allNeuronsStruct, allG(stimNum), allT(stimNum), results(stimNum, :), cellfun(@(x) x(stimNum, :), electrodeData, 'uniformOutput', 0), cellfun(@(x) x(stimNum, :), electrodeDataRaw, 'uniformOutput', 0), stimTimes(stimNum), cellfun(@(x) x(stimNum, :), sortCodes, 'uniformOutput', 0));
end


function testNoise_Callback(hObject, eventdata, handles)
% hObject    handle to testNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of testNoise as text
%        str2double(get(hObject,'String')) returns contents of testNoise as a double
handles = checkNewInput(hObject, handles, 'testNoiseVarData');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function testNoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to testNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in useAdaptive.
function useAdaptive_Callback(hObject, eventdata, handles)
% hObject    handle to useAdaptive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useAdaptive
checkDefaults(handles);



function startStepFac_Callback(hObject, eventdata, handles)
% hObject    handle to startStepFac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startStepFac as text
%        str2double(get(hObject,'String')) returns contents of startStepFac as a double
handles = checkNewInput(hObject, handles, 'startStepFacVarData');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function startStepFac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startStepFac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = setMaxStrengthFromIrradiance(handles)
% CHECKED %
% This function finds the maximum strength that is allowable, considering
% the irradiance boundary and the laser equation values from TDT

% Calculate the maximum power in mW
maxIrr = handles.maxAllowableIrradiance; %*(pi*(handles.laserRadius^2));

% Apply the irradiance boundary to the upper G bound
if ~get(handles.testing, 'Value') && handles.TDTConnected
    %     A = handles.DA.GetTargetVal('RZ2.LasEquA');
    %     B = handles.DA.GetTargetVal('RZ2.LasEquB');
    %     C = handles.DA.GetTargetVal('RZ2.LasEquC');
    %     D = handles.DA.GetTargetVal('RZ2.LasEquD');
    %     laserEQ = [A B C D];
    
    laserEQ = handles.curLaserParameters;
    
    % TODO: If this is successful, save the laser parameters %
else
    laserEQ =  [-9.9757   75.5625  201.5847  -10.7840];
end

% Save this laserEq to the neuronStruct of both neurons
for neuronNum = 1:handles.numNeurons
    handles.allNeuronsStruct{neuronNum}.laserEQ = laserEQ;
end

% Try different starting positions for fzero to get a positive
% maxPs
handles.maxIrradianceStrength = handles.maxStrength; % Set the value as the standard maximum strength
if handles.fIrrFromControlV(handles.maxStrength, laserEQ) >= maxIrr
    % If the laser at maximum power is over the irradiance maximum,
    % then find the upper boundary for the laser strength
    maxP = 0;
    ready = false;
    for x0 = 1:handles.maxStrength
        maxP = handles.fIrrFromControlV(maxIrr, laserEQ); % Maximum power in arbitrary laser strength units
        if 0 < maxP && maxP < handles.maxStrength
            ready = true;
            break;
        end
    end
    %
    if ready
        handles.maxIrradianceStrength = maxP; % If a positive maxP was found, then use it (if it is lower than the max strength)
        
        % If no positive maxP was found, it means that, in the range of
        % strengths be used, there is no stimulation that will go above
        % the irradiance boundary.
    end
end


% --- Executes on button press in useVerification.
function useVerification_Callback(hObject, eventdata, handles)
% hObject    handle to useVerification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useVerification


% --- Executes on button press in useAdaptingEpochs.
function useAdaptingEpochs_Callback(hObject, eventdata, handles)
% hObject    handle to useAdaptingEpochs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useAdaptingEpochs
checkDefaults(handles);


function adaptingRunsPerEpoch_Callback(hObject, eventdata, handles)
% hObject    handle to adaptingRunsPerEpoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of adaptingRunsPerEpoch as text
%        str2double(get(hObject,'String')) returns contents of adaptingRunsPerEpoch as a double
handles = checkNewInput(hObject, handles, 'runsPerBlockVarData', true);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function adaptingRunsPerEpoch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adaptingRunsPerEpoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in characterizeNeurons.
function characterizeNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to characterizeNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable the GUI
disableGUI(handles);

% try
handles = characterizeNeuronCallback(handles);
% catch m
% Re-enable the GUI
%     handles = updateGUI(handles);
%     guidata(hObject, handles); % Update handles
%     rethrow(m);
% end


% Re-enable the GUI
handles = updateGUI(handles);

guidata(hObject, handles); % Update handles


% --- Executes on button press in fitSDCurve.
function fitSDCurve_Callback(hObject, eventdata, handles)
% hObject    handle to fitSDCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable the GUI
disableGUI(handles);

% neuronNum = get(hObject, 'UserData'); % Number of this neuron
handles = fitSDCurveCallback(handles);

% Re-enable the GUI
handles = updateGUI(handles);

guidata(hObject, handles); % Update handles


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Outputs from this function are returned to the command line.
function varargout = TDT_Control_OutputFcn(hObject, eventdata, handles)
% CHECKED
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Resize the figure
set(hObject, 'Visible', 'on');
%set(hObject, 'Position', [331 12 240 64]); % For carbostar experiments
set(hObject, 'Position', [75 12 240 64]); % For tetrode experiments

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in fitAllCurve.
function fitAllCurve_Callback(hObject, eventdata, handles)
% hObject    handle to fitAllCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Characterize all data
handles = fitSDCurveCallback(handles, 0);
guidata(hObject, handles);

function TDTComError
% A function to indicate that an error occurred durring communication with
% the TDT system
error('Connection with TDT system is not active.');


% --- Executes on button press in getLaserEq.
function getLaserEq_Callback(hObject, eventdata, handles)
% hObject    handle to getLaserEq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new laser control parameters from TDT
handles = setLaserParameters(handles);
handles = setControlVoltages(handles);
guidata(hObject, handles);


function laserParamA_Callback(hObject, eventdata, handles)
% hObject    handle to laserParamA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'String', num2str(handles.curLaserParameters(1)));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function laserParamA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to laserParamA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function laserParamB_Callback(hObject, eventdata, handles)
% hObject    handle to laserParamB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'String', num2str(handles.curLaserParameters(2)));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function laserParamB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to laserParamB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function laserParamC_Callback(hObject, eventdata, handles)
% hObject    handle to laserParamC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'String', num2str(handles.curLaserParameters(3)));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function laserParamC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to laserParamC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function laserParamD_Callback(hObject, eventdata, handles)
% hObject    handle to laserParamD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'String', num2str(handles.curLaserParameters(4)));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function laserParamD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to laserParamD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in optServeConnect.
function optServeConnect_Callback(hObject, eventdata, handles)
% hObject    handle to optServeConnect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable the button
set(hObject, 'Enable', 'off');

% Toggle the state of the connection
if (handles.optServeConnected)
    % If there is a connection, end it
    handles = endOptServeComm(handles);
else
    % If there is no connection, start one
    handles = startOptServeComm(handles);
end

% Re-enable the button
set(hObject, 'Enable', 'on');

guidata(hObject, handles);

% --- Executes on button press in optServeCheck.
function optServeCheck_Callback(hObject, eventdata, handles)
% hObject    handle to optServeCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = checkOptServeConnectivity(handles);

% Update handles
guidata(hObject, handles);



function yPlotExistingFilePath_Callback(hObject, eventdata, handles)
% hObject    handle to yPlotExistingFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yPlotExistingFilePath as text
%        str2double(get(hObject,'String')) returns contents of yPlotExistingFilePath as a double

% --- Executes during object creation, after setting all properties.
function yPlotExistingFilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yPlotExistingFilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in controllerYPlotModify.
function controllerYPlotModify_Callback(hObject, eventdata, handles)
% hObject    handle to controllerYPlotModify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if the controller file exists
fileString = get(handles.yPlotExistingFilePath, 'String');
if (length(fileString) > length(handles.yPlotFileExtension)) && strcmp(fileString((end - length(handles.yPlotFileExtension) + 1):end), handles.yPlotFileExtension)
    % If the inputted string has the correct extension
    if ~exist(fileString, 'file')
        % If the file does not exist, then exit
        warndlg(sprintf('The file %s does not exist!', fileString), 'File does not exist');
        return;
    end
end

% Set the y-axis range for every plot in the Explore controller in TDT
newY = get(handles.controllerYPlotLevel, 'String'); % Get the y axis level that will be set to every pile plot

% Copy the file and make a backup, in case something goes horribly wrong
copyfile(fileString, ['C:\TDT\OpenEx\MyProjects\Acute Control PCA\' handles.yPlotFile '_bak.xpc']);

% Extract the information from the yPlot Controller
fileCell = loadTextFile(fileString);

% Go through each line of the file, and modify those which need to be
% modified
sectionTitle = '%PilePlot';

paramFirstSection = 'NAME=Y-Axis Range;VALUE=';
pollFirstSection = 'NAME=Poll Period;VALUE=';
line = 1;
while line < length(fileCell)
    thisLine = fileCell{line};
    if firstPartOfLineMatch(thisLine, sectionTitle)
        % line = line + 79;
        % We have passed the title section of one of the Pile Plots. Seach
        % for the parameter section
        while line < length(fileCell)
            % Make sure we don't go past the end of the file...
            thisLine = fileCell{line};
            if firstPartOfLineMatch(thisLine, pollFirstSection);
                % Modify the value
                fileCell{line} = [pollFirstSection num2str(100) ';GROUPID=Polling;'];
                break; % Break out of this inner while loop
            end
            
            if firstPartOfLineMatch(thisLine, paramFirstSection);
                % Modify the value
                fileCell{line} = [paramFirstSection newY ';GROUPID=Scaling;'];
                break; % Break out of this inner while loop
            end
            
            line = line + 1;
        end
    end
    
    line = line + 1;
end

% Write the file to disk
writeTextFile(fileCell, fileString);

function doesMatch = firstPartOfLineMatch(line, matchPortion)
portionSize = length(matchPortion);
doesMatch = length(line) >= portionSize && strcmp(line(1:portionSize), matchPortion);


% --- Executes on button press in browseYPlot.
function browseYPlot_Callback(hObject, eventdata, handles)
% hObject    handle to browseYPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
startPath = get(handles.yPlotExistingFilePath, 'String');

[fileName, pathName] = uigetfile({'*.xpc', 'Controller File (.xp)'}, 'Choose Controller File', startPath);

if ischar(fileName) && ischar(pathName)
    set(handles.yPlotExistingFilePath, 'String', [pathName fileName]); % Set the GUI editbox to be this new file name
    guidata(hObject, handles); % Update handles
end


function controllerYPlotLevel_Callback(hObject, eventdata, handles)
% hObject    handle to controllerYPlotLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of controllerYPlotLevel as text
%        str2double(get(hObject,'String')) returns contents of controllerYPlotLevel as a double
handles = checkNewInput(hObject, handles, 'yPlotVarData');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function controllerYPlotLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to controllerYPlotLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in modCSFTwoElectrode.
function modCSFTwoElectrode_Callback(hObject, eventdata, handles)
% hObject    handle to modCSFTwoElectrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of modCSFTwoElectrode

allH = [handles.modCSFTwoElectrode handles.modCSFTetrode];
set(allH, 'Value', 0);
set(hObject, 'Value', 1);
handles.modCSFIndex = 1;

set(handles.controllerExistingFilePath, 'String', handles.controllerFilePath{handles.modCSFIndex});
handles = checkControllerFileStatus(handles);

guidata(hObject, handles);


% --- Executes on button press in modCSFTetrode.
function modCSFTetrode_Callback(hObject, eventdata, handles)
% hObject    handle to modCSFTetrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of modCSFTetrode

allH = [handles.modCSFTwoElectrode handles.modCSFTetrode];
set(allH, 'Value', 0);
set(hObject, 'Value', 1);
handles.modCSFIndex = 2;

set(handles.controllerExistingFilePath, 'String', handles.controllerFilePath{handles.modCSFIndex});
handles = checkControllerFileStatus(handles);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tdtCharBlockName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tdtCharBlockName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = tsText(hObject, handles, neuronNum, tsNum)
newValRaw = str2double(get(hObject, 'String'));
curVal = handles.currentTriggerTimespan(neuronNum, tsNum);
newTVal = returnValidValue(newValRaw, [curVal 0 handles.maxDuration], false);

if ~isempty(handles.selectedStims{neuronNum})
    thisLoc = handles.selectedStims{neuronNum}(end);
    
    if thisLoc >= 1
        curTrigVal = handles.allNeuronsStruct{neuronNum}.triggerTimespan(thisLoc, tsNum);
        
        % Get the other timespan limit's value
        otherTsNum = mod(tsNum, 2) + 1;
        otherTsVal = handles.currentTriggerTimespan(neuronNum, otherTsNum);
        
        % Record the value for later
        if ~isnan(otherTsVal)
            % If there is already a selection in place for the other timespan
            % limit, make sure that the timespan is in order (lower on the
            % left, upper on the right)
            if (tsNum < otherTsNum && otherTsVal < newTVal) || (tsNum > otherTsNum && otherTsVal > newTVal)
                % Reverse the order
                handles.currentTriggerTimespan(neuronNum, [tsNum otherTsNum]) =  [otherTsVal newTVal];
            else
                handles.currentTriggerTimespan(neuronNum, tsNum) = newTVal;
            end
        else
            % If there is no selection already in place, then only assign a new
            % value if this is different than the currently set timespan value
            if (newTVal ~= curTrigVal)
                handles.currentTriggerTimespan(neuronNum, tsNum) = newTVal;
            else
                % If the value is the same as the current real value, clear the
                % "new" value
                handles.currentTriggerTimespan(neuronNum, tsNum) = NaN;
            end
        end
        
        % Update the timespan boxes
        updateTimespanBoxes(handles);
        
        % Update the voltage traces
        handles = drawVoltageTraces(handles, neuronNum);
    end
end


function tsTextL_2_Callback(hObject, eventdata, handles)
% hObject    handle to tsTextL_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tsTextL_2 as text
%        str2double(get(hObject,'String')) returns contents of tsTextL_2 as a double
handles = tsText(hObject, handles, 2, 1);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tsTextL_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tsTextL_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tsTextH_2_Callback(hObject, eventdata, handles)
% hObject    handle to tsTextH_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tsTextH_2 as text
%        str2double(get(hObject,'String')) returns contents of tsTextH_2 as a double
handles = tsText(hObject, handles, 2, 2);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tsTextH_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tsTextH_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tsTextL_1_Callback(hObject, eventdata, handles)
% hObject    handle to tsTextL_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tsTextL_1 as text
%        str2double(get(hObject,'String')) returns contents of tsTextL_1 as a double

handles = tsText(hObject, handles, 1, 1);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tsTextL_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tsTextL_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tsTextH_1_Callback(hObject, eventdata, handles)
% hObject    handle to tsTextH_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tsTextH_1 as text
%        str2double(get(hObject,'String')) returns contents of tsTextH_1 as a double

handles = tsText(hObject, handles, 1, 2);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tsTextH_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tsTextH_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
