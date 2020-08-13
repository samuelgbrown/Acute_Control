function [tOutAll, yOutAll] = solveIAFNeuronOutputResults(tF, gFun, maxSpikes, varargin)
% Integrate the IAF ode and output spike times
%
% This function will solve the ode equation which defines an Integrate And
% Fire neuron.  tF, or final time, is the number of milliseconds that the
% simulation should run for.  g(t) is a function which accepts the time in
% ms, and outputs the neuron's membrane conductance.  MaxSpikes represents
% the number of action potentials the ode waits for. If set to any non-zero
% positive integer, the ode will wait for that number of spikes to occur
% before stopping the simulation.  If set to 0, the ode will not stop until
% the time has expired.

% Define global parameters, to be used by the ode solver
global alpha beta E g vThresh noise lastT lastDiffT

%% Define user inputs
if nargin >= 4 && length(varargin{1}) == 4
    % If the user added parameter input of the correct format
    aNew = varargin{1}(1);
    bNew = varargin{1}(2);
    vNew = varargin{1}(3);
    ENew = varargin{1}(4);
else
    % If not, use the default values
    aNew = 1;
    bNew = 1.5;
    vNew = .5;
    ENew = 1;
end

%% Define parameters
% Static parameters\
alpha = aNew;
beta = bNew;
vThresh = vNew;
E = ENew;

% Input function (pulse train)
if isempty(gFun)
    % If no conductance input was specified
    amplitude = .7; % Amplitude of input
    tHigh = 10; % Number of ms each pulse to be on
    pulsePeriod = 100; % Period of each pulse
    g = @(t) amplitude*(mod(t, pulsePeriod) <= tHigh); % Define the conductance input function
else
    g = gFun; % Use the user's input function
end

%% Set up simulation
% Starting parameters
doPlot = false; % Should the function plot the ode output?
v0 = 0; % Starting voltage
t0 = 0; % Initial start-time (will change as action potentials occur, and the ode solver is restarted
if isempty(tF)
    % If the user does not supply a start time
    tF = 300; % Ending time
end

% Create empty output matrices
tOutAll = [];
yOutAll = [];

% Set up ode settings
refine = 4; % Increase the quality of integration
options = odeset('Events',@events, 'OutputSel',1, 'Refine', refine); % ,'OutputFcn', @odeplot);

%% Run simulation
numSpikes = 0; % Running total of the number of spikes that have occurred
while (true)
    % Set up timespan for this solving
    tSpan = [t0 tF];
    
    % Solve until next terminal event (a neuron spike), so that the solver can
    % be reset
    [tOut, yOut] = ode15s(@IAFNeuron, tSpan, v0, options);
    
    % Store outputs
    tOutAll = [tOutAll; tOut]; % Add the new times onto the tOutAll vectors
    yOutAll = [yOutAll; yOut]; % Add the new outputs onto the yOutAll vectors
    
    % Check ending conditions
    % Check if time has expired
    if tOutAll(end) >= tF
        break
    end
    
    % If time did not expire, but this line is reached, then a spike must
    % have occurred
    numSpikes = numSpikes + 1;
    
    % If there exists a maximum number of spikes, and the  maximum number
    % of spikes has occurred, break out of the loop
    if (numSpikes >= maxSpikes) && (maxSpikes ~= 0)
        break;
    end
    
    % Prepare the solver for the next round (directly after the action
    % potential
    options = odeset(options,'InitialStep',tOut(end)-tOut(end - 1),...
        'MaxStep',tOut(end)-tOut(1));
    
    t0 = tOut(end); % Use the last timestep as the first timestep of the next round
end

%% Plot results
if doPlot
    % Plot the outputs of the integration, including action potentials
    plot(tOutAll, yOutAll, 'b', tOutAll, g(tOutAll), 'r--', tEventsAll, yEventsAll, 'cx');
    legend('Neuron output', 'Stimulation', 'Action Potentials');
    ylabel('Voltage');
    xlabel('Time');
    title('Integrate-and-fire neuron');
end

end

%% Define event function
function [value, isterminal, direction] = events(t, y)
global vThresh
value = y(1) - vThresh; % The value of y that will trigger the event
isterminal = true; % End integration (will be started again in the loop)
direction = 1; % Only fire an event when the voltage increases through vThresh (not decreases)
end