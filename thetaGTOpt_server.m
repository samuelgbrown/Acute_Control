% This script can be run on a powerful "workhorse" computer to offload some
% of the work from the acquisition computer.  It is responsible for running
% the theta and GT optimizations that must be performed in real time for
% the Acute Adaptive Control experiment.
%
% The script will a) attempt to find the optimal theta that fits the real
% data acquired from an experiment to a lookup table of Fokker-Planck
% calculated neuron firing probabilities, and b) attempt to minimize a cost
% function across GT space using this calculated neuron model, to find
% stimulation parameters to exert optimal control over the neurons.

%% Prepare the script
% Set up parameters
version = 1; % The version number of the server script
useParallel = true; % Should parallel computation be used to speed up optimization (on problems where it will help)?
maxAllowableIrradiance = 80; % (mW/mm^2) Irradiance boundary of 80mW/mm^2
FPTableFile = 'logisticFPTable.mat'; % The name of the file with the FP lookup table in it

% Load the Fokker-Planck data
load(FPTableFile);
allStrsRaw = allStrs;
betaRaw = beta;

% Prepare the error calculation function
fErr = @FPErrorCalc;

% Set up optimization options
algorithm = 'sqp'; % 'interior-point;
thetaOpt  = optimoptions('fmincon', 'Algorithm', algorithm, 'Display', 'none'); % Optimization options;
GTOpt  = optimoptions('fmincon', 'Algorithm', algorithm, 'Display', 'none'); % Optimization options;

while true
    %% Begin communication with the client
    % Set up TCP/IP object
    t = tcpip('0.0.0.0', 30000, 'NetworkRole', 'server');
    t.OutputBufferSize = 2^18;
    t.InputBufferSize = 2^18;
    
    % Open communication (wait for client to connect)
    fprintf('\n\nWaiting for client connection...\n');
    fopen(t);
    
    %% Request some information from client
    % The irradiation ratio (ratio between the maximum irradiance measured at
    % the tip of the laser at highest power, and the highest strength that the
    % FP lookup table included) must be retreived from the client to complete
    % the calculations.  This allows all calculations to be done in units of
    % irradiance.
    fprintf('Requesting irradiation ratio...\n');
    irradiationRatio = tcpMatRead(t, true);
    
    % Convert values from arbitrary units to irradiance (mW/mm^2)
    allStrs = allStrsRaw*irradiationRatio;
    beta = betaRaw/irradiationRatio;
    
    fprintf('Irradiation ratio received!  Successfully connected to client!\n\n');
    
    % Set up boundaries on optimizations
    thetaLB = [0 0 .002]; % Lower bounds for t in the optimization
    % Note: Theta upper bound (thetaUB) will be set after data set is received,
    % because it depends on the max stimulation strength in the dataset
    
    GTLB = [0 0];
    GTUB = [max(allStrs) max(allDurs)];
    
    % Prepare standard starting points for optimization
    gPts = linspace(0, GTUB(1), 5); %[5 20 50];
    tPts = linspace(0, GTUB(2), 5); %[2 7 15];
%     gPts = gPts(2:3);
    tPts = tPts(2:end);
    [gX, tX] = ndgrid(gPts, tPts);
    gtX = [gX(:) tX(:)];
    
    % Prepare the approximation function
    approxFun = @(GT, t) lininterpn(allStrs, allDurs, allAlphas, allSigmas, PAll, GT(:,1)*(t(2)/beta), GT(:,2), t(1), t(3));
    
    % Prepare the cost function
    costFunLambda = 1e-5; %.01;
    costFun = @(P1, P2, GT) -(P1 * (1 - P2)) + costFunLambda*(GT(1).^2); % Cost function to optimize using control
    optCostFun = @(GT, t1, t2) costFun(approxFun(GT, t1), approxFun(GT, t2), GT);
    
    % Start a parallel pool again, if needed
    if useParallel
        if isempty(gcp)
            % Make a parallel pool that will not timeout for 12 hours
            parpool('IdleTimeout', 12*60);
        end
        
        % Send constant values to the parallel workers, so we don't need to
        % send the large arrays to them each time we start a computation
        allStrsPar = parallel.pool.Constant(allStrs);
        allDursPar = parallel.pool.Constant(allDurs);
        allAlphasPar = parallel.pool.Constant(allAlphas);
        allSigmasPar = parallel.pool.Constant(allSigmas);
        PAllPar = parallel.pool.Constant(PAll);
        
        % Define the approximation function in terms of the parallel pool
        % constants that we sent
        approxFunPar = @(GT, t) lininterpn(allStrsPar.Value, allDursPar.Value, allAlphasPar.Value, allSigmasPar.Value, PAllPar.Value, GT(:,1)*(t(2)/beta), GT(:,2), t(1), t(3));
    else
        % If not using the parallel pool, just use the same approximation
        % function as normal when we would otherwise use the special
        % parallel function
        approxFunPar = approxFun;
    end
    
    %% Start the loop to wait for new data from the acquisition computer
    fprintf('Waiting for data...\n');
    while strcmp(t.Status, 'open')
        %% Wait for new information to become available
        if t.BytesAvailable
            %% Read the data to be fit
            % The order of these reads is very important, and must be
            % synchronized between the server and client! (Could use a key/value
            % system, but a simple stream is probably more time efficient)
            %
            % N = number of data points
            % M = size of theta (number of neuron parameters)
            
            fprintf('Reading data...\n');
            
            newVal = tcpMatRead(t, true); % Read in the newest value
            if length(newVal) == 1 && newVal == -1
                % If a "-1" was received, assume that the connection is being terminated
                break;
            end
            
            GTH = newVal; % N x 4 matrix, columns are G, T, HA, HB
            lastThetas = tcpMatRead(t, true); % 2 x M matrix, rows are each neuron's previous thetas
            lastGTs = tcpMatRead(t, true); % 2 x 2 matrix, rows are the previous optimal GT's for each neuron
            
            % Clear the buffer from the client of all other data (just in
            % case)
            if t.BytesAvailable > 0
                fread(t, t.BytesAvailable);
            end
            
            % Extract the individual vectors
            GT = GTH(:, 1:2);
            HA = GTH(:, 3);
            HB = GTH(:, 4);
            
            lastThetaA = lastThetas(1, :);
            lastThetaB = lastThetas(2, :);
            
            lastGTA = lastGTs(1, :);
            lastGTB = lastGTs(2, :);
            
            % Assign failsafe values for the optimization targets
            bestThetaA = lastThetaA;
            bestThetaB = lastThetaB;
            bestGTA = lastGTA;
            bestGTB = lastGTB;
            
            % Start a parallel pool again, if needed
            if useParallel && isempty(gcp)
                % Make a parallel pool that will not timeout for 12 hours
                parpool('IdleTimeout', 12*60);
            end
            
            % Start a timer
            tic;
            
            %% Prepare final parameters
            % Calculate the upper bound for theta (must be done in loop because
            % it depends on G)
            thetaUB = [max(allAlphas) min(beta*(max(allStrs)/max(GT(:, 1))), maxAllowableIrradiance) max(allSigmas)]; % TODO: Beta restriction may be too severe
            
            % Prepare standard starting points for optimization (must be done
            % in loop because it depends on thetaUB)
            alphaPts = linspace(thetaLB(1), thetaUB(1), 4); % [.05 .2];
            betaPts = linspace(thetaLB(2), thetaUB(2), 4); % [.0005 .01];
            sigmaPts = linspace(thetaLB(3), thetaUB(3), 3); % [.0005 .01];
%             alphaPts = alphaPts(2:3);
%             betaPts = betaPts(2:3);
%             sigmaPts = .05; ;
            [aX, bX, sX] = ndgrid(alphaPts, betaPts, sigmaPts);
            thetaX = [aX(:) bX(:) sX(:)];
            
            %% Form the starting points to be used during this optimization run
            thisThetaXA = [thetaX;lastThetaA];
            thisThetaXB = [thetaX;lastThetaB];
            thisGTXA = [gtX;lastGTA];
            thisGTXB = [gtX;lastGTB];
            
            %% Calculate the optimal thetas for both neurons
            fprintf('Optimizing thetas...\n');
            
            try
                bestThetaA = multipleInitialPoint(@(x) fmincon(@(t) fErr(approxFunPar(GT, t), HA), x, [], [], [], [], thetaLB, thetaUB, [], thetaOpt), thisThetaXA, useParallel);
            catch e
                fprintf('Error (thetaA): %s\n\n', e.message);
            end
            
            try
                bestThetaB = multipleInitialPoint(@(x) fmincon(@(t) fErr(approxFunPar(GT, t), HB), x, [], [], [], [], thetaLB, thetaUB, [], thetaOpt), thisThetaXB, useParallel);
            catch e
                fprintf('Error (thetaB): %s\n\n', e.message);
            end
            
            t_optTheta = toc;
            
            %% Calculate the optimal GT's for both neurons
            fprintf('Optimizing GTs...\n');
            try
                bestGTA = multipleInitialPoint(@(x) fmincon(@(T) optCostFun(T, bestThetaA, bestThetaB), x, [], [], [], [], GTLB, GTUB, [], GTOpt), thisGTXA, false);
            catch e
                fprintf('Error (GTA): %s\n\n', e.message);
            end
            
            try
                bestGTB = multipleInitialPoint(@(x) fmincon(@(T) optCostFun(T, bestThetaB, bestThetaA), x, [], [], [], [], GTLB, GTUB, [], GTOpt), thisGTXB, false);
            catch e
                fprintf('Error (GTB): %s\n\n', e.message);
            end
            
            t_optGT = toc - t_optTheta;
            
            %% Assemble the data to be sent back to the client
            thisThetas = [bestThetaA;bestThetaB];
            thisGTs = [bestGTA;bestGTB];
            calcTime = [t_optTheta t_optGT];
            
            %% Send the data back to the client
            fprintf('Sending back data...\n');
            % Remember, the order is very important!
            tcpMatWrite(t, thisThetas);
            tcpMatWrite(t, thisGTs);
            tcpMatWrite(t, calcTime);
            
            fprintf('Done!\n\nWaiting for new data...\n');
        end
    end
    
    %% Disconnect from this instance
    fprintf('Disconnecting...\n');
    fclose(t);
end
