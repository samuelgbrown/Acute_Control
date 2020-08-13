% This script will prepare the workspace to use the Fokker-Planck IAF
% model, more specifically the lookup table and associated interpolation
% software to find a firing probability for a neuron given its
% characterization and the parameters of the stimulation affecting it.

%% Parameters
% Irradiance conversion options
doIrradiance = true; % Should the FP model be converted to irradiance units?
maxAllowableIrradiance = 80;
laserEq = [-9.9757   75.5625  201.5847  -10.7840];

%% Prepare the workspace
FPTableFile = 'logisticFPTable.mat'; % The name of the file with the FP lookup table in it

load(FPTableFile);

if doIrradiance
    % Convert to irradiance values
    allStrsRaw = allStrs;
    betaRaw = beta;
    
    laserRadius = .1;
    powToIrr = 1/((laserRadius^2)*pi*1000);
    fIrrFromControlV = @(controlV, laserEq) powToIrr*polyval(laserEq, controlV);
    
    allStrs_RawMax = allStrsRaw(end);
    irrAtMaxG = fIrrFromControlV(allStrs_RawMax, laserEq); % The laser irradiance at 3V (the max of allStrs at time of coding) across the previous experiments.

    irrRatio = irrAtMaxG/allStrs_RawMax; % The ratio by which allStrs needs to be adjusted, and the inverse of that by which b0 must be adjusted
    
    % Perform the conversion
    allStrs = allStrsRaw*irrRatio; % Convert allStrs to be in units of irradiance
    beta = betaRaw/irrRatio; % Convert beta to be in units of irradiance
end

% Define the function that will perform the interpolation on the FP table
approxFun = @(GT, t) lininterpn(allStrs, allDurs, allAlphas, allSigmas, PAll, GT(:,1)*(t(2)/beta), GT(:,2), t(1), t(3));