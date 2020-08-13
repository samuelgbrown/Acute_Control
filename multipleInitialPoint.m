function [yMin, hMin] = multipleInitialPoint(f, x0, varargin)
% This function can be used to apply multiple starting conditions to an
% optimization function, and return the most optimal result.
%
% f: Function handle of the form: [y, h] = f(x)
%       , where y is the result of the optimization, h is the value to be
%       minimized by the optimization, and x is an initial condition.
%
% x0: A NxM matrix of initial conditions, where N is the number of
% conditions to test, and M is the size of the initial condition vector
% accepted by f.
%
% doParallel: a logical that dictates whether or not the function should
% use a parallel worker pool to optimize the calculation

%% Prepare the function
doParallel = false; % Should the optimizations be run in parallel?
numIterations = size(x0, 1); % Extract the number of iterations of the optimization

if nargin > 2
    doParallel = varargin{1};
end

if doParallel && isempty(gcp)
    parpool;
end

%% Loop through all iterations
hMinAll = zeros(numIterations, 1);
yMinAll = zeros(size(x0));
if doParallel && numIterations > 1
    %% Use parallel computation
    %     ticBytes(gcp);
    parfor iterNum = 1:numIterations
        [yMinAll(iterNum, :), hMinAll(iterNum)] = f(x0(iterNum, :)); % Evaluate the function
    end
    %     tocBytes(gcp);
else
    %% No not use parallel computation
    for iterNum = 1:numIterations
        [yMinAll(iterNum, :), hMinAll(iterNum)] = f(x0(iterNum, :)); % Evaluate the function
    end
end

%% Find the best result, and output it
[hMin, minInd] = min(hMinAll); % Get the index in hMinAll that has the lowest cost (best optimization)
yMin = yMinAll(minInd, :); % Return the index in yMinAll that represents the the lowest cost