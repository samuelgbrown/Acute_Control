function [PFinal, POverTime, allP, t, V] = FokkerPlanckProbSim(G, tStimEnd, alpha, beta, vThresh, varargin)
% This function will simulate the probability distribution over time of
% aFokker
% leaky integrate and fire neuron subject to an optical stimulation.  In
% particular, it will measure the probability mass that has crossed the
% voltage threshold (causing a "spike"), integrate it, and return it as the
% probability that the given neuron would fire under the given conditions.
% The parameters are:edit prepare
%
% G - The strength of stimulation (unit TBD)
% T - The duration of stimulation (in ms)
% alpha - alpha in the leaky integrate and fire neuron equation
% beta - beta in the leaky integrate and fire neuron equation
% vThresh - the voltage threshold for the neuron
% doPlot (optional) - display plots that informative of the calculation?
% [numV numT] (optional) - the number of points in V and T there are in the
% calculation
sigma = .2; % Standard deviation of noise
doPlot = false;
numV = 301; % number of grid points in V
numT = 5001; % number of time steps to be iterated
tSimEnd = tStimEnd;
if nargin > 5
    sigma = varargin{1};
    if nargin > 6
        doPlot = varargin{2};
        if nargin > 7
            numV = varargin{3}(1);
            numT = varargin{3}(2);
            if nargin > 8
                tSimEnd = varargin{4};
            end
        end
    end
    
end

% close all;
% Specify unit parameters
% alpha = .5; % Alpha of unit
% beta = 1; % Beta of unit
% vThresh = 3; % Threshold for resetting
tauM = 1; % Membrane time constant
% tauM = 1.9; % This seems to make the values from this function match
% pretty well with the eulerMaruyama calculation
v0 = 0; % Initial voltage for all units
useGaussInit = false; % Use a gaussian initial condition (alternative: Use a burn in from P0: P(v0,1) = 1, otherwise 0)
useResetting = false; % Model resetting traces to 0 after hitting threshold

% Specify simulation parameters
VLower = -.67; % -1.5; TESTING % Smallest V
cutoffThresh = 1e-4; % Cutoff threshold for the probability being close to 1.  If the amount of mass left in the boundaries is less than cutoffThresh, it is assumed to be 0 (no traces are left below threshold, the neuron will always spike), and the simulation finishes

% Stimulation input parameters
tGOn = 0;  % Time at which stimulation will start
tGOff = tStimEnd; % Time at which stimulation will end

% Prepare simulation
% if useResetting
VUpper = vThresh;
% end
dV = (VUpper - VLower)/(numV - 1);
trueTF = tSimEnd/1; % True tFinal, in seconds
dt = trueTF/(numT - 1); % Change in t for each time step
V = (VLower):dV:VUpper;   % Vector of V values, to be used for plotting
t = 0:dt:trueTF; % Create a t-vector

P = zeros(numV,numT); % Iinitialize everything to zero

% Build input function
g = @(t) G*((tGOn < t) & (t < tGOff));
allG = g(t); % Precalculate the G vector

% Specify initial conditions
[~, v0Ind] = min(abs(V - v0)); % Find index of initial v0
if useGaussInit
    vInds = 1:numV; 
    P(:,1) = (1/(sigma*sqrt(2*pi)))*exp((-1/2)*(((V - v0)/sigma).^2));
else
%     P(v0Ind ,1) = 1; % All units start at V ~ 0 (Cannot be exactly on boundary, or it will get absorbed)
    P(:,1) = [3.31720085761631e-24;1.06948368856159e-23;2.18107909959184e-23;3.85675161769250e-23;6.38192668422996e-23;1.01839995162329e-22;1.59017115317698e-22;2.44878130295965e-22;3.73606814351045e-22;5.66277959803984e-22;8.54144704319526e-22;1.28346461732728e-21;1.92257076026823e-21;2.87221146000256e-21;4.28063978922048e-21;6.36559085237779e-21;9.44622731281109e-21;1.39894821631462e-20;2.06771536325887e-20;3.05028517543568e-20;4.49117993743101e-20;6.60020297634119e-20;9.68133678511835e-20;1.41741668663190e-19;2.07131421501825e-19;3.02121831449772e-19;4.39852389456028e-19;6.39177244410079e-19;9.27098245301428e-19;1.34221102374201e-18;1.93957615788042e-18;2.79759097025492e-18;4.02766424060933e-18;5.78780983334291e-18;8.30170840066277e-18;1.18853831280929e-17;1.69844560072677e-17;2.42260767994402e-17;3.44911578780491e-17;4.90146532376879e-17;6.95244894243405e-17;9.84336591038307e-17;1.39105285227689e-16;1.96217637832467e-16;2.76265780357532e-16;3.88249595298450e-16;5.44615528773778e-16;7.62543009227673e-16;1.06569813214467e-15;1.48661949190479e-15;2.06995744669490e-15;2.87686316889106e-15;3.99092351566820e-15;5.52616983236966e-15;7.63786341482830e-15;1.05369906924181e-14;1.45097043561774e-14;1.99433447809683e-14;2.73611977519244e-14;3.74688265006765e-14;5.12157056556582e-14;6.98770424854427e-14;9.51621741900412e-14;1.29357877692808e-13;1.75517445960471e-13;2.37709631611467e-13;3.21345736012765e-13;4.33608484593772e-13;5.84013196953908e-13;7.85140535585585e-13;1.05359130319498e-12;1.41122798456222e-12;1.88678577800699e-12;2.51795870698677e-12;3.35409524983731e-12;4.45967426482152e-12;5.91877765389443e-12;7.84083345701054e-12;1.03679746623015e-11;1.36844483314257e-11;1.80286207855486e-11;2.37082625697431e-11;3.11199677596626e-11;4.07737731955037e-11;5.33243032127375e-11;6.96100849318133e-11;9.07030707831136e-11;1.17970883717252e-10;1.53154884355898e-10;1.98467869342484e-10;2.56716071342178e-10;3.31451172909078e-10;4.27159303238013e-10;5.49495498885267e-10;7.05573923692866e-10;9.04326313647502e-10;1.15694370189213e-09;1.47741955961197e-09;1.88321614249958e-09;2.39608015396219e-09;3.04303893182439e-09;3.85761435543970e-09;4.88129869018368e-09;6.16534478751122e-09;7.77293260923237e-09;9.78178513013302e-09;1.22873194974533e-08;1.54064341113709e-08;1.92820492781374e-08;2.40885385351492e-08;3.00382099297364e-08;3.73890217400270e-08;4.64537456578910e-08;5.76108226196046e-08;7.13171925815386e-08;8.81234199047629e-08;1.08691480932117e-07;1.33815630092189e-07;1.64446815693952e-07;2.01721176686825e-07;2.46993217198376e-07;3.01874326694650e-07;3.68277390110696e-07;4.48468314129742e-07;5.45125382668309e-07;6.61407446107001e-07;8.01032044257626e-07;9.68364661622870e-07;1.16852041408343e-06;1.40747956736074e-06;1.69221838776664e-06;2.03085692284661e-06;2.43282540117244e-06;2.90905102398454e-06;3.47216699346427e-06;4.13674567956162e-06;4.91955786552611e-06;5.83986002718983e-06;6.91971158783351e-06;8.18432404396350e-06;9.66244377204404e-06;1.13867701963973e-05;1.33944108181406e-05;1.57273743681036e-05;1.84331030470997e-05;2.15650444487883e-05;2.51832633180336e-05;2.93550927759736e-05;3.41558240374572e-05;3.96694329534585e-05;4.59893409279876e-05;5.32192068847008e-05;6.14737459932398e-05;7.08795698115316e-05;8.15760413621547e-05;9.37161374450719e-05;0.000107467309204724;0.000123012330628627;0.000140550123272084;0.000160296544097352;0.000182485121906601;0.000207367726460524;0.000235215153035687;0.000266317603913934;0.000300985047149261;0.000339547431956876;0.000382354739252999;0.000429776845281019;0.000482203175929544;0.000540042129320299;0.000603720244557680;0.000673681095224413;0.000750383887313732;0.000834301742838951;0.000925919652382195;0.00102573208235587;0.00113424022576630;0.00125194888879393;0.00137936300953479;0.00151698380976927;0.00166530458561143;0.00182480614830972;0.00199595193226787;0.00217918279347418;0.00237491152789418;0.00258351714591229;0.00280533894550567;0.00304067043339310;0.00328975314980759;0.00355277045867050;0.00382984137067062;0.00412101447194017;0.00442626203554024;0.00474547439668865;0.00507845467545780;0.00542491393241917;0.00578446684330702;0.00615662797812086;0.00654080876810621;0.00693631524068728;0.00734234659763537;0.00775799470553047;0.00818224455992306;0.00861397577556879;0.00905196514475873;0.00949489029419900;0.00994133445823100;0.0103897923725797;0.0108386772784495;0.0112863290118576;0.0117310231378239;0.0121709810736649;0.0126043811304142;0.0130293703865796;0.0134440772943042;0.0138466249047895;0.0142351445878208;0.0146077901096383;0.0149627519244641;0.0152982715279071;0.0156126557154075;0.0159042905859947;0.0161716551310060;0.0164133342491358;0.0166280310332599;0.0168145781809018;0.0169719483889161;0.0170992636038429;0.0171958030123063;0.0172610096705942;0.0172944956889474;0.0172960459038550;0.0172656199905168;0.0172033529872861;0.0171095542240445;0.0169847046667376;0.0168294527104011;0.0166446084725970;0.0164311366579437;0.0161901480820571;0.0159228899594487;0.0156307350744956;0.0153151699672860;0.0149777822767704;0.0146202473920703;0.0142443145689052;0.0138517926718470;0.0134445357044630;0.0130244282884208;0.0125933712493363;0.0121532674616894;0.0117060080976312;0.0112534594151529;0.0107974502100803;0.0103397600439267;0.00988210834603008;0.00942614447387414;0.00897343880031456;0.00852547488087158;0.00808364273856661;0.00764923328823596;0.00722343390708663;0.00680732514369480;0.00640187854388996;0.00600795555919728;0.00562630749188051;0.00525757642026405;0.00490229703901167;0.00456089934146870;0.00423371206507158;0.00392096681620584;0.00362280278773197;0.00333927198065927;0.00307034484106631;0.00281591622425991;0.00257581160023324;0.00234979341761314;0.00213756754735356;0.00193878973230369;0.00175307197431966;0.00157998879666026;0.00141908332587306;0.00126987314410469;0.00113185586963240;0.00100451443029254;0.000887322001268538;0.000779746585296292;0.000681255219662876;0.000591317800341102;0.000509410519156413;0.000435018914975134;0.000367640544498475;0.000306787282320696;0.000251987263449865;0.000202786484493433;0.000158750082185623;0.000119463309895062;8.45322342220047e-05;5.35841748037993e-05;2.62679110288241e-05;0];
end

% Set boundaries
P(end,1) = 0; % Dirichlet

% Normalize initial condition
P(:,1) = P(:,1)/sum(P(:,1));


% Central difference
% %iterate difference equations
% for j=0:(numT - 1)
%     for i=2:numV-1
%         if j == 0
%                         P(:,j + 1) = (1/(stdNoise*sqrt(2*pi)))*exp((-1/2)*(((V - v0)/stdNoise).^2));
% %             P(v0Ind,j+1) = 1;
%         else
%             dPdV = (P(i+1, j) - P(i-1, j))/(2*dV); % First derivative wrt V
%             dPdV2 = (P(i+1, j) - 2*P(i, j) + P(i-1, j))/(dV^2); % Second derivative wrt V
%
%             P(i,j+1) = P(i, j) + (dt/tauM)*(-dPdV*(-alpha*V(i) + beta*G) + (dPdV2/2)*(stdNoise^2)); % New value of P (WHY DOES THE DP/DV TERM CREATE A POSITIVE VALUE???)
%             %             P(i,j+1) = P(i, j) + (dt/tauM)*((dPdV2/2)*(stdNoise^2)); % New value of P
% %             P(i,j+1) = P(i, j) + (dt/tauM)*(-dPdV*(-alpha*V(i) + beta*G)); % New value of P
%         end
%     end
%
% %     P([1 end], j+1) = 0; % Set boundary condition
%     P(1,j+1) = P(2,j+1);          % P(1,j+1) found from no-flux condition
%     P(numV,j+1) = P(numV-1,j+1);  % P(numV,j+1) found from no-flux condition
%
%     P(:, j+1) = P(:, j+1)/sum(P(:, j+1)); % Normalize probability so that each time-step's total probability sums to 1
% end

% Crank-Nicolson
% Prepare constants for the CN calculation
w = dt/(2*tauM*dV); % Scalar constant
x = (sigma^2)/dV; % Scalar constant
y = 2*w*x; % Scalar constant (Neumann boundary condition)

uVec = ones(numV, 1); % Vector of unit value, used to convert a scalar into a vector of size numVx1
if useResetting
    injectVec = zeros(numV,1);
    injectVec(v0Ind) = 1; % Vector of all 0's except a 1 at v0, to be multiplied by inject to represent the mass of traces that have been reset
end
totalMassCrossed = zeros(numT,1);
%iterate difference equations
for j=1:(numT - 1)
    %% Update user
    if mod(j, round((numT - 1)/10)) == 0
%         fprintf('Starting timestep %d of %d...\n', j, numT - 1);
    end
    
    %% Mass calculation
    % Get stimulation strength
%     thisT = tSimEnd*(j/(numT - 1));
    G = allG(j + 1);
    
    % Build two tri-diagonal matrices
    gammaLower = -alpha*V(1:(end - 1)) + beta*G;
    gammaUpper = -alpha*V(2:end) + beta*G;
    B = diag((1 - y)*uVec, 0) + diag(w*(gammaLower + x), -1) + diag(w*(x - gammaUpper), 1);
    A = diag((1 + y)*uVec, 0) + diag(-w*(gammaLower + x), -1) + diag(w*(gammaUpper - x), 1);
    
    % Modify tri-diagonals for the lower boundary condition P_1^j+1 * (1 +
    % y) - y * P_2^j+1 = -y * P_1^j + y * P_2^j
    A(1,[1 2]) = [(1 + y) -y];
    B(1,[1 2]) = [-y y];
    
    % Solve the implicit equation
    P(:,j + 1) = A\B*P(:,j);
    
    % Upper boundary condition and resetting
    totalMassCrossed(j) = P(numV,j+1); % Amount of mass that has crossed (or is about to cross) the voltage threshold boundary
    
    if useResetting
        % Reinject the mass that was lost at the absorbing boundary
        % (voltage threshold)
        P(:,j+1) = P(:,j+1) + injectVec*totalMassCrossed(j); % Add the injected mass back into v0
    end
    
    % Use absorbing (Dirichlet) boundary for upper boundary (threhold)
    P(numV, j+1) = 0; % Set boundary condition (0 value) at high V boundary
    
    % Test stopping condition (to save computation cycles)
    if ~useResetting && (sum(P(:, j + 1)) < cutoffThresh)
        % The computation is close enough to 1 to assume that the neuron
        % will always spike
        if j + 2 <= numT
            P(:,(j + 2):end) = 0;
        end
        break;
    end
    
    %     % Use reflecting boundary at lower boundary (assume there is VERY
    %     % little mass near it: if there isn't very little mass, then lower the
    %     % boundary)
    %     P(1,j+1) = P(2,j+1);          % P(1,j+1) found from no-flux condition at low V boundary
    %     P(numV,j+1) = P(numV-1,j+1); % P(numV,j+1) found from no-flux condition
    
    %     P(:, j+1) = P(:, j+1)/sum(P(:, j+1)); % Normalize probability so that each time-step's total probability sums to 1
end
% Pfinal = sum(totalMassCrossed);
POverTime = 1 - sum(P,1); % Get the fraction of the total amount of mass still left inside of the boundaries at each time point (initial amount of mass is 1, so no division necessary for the fraction)
PFinal = POverTime(end); % Get the probability at the last time point
allP = P;

% fprintf('Finished.\n\n');

if doPlot
    figure;
    mesh(t, V, P);
    hold on;
    % plot(V,P(:,1));
    % plot(V,P(:,11));
    % plot(V,P(:,101));
    % plot(V,P(:,1001));
    % plot(V,P(:,2001));
    xlabel('Time');
    ylabel('V');
    zlabel('P(V,t)');
    
    %calculate approximation to the integral of c from V=0 to V=1
    s = zeros(1,numT);
    for j=1:numT
        s(j) = sum(P(1:numV-1,j))*dV;
    end
    
    figure;
    plot(t,s);
    xlabel('t');
    ylabel('c_{total}');
    title('Total mass over time');
    axis([0 trueTF 0 1.1*max(s)]);
    
    figure;
    plot(t, POverTime);
    xlabel('t')
    ylabel('Probability');
    title('Probability over time');
end