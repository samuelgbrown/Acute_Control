function dv = IAFNeuron(t, v)
% Integrate and fire neuron model
%
% This simple function represents the integrate and fire neuron model,
% shows how intracellular voltage evolves as a function of time, given the
% neuron's parameters (alpha, beta, and E), and some optogenetic-based
% conductance input.  It can be integrated using ode15s/ode45.

% Bring in parameters from calling function
global alpha beta E g
% Alpha is the decay rate, beta is the sensitivity to g, E is the
% characteristic potential of the ChR2 protein, and g is the conductance
% conferred by the protein (increases when exposed to light.
dv = -alpha*v + g(t)*beta*(E - v);