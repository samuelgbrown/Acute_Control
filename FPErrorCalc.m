function e = FPErrorCalc(x, y)
% This function calculates the squared error of the inputs (specifically
% written to be used in the neuron control experiment, but, I mean...I
% guess you could totally just use this anywhere else...it's not like this
% is a particularly special implementation of the squared error
% calculation, aside from the NaN handling, but even that is probably
% pretty easily abstracted to plenty of other possible use-cases.  So, I
% suppose the question really remains...why did I write this function?
% Well, the first reason is that I didn't really want to put a function at
% the end of the script that I'm writing to test all of this stuff, because
% Matlab does this annoying thing that, when you put a function at the end
% of a script, it'll work totally fine except for the fact that you can't
% use the "Run Section" button.  It's not like I even use the "Run Section"
% feature all that much, but...honestly, it kind of just pisses me off to
% not have access to a feature that I should be able to use, aside from one
% small stupid thing...and if the only difference between the way I was
% going to do it and how it actually needs to be implemented is writing a
% single file, well...I can't really rationalize not doing it given how
% cheap hard-drive space is.  And all that, given the fact that I just
% needed a bit of a weird specific implementation thing with the whole
% NaN's stuff, well...yeah, I guess that's why I wrote this.  Still not
% really a good reason why I named it "FPSquaredError", as in
% "Fokker-Planck Squared Error", because it's absolutely, definitely not
% specifically useful for Fokker-Planck data...in fact, I may need to just
% try and use a likelihood for the FP data instead of a squared error...in
% which case I would be doubly wrong! You know, about the name AND what
% kind of error function to use.  Either way, I suppose this has all just
% been a long-winded way of saying that, like...there's nothing super
% special about this function. Except the NaN handling.  You know.  Aside
% from that...well...yeah. Nothing much).
%
% Ok, I just added a multiplier thing to penalize NaN's being produced.
% There...I guess that's something interesting, right?

% x is the FP model
% y is the real data

useSQE = true;
tol = 1e-5; % (Added 1/4/19)  The amount by which a 0 or 1 will be changed if exactly equal to a boundary

% Check inputs for NaN values
% invalid = any(isinf(x) | isnan(x) | isinf(y) | isnan(y) | x == 1 | x == 0, 2);
invalid = any(isinf(x) | isnan(x) | isinf(y) | isnan(y), 2);
x(invalid, :) = [];
y(invalid, :) = [];

if useSQE
    % Add a penalty for NaN?
%     mBase = .5;
%     mult = mBase*sum(invalid); % Multiply the cost for each value that is NaN.
    
    % Sum together the squared residuals
    e = sum((x - y).^2); % + mult;
else
    % OLD
    % Check for 0 or 1 values in the FP model (super extreme cases)
    x(x == 0) = tol;
    x(x == 1) = 1 - tol;
    
    e = y.*log(x) + (1 - y).*log(1 - x);
    e = -sum(e); % Want to maximize the likelihood, so negate this so we minimize its negative
    
    % NEW
    %     e = (x.^y ).*((1 - x).^(1 - y));
    %     e = -sum(e); % Want to maximize the likelihood, so negate this so we minimize its negative
    %     e = -prod(e);
    
    % NEW TEST (This is mathematically equivalent to the "new" formulation,
    % and generates the same results. It SHOULD generate the same results
    % as the original but does not)
    %     e = exp(e);
    %     e = -prod(e);
    
    %     e = e(~isnan(e)); % Take only the portions of the sum that are not NaN
    
end

if any(invalid)
    fprint('FPError-Foobar\n');
end

% Ensure that e is not empty
if isempty(e)
    e = Inf;
end

return;