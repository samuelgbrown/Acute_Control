function v = lininterpn(varargin)
% Modification assumes that the matrix to be interpolated is 4-D and has
% dimensions of equal size!!!  Made to be used with Fokker-Planck table by
% Sam Brown
% 
% "Fast linear Interpolation", by Jeffrey Wu, Mathworks File Exchange, v
% 1.3 From
% https://www.mathworks.com/matlabcentral/fileexchange/28376-faster-linear-interpolation
% Bounding modifications by Sam (12/3/18)
%
% linear interpolation - input 1-dimensional arrays X1, X2, ... , Xn, a n-dimensional array V, and x1, x2, ..., xn query values
% assumes Xi's are in increasing order
%
% Differences from matlab built-in :
%       much, much faster
%       if coordinate is exactly on the spot, doesn't look at neighbors.  e.g. interpolate([blah, blah2], [0, NaN], blah) returns 0 instead of NaN
%       extends values off the ends instead of giving NaN
%

tol = 1e-6; % Hardcoded tolerance value (for testing, should be input later?)
bds = [0 1]; % Hardcoded boundary values (for testing, should be input later?)

if (mod(length(varargin), 2) ~= 1), error('Invalid number of inputs.  Should be odd'); end
n = (length(varargin) - 1) / 2;
isFirst = true;
for i = 1:1:n
    if (length(varargin{i}) ~= size(varargin{n+1}, i)), error('length(X%d) does not match size(V, %d)', i, i); end
    if ~(length(varargin{n+1+i}) == 1 || (isFirst || length(varargin{n+1+i}) == y))
        error('Query value x%d should be just one value or the same size as the other query values', i);
    end
    if isFirst
        y = length(varargin{n+1+i});
    end
    if length((varargin{n+1+i})) == 1
        varargin{n+1+i} = repmat(varargin{n+1+i}, y, 1);
    end
    isFirst = false;
end
pindices = zeros(y, n);
oindices = zeros(y, n);
slopes = ones(y, n);
for i = 1:1:n
    x = varargin{n+1+i};
    X = varargin{i};
    
    %     pindex = find((x >= X), 1, 'last');
    %     oindex = find((x <= X), 1, 'first');
    
    map1=bsxfun(@ge,x(:),X(:).');
    interpUnder = all(~map1, 2); % Are there any all 0 rows, implying that interpolation exceeds bounds?
    [~, pindexPre]=max(fliplr(map1),[],2);
    pindex = size(map1, 2) - pindexPre + 1;
    
    map2=bsxfun(@le,x(:),X(:).');
    interpOver = all(~map2, 2); % Are there any all 0 rows, implying that interpolation exceeds bounds?
    [~, oindex]=max(map2,[],2);
    
    pindex(interpUnder) = oindex(interpUnder);
    oindex(interpOver) = pindex(interpOver);
    
    interpProper = (pindex ~= oindex);
    
    %     if any(interpUnder)
    %         warning('interpolating before beginning in dimension %d', i);
    %         pindex = oindex;
    %     end
    %     if any(interpOver)
    %         warning('interpolating after end in dimension %d', i);
    %         oindex = pindex;
    %     end
    if any(interpProper)
        Xp = X(pindex(interpProper));
        Xo = X(oindex(interpProper));
        thisX = x(interpProper);
        slopes(interpProper, i) = (thisX(:) - Xp(:)) ./ (Xo(:) - Xp(:));
    end
    
    pindices(:, i) = pindex;
    oindices(:, i) = oindex;
end
V = varargin{n+1};
% vSize = size(V, 1); % Size of the dimensions of V (which should be all the same size)
v = zeros(y, 1);
for indexgetter = 1:1: 2^n
    multiplier = ones(y, 1);
    %     indices2 = cell(1, n);  % Each cell represents a single dimension containing y elements
    indices = zeros(y, n);
    for i = 1:1:n
        index = mod(indexgetter, 2);
        indexgetter = (indexgetter - index)/2;
        if index == 0
            %             indices2{i} = pindices(:, i);
            indices(:, i) = pindices(:, i);
            multiplier = multiplier .* (1 - slopes(:, i));
        else
            %             indices2{i} = oindices(:, i);
            indices(:, i) = oindices(:, i);
            multiplier = multiplier .* slopes(:, i);
        end
    end
    
    %     vInds = sub2ind_nocheck(size(V), indices2{:}); % Use a custom sub2ind script that doesn't check for errors because wE'RE C R A Z Y (and also want this to run fast no h8 plz)
    %     vInds = sub2ind(size(V), indices{:});
    
    indices = indices - [zeros(y, 1) ones(y, 3)];
    %     vInds = sum(bsxfun(@power, 45*ones(y, 1), 0:3).*indices, 2);
%     vInds = indices*[1; vSize; vSize^2; vSize^3];
    vInds = indices*[1; 45; 2025; 91125];
    
    v = v + V(vInds) .* multiplier;
end

% Enforce bounding, according to a tolerance
v(v > (bds(2) - tol)) = bds(2);
v(v < (bds(1) + tol)) = bds(1);

end

