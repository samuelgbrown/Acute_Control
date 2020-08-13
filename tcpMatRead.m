function Y = tcpMatRead(t, varargin)
% This function will read a matrix from the client.  The protocol is as
% follows:
%
% N D1 D2 ... D_N A1 A2 ... A_M
% where N is the number of dimensions, DX is the size of dimension X of the
% matrix ,and AX is the value at the index X of the vectorization of Y.

if nargin > 1
   waitForData = varargin{1}; 
end

Y = [];

while waitForData && t.BytesAvailable == 0
    % Wait for new data
end
    
if t.BytesAvailable > 0
    % Get the number of dimensions, and their respective sizes
    N = fread(t, 1); % Warning: Can only handle matrices up to 255 dimensions.  Oh no.
    D = fread(t, N, 'double');
    %     D = fread_double(t, N);
    
    % Read the matrix data
    A = fread(t, prod(D), 'double');
    %     A = fread_double(t, prod(D));
    
    % Arrange the matrix data properly (if the value is not simply a
    % column vector or scalar)
    if N > 1
        Y = reshape(A, D');
    else
        Y = A;
    end
end

end

function y = fread_double(t, N)
% Read doubles from the stream t.
% N is the number of doubles to read

bpd = 8; % Bytes per double

y = typecast(uint8(fread(t, N*bpd)), 'double');
end