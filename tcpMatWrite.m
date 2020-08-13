function tcpMatWrite(t, Y)
% This function will write a matrix to the server.  The protocol is as
% follows:
%
% N D1 D2 ... D_N A1 A2 ... A_M
% where N is the number of dimensions, DX is the size of dimension X of the
% matrix ,and AX is the value at the index X of the vectorization of Y.

% Find the size of each dimension of the matrix
D = size(Y);
if length(D) == 2 && D(2) == 1
    % If there are only two dimensions, and the size of the second is
    % 1...then there is really only 1 dimension
    D = D(1);
end

% Send the number of dimensions
fwrite(t, length(D));

% Send the size of each dimension
fwrite(t, D, 'double');
% fwrite_double(t, D);

% Finally, send the actual data itself
fwrite(t, Y(:), 'double');
% fwrite_double(t, Y(:));
end

function fwrite_double(t, y)
% Write the column array y as doubles to the stream t.

fwrite(t, typecast(y(:), 'uint8'));
end