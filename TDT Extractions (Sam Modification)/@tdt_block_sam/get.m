% get.m  Get tdt_block properties
%   Jason Ritt, BU, 5/12/11
%
%   get(DT)
%
% Displays some basic information about the tdt_block object DT.
% Very primitive.
%
% Must be at least ver 2.1 object

function get(DT,varargin)

if nargin<1
  disp('No object specified')
  % Can't happen
  return
end

if nargin>1
  warning('tdt_block: get method requests not yet implemented')
end

if nargout>0
  error('tdt_block: get method outputs not yet implemented')
end
  
display(DT)
