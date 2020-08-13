% display.m  Get tdt_block properties
%   Jason Ritt, BU, 5/12/11
%
%   display(DT)
%
% Displays some basic information about the tdt_block object DT.
% Mostly here to be called by Matlab whenever ans is a
% tdt_block object
%
% Must be at least ver 2.1 object

function display(DT)
  
C = DT.contents;
num_data = length(C);

disp(' ')
disp('/---\')
disp(' ')
disp(['tdt_block, ver ' DT.ver ', associated to files:'])
disp(' ')
disp(['    ' DT.tsq_file])
disp(['    ' DT.tev_file])
disp(' ')
disp(['with ' num2str(num_data) ' stores:'])
disp(' ')


for k=1:num_data
  disp(['    ' C(k).tag])
  disp(['    ' get_evname(C(k).type)])
  disp(['    ' num2str(C(k).samprate) ' Hz'])
  disp(['    ' 'Data type: ' get_dname(C(k).format)])
  switch get_evname(C(k).type)
    case 'STRON'
      % do nothing more
    case 'MARK'
      % do nothing more
    otherwise
      % SCALAR SNIP STREAM
      disp(['    ' 'Channels: ' num2str(reshape(C(k).hash.channels,1,[]))])
      tb = min(C(k).hash.timerange(:,1));
      te = max(C(k).hash.timerange(:,2));
      disp(['    ' 'Outer time interval: ' sprintf('%f %f sec',[tb te])])
      disp(['    ' '  Duration: ' sprintf('%f sec',te-tb)])
  end
  
  disp(' ')
end

disp(' ')
disp(['Zero time (MARK) is ' num2str(DT.metadata.timezero)])

disp(' ')
disp('\---/')
disp(' ')


