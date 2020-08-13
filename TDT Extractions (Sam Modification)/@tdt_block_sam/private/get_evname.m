function evname_out = get_evname(evnum)

%% TDT encoding arrays

EVNUMS = [hex2dec('0');
          hex2dec('101');
          hex2dec('102');
          hex2dec('201');
          hex2dec('8101');
          hex2dec('8201');
          hex2dec('8801')];

EVNAMES = {'UNKNOWN';
           'STRON';
           'STROFF';
           'SCALAR';
           'STREAM';
           'SNIP';
           'MARK'};

% DNUMS = 0:5;
% 
% DNAMES = {'FLOAT';
%          'LONG';
%          'SHORT';
%          'BYTE';
%          'DOUBLE';
%          'QWORD'};
       
%% Reverse lookup

if length(evnum)==1

  f = find(evnum==EVNUMS);
  if ~isempty(f)
    evname_out = EVNAMES{f};
  else
    warning('Unknown Store type');
    evname_out = '';
  end
  
elseif length(evnum)>1

  evname_out = cell(length(evnum),1);
  for k=1:length(evname_out)
    f = find(evnum(k)==EVNUMS);
    if ~isempty(f)
      evname_out{k} = EVNAMES{f};
    else
      warning('Unknown Store type');
      evname_out{k} = '';
    end    
  end % for k

else
  
  warning('No input!!!')
  evname_out = '';

end % if length(evnum)


