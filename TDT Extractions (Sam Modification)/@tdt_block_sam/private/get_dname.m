function dname_out = get_dname(dnum)

%% TDT encoding arrays

% EVNUMS = [hex2dec('0');
%           hex2dec('101');
%           hex2dec('102');
%           hex2dec('201');
%           hex2dec('8101');
%           hex2dec('8201');
%           hex2dec('8801')];
% 
% EVNAMES = {'UNKNOWN';
%            'STRON';
%            'STROFF';
%            'SCALAR';
%            'STREAM';
%            'SNIP';
%            'MARK'};

DNUMS = 0:5;

DNAMES = {'FLOAT';
         'LONG';
         'SHORT';
         'BYTE';
         'DOUBLE';
         'QWORD'};
       
%% Reverse lookup

if length(dnum)==1

  f = find(dnum==DNUMS);
  if ~isempty(f)
    dname_out = DNAMES{f};
  else
    warning('Unknown data type');
    dname_out = '';
  end
  
elseif length(dnum)>1

  dname_out = cell(length(dnum),1);
  for k=1:length(dname_out)
    f = find(dnum(k)==DNUMS);
    if ~isempty(f)
      dname_out{k} = DNAMES{f};
    else
      warning('Unknown data type');
      dname_out{k} = '';
    end    
  end % for k

else
  
  warning('No input!!!')
  dname_out = '';

end % if length(dnum)


