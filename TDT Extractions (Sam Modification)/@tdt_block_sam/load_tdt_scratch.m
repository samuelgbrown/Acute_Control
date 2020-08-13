% Direct access to TDT tank files

data_dir = '/Volumes/RITTLABPORT/Turning_20110508/DataTanks/Turning_20110508_DT1_050811/after_turn';

file_base = 'Turning_20110508_DT1_050811_after_turn';

% TSQ is (fixed format) headers, TEV is (variable size) actual data


EVTYPE = struct('UNKNOWN',hex2dec('0'), ...
                'STRON',hex2dec('101'), ...
                'STROFF',hex2dec('102'), ...
                'SCALAR',hex2dec('201'), ...
                'STREAM',hex2dec('8101'), ...
                'SNIP',hex2dec('8201'), ...
                'MARK',hex2dec('8801'));
              
DFORM = struct('FLOAT',0, ...
               'LONG',1, ...
               'SHORT',2, ...
               'BYTE',3, ...
               'DOUBLE',4, ...
               'QWORD',5);
             

%% Load TSQ file into header

HEADER_SIZE = 40; % bytes

FID_TSQ = fopen([data_dir filesep file_base '.tsq'],'r');
if FID_TSQ<1
  error('Failed to open data tank header file!')
end

TANK_HEAD = struct('size',[],'type',[],'code',[],'channel',[], ...
               'sortcode',[],'timestamp',[],'offset',[],'format',[],'sampfreq',[]);

fseek(FID_TSQ,0,'bof');
TANK_HEAD.size = fread(FID_TSQ,inf,'*uint32',HEADER_SIZE-4);

fseek(FID_TSQ,4,'bof');
TANK_HEAD.type = fread(FID_TSQ,inf,'*uint32',HEADER_SIZE-4);

fseek(FID_TSQ,4+4,'bof');
TANK_HEAD.code = reshape(fread(FID_TSQ,inf,'4*uint8=>char',HEADER_SIZE-4),4,[])';

fseek(FID_TSQ,4+4+4,'bof');
TANK_HEAD.channel = fread(FID_TSQ,inf,'*uint16',HEADER_SIZE-2);

fseek(FID_TSQ,4+4+4+2,'bof');
TANK_HEAD.sortcode = fread(FID_TSQ,inf,'*uint16',HEADER_SIZE-2);

fseek(FID_TSQ,4+4+4+2+2,'bof');
TANK_HEAD.timestamp = fread(FID_TSQ,inf,'double',HEADER_SIZE-8);

fseek(FID_TSQ,4+4+4+2+2+8,'bof');
TANK_HEAD.offset = fread(FID_TSQ,inf,'uint64',HEADER_SIZE-8);

fseek(FID_TSQ,4+4+4+2+2+8+8,'bof');
TANK_HEAD.format = fread(FID_TSQ,inf,'*uint32',HEADER_SIZE-4);

fseek(FID_TSQ,4+4+4+2+2+8+8+4,'bof');
TANK_HEAD.sampfreq = fread(FID_TSQ,inf,'float',HEADER_SIZE-4);

FID_TSQ = fclose(FID_TSQ);
if FID_TSQ
  warning('Problem closing data tank header file')
end


%% Parse header

% How many data streams do we have?

%tags = double(unique(TANK_HEAD.code,'rows')); % Not all printable

% Looks like first tag is always "unknown" and contains entire data file
% Second tag is "mark" and is 20 bytes (?)
% Last tag is also "mark" and 20 bytes
%
% Remaining are data

% Determine which tags are data

F_stream = find(TANK_HEAD.type==EVTYPE.STREAM);
tags_stream = unique(TANK_HEAD.code(F_stream,:),'rows');

num_stream = size(tags_stream,1);


% How many channels in each stream?

chans = cell(num_stream,1);
num_chans = nan(num_stream,1);
for k=1:num_stream
  f = find(strcmp(tags_stream(k,:),cellstr(TANK_HEAD.code)));
  chans{k} = unique(TANK_HEAD.channel(f));
  num_chans(k) = length(chans{k});
end

%% Get a stream

tag_num = 1;
chan_num = 5;

ft = find(strcmp(tags_stream(tag_num,:),cellstr(TANK_HEAD.code)));
fc = find(chan_num==TANK_HEAD.channel(ft));

blk_size = double(TANK_HEAD.size(ft(fc)));
blk_locs = TANK_HEAD.offset(ft(fc));
blk_offs = [0; cumsum(blk_size)] + 1;

FID_TEV = fopen([data_dir filesep file_base '.tev'],'r');
if FID_TEV<1
  error('Failed to open data tank data file!')
end

DATA = nan(sum(blk_size),1);
for k=1:length(blk_size)
  fseek(FID_TEV,blk_locs(k),'bof');
  DATA(blk_offs(k) + (0:(blk_size(k)-1))) = fread(FID_TEV,blk_size(k),'long');
end

fclose(FID_TEV)

%% --- Grab from TEV file  OBSOLETE

Raw1_inds = find(strcmp('Raw1',cellstr(char(H_code))));
Raw1_size = H_size(Raw1_inds);
Raw1_locs = H_offset(Raw1_inds);
Raw1_btimes = H_timestamp(Raw1_inds);

if size(unique(Raw1_size))~=[1 1]
  error('Changing stream size')
else
  Raw1_size = Raw1_size(1);
end


FID_TEV = fopen([data_dir filesep file_base '.tev'],'r');

DATA = nan(200*Raw1_size,1);

for k=1:200
  DATA((k-1)*Raw1_size+(1:Raw1_size)) = fread(FID_TEV,Raw1_size,'int16');
end


fclose(FID_TEV);



