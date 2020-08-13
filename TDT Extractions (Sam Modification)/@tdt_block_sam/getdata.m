% getdata.m  tdt_block method for getting raw data from block
%   Jason Ritt, BU 5/12/11
% 
% VER 2.0, 1/6/2012
% VER 2.1, 1/6/2017 (Samuel Brown, BU)
%
%   data = getdata(DT,TAG,param1,value1,...)
%
% DT is a tdt_block object.  TAG is a string containing the four letter name
% of the store.  The parameter value pairs follow usual Matlab convention and
% are optional.
%
% The output data depends on the data requested.  For streams, data will be a structure
% with two fields:
%
%   data.times  A vector of timestamps
%   data.vals   A vector of stream values
%
% Parameter value pairs (defaults in parens):
% 
%   channel (1)                 Channel within stream to get
%   interval ([-inf inf])       Interval of times over which to load
%   convert (false)             Convert to double if true, leave in native format if false
%
% NOTES: Only streams currently implemented
%        Not all params available


%--- TODO: implement interval, use offsets only for requested blocks
%--- TODO: boolean to decide when to make time vector
%--- TODO: load multiple channels
%--- TODO: Sophisticated handling of grouped snips (tetrodes)

%--- TODO: Should 10 bytes be removed at beg or end of stream recs?


function data = getdata(DT,req_tag,varargin)

load_struct = parse_inputs(varargin);

%channel = varargin{1};

ft = find(strcmp({DT.contents.tag},req_tag),1);
if isempty(ft)
  error('Store name not found')
end

CURCONT = DT.contents(ft);
time_zero = DT.metadata.timezero;

%% Constants

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

TSQ_HEADER_SIZE = 40; % bytes

%% Open TEV and TSQ files

FID_TEV = fopen(DT.tev_file,'r');
if FID_TEV<1
  error('Failed to open block data file!')
end

FID_TSQ = fopen(DT.tsq_file,'r');
if FID_TSQ<1
  error('Failed to open block header file!')
end

%try  % guard againt unclosed file

%% Load the appropriate store type

evnum = CURCONT.type;

switch evnum

  
  %------ Streams
  
  case EVTYPE.STREAM
    
    % Find blocks
        
    req_chan = find(CURCONT.hash.channels==load_struct.channel);
    if isempty(req_chan)
      error('Not a valid channel for this store')
    end

    samprate = CURCONT.samprate;
    cur_inds = CURCONT.hash.tsqinds{req_chan}; % 1-offset TSQ record index
    num_blocks = length(cur_inds);

    % Determine fread data type
    
    switch CURCONT.format
      case DFORM.FLOAT
        dtype = '*single';
        ptype = 'single';
        fudgeFactor = 1;
      case DFORM.LONG
        dtype = '*int32';
        ptype = 'int32';
        fudgeFactor = 1;
     case DFORM.SHORT
        dtype = '*int16';
        ptype = 'int16';
        fudgeFactor = 2;
      case DFORM.BYTE
        dtype = '*uint8';
        ptype = 'uint8';
        fudgeFactor = 1;
      case DFORM.DOUBLE
        dtype = 'double';
        ptype = 'double';
        fudgeFactor = 1;
      case DFORM.QWORD
        dtype = 'int64=>double';
        ptype = 'double';
        fudgeFactor = 1;
      otherwise
        error('Store has unknown data type')
    end

    % Get block sizes

    if CURCONT.hash.fixedsizeQ(req_chan)
      blk_size = fudgeFactor*(double(CURCONT.hash.size(req_chan))-10); % STREAMs lose 10
    else
      % Must now collect all sizes
      disp('Non-constant block size; getting sizes...')
      blk_size = nan(num_blocks,1);
      for k=1:num_blocks
        fseek(FID_TSQ,TSQ_HEADER_SIZE*(cur_inds(k)-1),'bof');
        blk_size(k) = fudgeFactor*(fread(FID_TSQ,1,'uint32')-10);
      end
      d_offs = [0; cumsum(blk_size(1:end-1))];
      disp('...done')
    end
    
    %SBT: Select the block numbers that include the time range you want.
    if isstruct(load_struct) && isfield(load_struct,'interval')
        startBlock = max(1,floor(samprate * load_struct.interval(1) / blk_size));
        endBlock = ceil(samprate * load_struct.interval(2) / blk_size);
        startT = load_struct.interval(1); % Start time
        endT = load_struct.interval(2); % End time
    else
        startBlock = 1;
        endBlock = num_blocks;
        startT = -Inf; % Start time
        endT = Inf; % End time
    end
       
    data = struct('vals',[],'times',[]);
    
    
    % Pre-allocate
    
    if CURCONT.hash.fixedsizeQ(req_chan)
      %SBT - Account for interval size specified
      data.vals = nan((endBlock - startBlock + 1)*blk_size,1,ptype);
      data.times = nan(size(data.vals));
    else
      data.vals = zeros(sum(blk_size(:)),1,ptype);
      data.times = double(data.vals);
    end

    
    % Loop: get TEV offset from TSQ, read TEV block
    allK = startBlock:endBlock;
    curDataInd = 1; % Current index in the data.vals or data.times vector
    for kInd=1:length(allK) %SBT
      k = allK(kInd);
      % Get data block loc in TEV file from TSQ file
      
      fseek(FID_TSQ,4+4+4+2+2 + TSQ_HEADER_SIZE*(cur_inds(k)-1),'bof');
      cur_time = fread(FID_TSQ,1,'double') - time_zero;
      cur_loc = fread(FID_TSQ,1,'uint64');
      
      if k < num_blocks
          fseek(FID_TSQ,4+4+4+2+2 + TSQ_HEADER_SIZE*(cur_inds(k + 1)-1),'bof');
          next_time = fread(FID_TSQ,1,'double') - time_zero;
      else
          next_time = cur_time + blk_size/samprate;
      end
      
      % Get data from TEV file
      
      fseek(FID_TEV,cur_loc,'bof');
      
      if CURCONT.hash.fixedsizeQ(req_chan)
          %disp('Fixed')
          % Get all data from this block
          newVals = fread(FID_TEV,blk_size,dtype);
          newTimes = cur_time + ((next_time - cur_time)/blk_size)*(0:(blk_size - 1)); % New method for interpolating times
%           newTimes = cur_time+(0:blk_size-1)/samprate; % Old method for extrapolating times
          
          % Extract only data within the time bounds
          newDataInds = (startT <= newTimes) & (newTimes <= endT); % Indices in allTimes and allVals to be used
          numNewData = sum(newDataInds); % Number of new data points
          yermom = curDataInd:(curDataInd + numNewData - 1); % Indices in data.vals and data.times that new data is going into (kept *ahem* legacy name)
          data.vals(yermom) = newVals(newDataInds); % Get only those new vals that lie within the requested interval
          data.times(yermom) = newTimes(newDataInds); % Get only those new times that lie within the requested interval
          
          % Update the current index in data.vals and data.times
          curDataInd = curDataInd + numNewData;
      else
          %disp('Not Fixed')
          yermom = d_offs(kInd) + (1:blk_size(kInd));
          data.vals(yermom) = fread(FID_TEV,blk_size(k),dtype);
          data.times(yermom) = cur_time+(0:blk_size(k)-1)/samprate;
      end
      
      
    end % for k
    
    % Remove all nans from data.vals and data.times (should be a number of
    % trailing nans if an interval was defined)
    data.vals(isnan(data.vals)) = [];
    data.times(isnan(data.times)) = [];
    
    if max(abs(1-(1./diff(data.times))/samprate))>0.01
      warning('Not all timestamp increments within 1% of samprate')
    end
    
    
  %------ Snips
  
  case EVTYPE.SNIP
    
    % Find recs
        
    req_chan = find(CURCONT.hash.channels==load_struct.channel);
    if isempty(req_chan)
      error('Not a valid channel for this store')
    end

    cur_inds = CURCONT.hash.tsqinds{req_chan}; % 1-offset TSQ record index
    num_blocks = length(cur_inds);

    % Determine fread data type
    
    switch CURCONT.format
      case DFORM.FLOAT
        dtype = '*single';
        ptype = 'single';
        fudgeFactor = 1;
      case DFORM.LONG
        dtype = '*int32';
        ptype = 'int32';
        fudgeFactor = 1;
     case DFORM.SHORT
        dtype = '*int16';
        ptype = 'int16';
        fudgeFactor = 2;
      case DFORM.BYTE
        dtype = '*uint8';
        ptype = 'uint8';
        fudgeFactor = 1;
      case DFORM.DOUBLE
        dtype = 'double';
        ptype = 'double';
        fudgeFactor = 1;
      case DFORM.QWORD
        dtype = 'int64=>double';
        ptype = 'double';
        fudgeFactor = 1;
      otherwise
        error('Store has unknown data type')
    end

    % Pre-allocate (assume fixed block size)
    % Time stamps are for first sample of waveform

    blk_size = double(CURCONT.hash.size(req_chan))-10;
    
    data = struct('vals',zeros(blk_size,num_blocks,ptype), ...
                  'times',nan(num_blocks,1));
    
    % Loop: get TEV offset from TSQ, read TEV block
    
    for k=1:num_blocks
      
      % Get data block loc in TEV file from TSQ file
      
      fseek(FID_TSQ,4+4+4+2+2 + TSQ_HEADER_SIZE*(cur_inds(k)-1),'bof');
      cur_time = fread(FID_TSQ,1,'double') - time_zero;
      cur_loc = fread(FID_TSQ,1,'uint64');
      
      % Get data from TEV file
      
      fseek(FID_TEV,cur_loc,'bof');
      
      yermom = (k-1)*blk_size + (1:blk_size);
      data.vals(:,k) = fread(FID_TEV,blk_size,dtype);
      data.times(k) = cur_time;
      
    end % for k
    
    
  %------ Marks (unused)  
    
  case EVTYPE.MARK
    
    disp('MARK getdata not implemented')
    data = struct('vals',[],'times',[]);
    
  %------ Strons
  
  case EVTYPE.STRON
    
    data = CURCONT.hash.times; % time_zero already subtracted in hash
    
  %------ Scalars
  
  case EVTYPE.SCALAR

    % Find blocks
    
    req_chan = find(CURCONT.hash.channels==load_struct.channel);
    if isempty(req_chan)
      error('Not a valid channel for this store')
    end

    data.vals = CURCONT.hash.data{req_chan}.vals;
    data.times = CURCONT.hash.data{req_chan}.times;
    
  otherwise
    
    warning('Unknown store type')
    data = struct('vals',[],'times',[]);

end % switch type

%catch
%  warning('Something went wrong with data extraction')
%end

%% Close TEV and TSQ files

FID = fclose(FID_TEV);
if FID
  warning('Something wrong with data file close')
end

FID = fclose(FID_TSQ);
if FID
  warning('Something wrong with header file close')
end

end % getdata


%%  Subfunctions

function load_struct = parse_inputs(args)

if isempty(args)

  % default params
  load_struct = struct('channel',1,'req_interval',[-inf inf],'convertQ',false);
  % return

else

  % parse input
  if ~isa(args,'cell')
    error('Problem with input arguments')
  end

  if mod(length(args),2)
    error('  data = getdata(DT,req_tag,param1,value1,...)')
  end

  for k=1:2:length(args)
    switch args{k}
      case 'channel'
        load_struct.channel = args{k+1};
      case 'interval'
        load_struct.interval = args{k+1};
      case 'convert'
        load_struct.convert = args{k+1};
      case 'group'
        disp('Grouping of tetrode channels not yet implemented')
      otherwise
        if ~isa(args{k},'char')
          warning('Parameter values should be strings!!!!')
        else
          warning(['Unknown parameter: ' args{k}])
        end
    end
  end % for k
end % else (isempty(args))

end % parse_inputs










