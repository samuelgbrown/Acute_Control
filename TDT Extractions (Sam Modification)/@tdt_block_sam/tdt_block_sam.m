% tdt_block.m  Class for accessing TDT data blocks
%   Jason Ritt, BU, 5/12/11
%
% VER 2.0, 1/6/2012
% VER 2.1, 5/24/2012
% VER 2.2, 5/30/2017 (Samuel Brown, BU)
%
%   DT = tdt_block(TSQ_path)
%
% Returns a TDT data block (tdt_block) object pointing to the TSQ file given
% by the full path (as a string) TSQ_path.  It is assumed that the
% corresponding TEV file has the same base name and is in the same
% directory.
%
% Note this class accesses BLOCKS, not the full TANK as conceived by TDT.  Given the
% complete independence of data structure from block to block, there didn't seem to be much
% reason to implement at the higher level of abstraction.
%
% NOTES:
%
% --- ver 2.0 is incompatible with earlier version.  Version string added to
% top level structure (hard coded in mk_class subfunction).  As it changes
% the fields of the top level object, should also break any effort to load
% ulabelled ver 1.0 objects with this code.  This is intentional, since the
% approach is completely different.  Original object parsed entire TSQ file
% and kept it in memory.  The blocks got too large for this to work, so ver
% 2.0 checks the TSQ file, but then keeps a minimal reference set of
% information, and getdata is modified to first get location info from TSQ
% before going to TEV file.
%
% --- ver 2.1 subtracts the first MARK time from all other timestamps, so
% that they can have reasonable values.  In fact there were many changes
% to the object between 2.0 and 2.1, but this is the first that radically
% alters the behavior of loaded data, so using 2.0 code on 2.1 objects may
% have unpredictable results.


% TSQ record structure is
% Size        long     4   Size of data, with header, in longs.
% Type        long     4   Snip, stream, scalar, etc.
% Code        long     4   4-character event name (cast as long)
% Channel     u short  2
% Sort Code   u short  2   Sort code for snip data.
% Time Stamp  double   8   Event start time in seconds.
% Event Off   __int64  8   Offset in TEV file or
%  or Strobe  double  (8)    raw value.
% Data Format long     4   float, long, etc.
% Frequency   float    4   Sampling frequency


function DT = tdt_block_sam(varargin)
if nargin<1
  disp('GUI load of block: not implemented')
  %DT = []; % This will probably fail when Matlab tries to get class structure
  DT = mk_class('','',struct);
  return
elseif nargin~=1
  disp('Unknown input arguments')
  DT = mk_class('','',struct);
  return
end

if isa(varargin{1},'tdt_block_sam')
    DT = varargin{1};
    return
end

%% Get files

%TSQ_fullname = '/Users/jritt/Expts/RealTimeFeedback/RawData/TDT/3592_Partial/3592_20110810/3592_3592_20110810.tsq';

if ~ischar(varargin{1})
    error('tdt_block argument must be a string (file name)')
end

TSQ_fullname = varargin{1};
[TSQ_path, TSQ_basename, TSQ_ext] = fileparts(TSQ_fullname);

if ~strcmp('.tsq',TSQ_ext)
    error('Chosen file not a TSQ file')
end

TEV_fullname = fullfile(TSQ_path,[TSQ_basename '.tev']);


%% Define constants and helpers


EVTYPE = struct('UNKNOWN',hex2dec('0'), ...
    'STRON',hex2dec('101'), ...
    'STROFF',hex2dec('102'), ...
    'SCALAR',hex2dec('201'), ...
    'STREAM',hex2dec('8101'), ...
    'SNIP',hex2dec('8201'), ...
    'MARK',hex2dec('8801'));

% Note: TICK store type is STRON, but no documentation for it

DFORM = struct('FLOAT',0, ...
    'LONG',1, ...
    'SHORT',2, ...
    'BYTE',3, ...
    'DOUBLE',4, ...
    'QWORD',5);


%%
%%%
%%%- Check TSQ file, get top level information
%%%

% TSQ is (fixed format) headers
% TEV is (variable size) actual data, pointed to by TSQ

%  CONTENTS top level is num_tags struct, with fields
%    tag - 4 character code
%    type - numeric for STREAM, SNIP, etc
%    samprate - Hz
%    format - numeric DFORM data type
%    hash - structure with the pointers
%
%  hash is dependent on store type


%% Open TSQ file

HEADER_SIZE = 40; % bytes


FID_TSQ = fopen(TSQ_fullname,'r');
if FID_TSQ<1
    error('Failed to open data block header file!')
end

%try

% (1) Get all tag names, sort into uniques
% (2) Loop through tags, getting information
% (3) Set up CONTENTS structure
%
% Will be a little bit inefficient by loading all data for a particular
% querry then selecting on subindices; sacrificing executation time to
% save on memory and simplicity

%% Find all tags in TSQ file

fseek(FID_TSQ,4+4,'bof');
tags_all = reshape(fread(FID_TSQ,inf,'4*uint8=>char',HEADER_SIZE-4),4,[])';

tags_unique = unique(tags_all,'rows');
tags_unique = tags_unique(all(double(tags_unique')>=32)',:);
% Cleans out tags with non-printing (0000,1000,2000 show up)

num_tags = size(tags_unique,1);

for k=1:num_tags
    disp([num2str(k) ': ' tags_unique(k,:)])
end

% Find all TSQ record indices for each tag

disp('--> Index table into TSQ')

tsq_inds = cell(num_tags,1);
CSTR = cellstr(tags_all); % Stupid Matlab
for tag_num = 1:num_tags
    tsq_inds{tag_num} = find(strcmp(tags_unique(tag_num,:),CSTR));
end

clear tags_all CSTR

% tsq_inds{tag_num} replaces check on tag name


%% Find/check data types for each tag

disp('--> Data types')

fseek(FID_TSQ,4,'bof');
types_all = fread(FID_TSQ,inf,'*uint32',HEADER_SIZE-4);

C_type = nan(num_tags,1);
for tag_num = 1:num_tags
    cur_type = unique(types_all(tsq_inds{tag_num}));
    if length(cur_type(:))~=1
        warning(['Tag ' tags_unique(tag_num,:) ' has inconsistent tag type'])
        C_type(tag_num) = EVTYPE.UNKNOWN;
    else
        event_fieldnames = fieldnames(EVTYPE); % reverse hash
        for k=1:length(event_fieldnames)
            if cur_type==EVTYPE.(event_fieldnames{k})
                C_type(tag_num) = cur_type;
            end
        end
        if isnan(C_type(tag_num))
            warning(['Tag ' tags_unique(tag_num,:) ' has unknown tag type'])
            C_type(tag_num) = EVTYPE.UNKNOWN;
        end
    end
end
clear types_all cur_type event_fieldnames


%% Find/check sampling rate

disp('--> Sampling rate')

% Should be same for all channels within a tag

fseek(FID_TSQ,4+4+4+2+2+8+8+4,'bof');
rates_all = fread(FID_TSQ,inf,'float',HEADER_SIZE-4);

C_samprate = nan(num_tags,1);
for tag_num = 1:num_tags
    cur_rate = unique(rates_all(tsq_inds{tag_num}));
    if length(cur_rate(:))>1
        warning(['Tag ' tags_unique(tag_num,:) ' has inconsistent rate'])
        C_samprate(tag_num) = nan; % redundant
    else
        C_samprate(tag_num) = cur_rate;
    end
end
clear rates_all cur_rate

%% Find/check data format

disp('--> Data format')

fseek(FID_TSQ,4+4+4+2+2+8+8,'bof');
format_all = fread(FID_TSQ,inf,'uint32',HEADER_SIZE-4);

C_format = nan(num_tags,1);
for tag_num = 1:num_tags
    cur_format = unique(format_all(tsq_inds{tag_num}));
    if length(cur_format(:))>1
        warning(['Tag ' tags_unique(tag_num,:) ' has inconsistent data format'])
        C_format(tag_num) = nan; % redundant
    else
        C_format(tag_num) = cur_format;
    end
end
clear format_all cur_format


%% Get channels for each tag

disp('--> Channels')

fseek(FID_TSQ,4+4+4,'bof');
channel_all = fread(FID_TSQ,inf,'*uint16',HEADER_SIZE-2);

C_chan_nums = cell(num_tags,1);
num_chan = nan(num_tags,1);
for tag_num = 1:num_tags
    C_chan_nums{tag_num} = unique(channel_all(tsq_inds{tag_num}));
    num_chan(tag_num) = length(C_chan_nums{tag_num});
end

% do not clear channel_all !

%% Get first MARK to subtract from timestamps

fseek(FID_TSQ,40+4,'bof'); % Assumes MARK is 2nd entry
format_mark = fread(FID_TSQ,1,'uint32');
if format_mark ~= EVTYPE.MARK
    warning('Second TSQ rec is not of MARK type: timestamps not normalized');
    time_zero = 0;
else
    fread(FID_TSQ,8,'uint8');
    time_zero = fread(FID_TSQ,1,'double');
    if isnan(time_zero)
        warning('Something wrong with MARK time: timestamps not normalized');
        time_zero = 0;
    end
end

%% Assemble hashes

fseek(FID_TSQ,4+4+4+2+2,'bof');
timestamp_all = fread(FID_TSQ,inf,'double',HEADER_SIZE-8)-time_zero;

fseek(FID_TSQ,0,'bof');
size_all = fread(FID_TSQ,inf,'*uint32',HEADER_SIZE-4);

C_hash = cell(num_tags,1);

for tag_num = 1:num_tags
    disp(['Working on tag ' tags_unique(tag_num,:)])
    C_hash{tag_num} = struct;
    switch C_type(tag_num)
        
        %------ Streams
        
        case EVTYPE.STREAM
            C_hash{tag_num}.channels = C_chan_nums{tag_num};
            C_hash{tag_num}.timerange = nan(num_chan(tag_num),2);
            C_hash{tag_num}.tsqinds = cell(num_chan(tag_num),1);
            C_hash{tag_num}.size = nan(num_chan(tag_num),1);
            C_hash{tag_num}.fixedsizeQ = false(num_chan(tag_num),1);
            
            cur_inds = tsq_inds{tag_num}; % inds into _all for this tag
            
            for k=1:num_chan(tag_num)
                if ~isempty(cur_inds)
                    f = find(channel_all(cur_inds)==C_hash{tag_num}.channels(k));
                    
                    C_hash{tag_num}.timerange(k,:) = timestamp_all(cur_inds(f([1 end])));
                    C_hash{tag_num}.tsqinds{k} = cur_inds(f);
                    cur_sizes = size_all(cur_inds(f));
                    C_hash{tag_num}.fixedsizeQ(k) = length(unique(cur_sizes))==1;
                    if C_hash{tag_num}.fixedsizeQ(k)
                        C_hash{tag_num}.size(k) = double(unique(cur_sizes));
                    else
                        C_hash{tag_num}.size(k) = nan;
                    end
                else
                    warning('Empty indices')
                end
                
            end
            
            %------ Snips
            
        case EVTYPE.SNIP
            
            % NOTE: Snips are stored as individual channels, not as groups for
            % stereotrodes, tetrodes, etc.  So here channels are kept separate,
            % and getdata will by default load waveforms from only a single
            % channel.  We will add a "group" parameter, but this could get
            % complicated.  For example, there is nothing that forces all
            % waveforms in a tetrode to get captured with the same timestamp!
            
            C_hash{tag_num}.channels = C_chan_nums{tag_num};
            C_hash{tag_num}.timerange = nan(num_chan(tag_num),2);
            C_hash{tag_num}.tsqinds = cell(num_chan(tag_num),1);
            C_hash{tag_num}.size = nan(num_chan(tag_num),1);
            C_hash{tag_num}.fixedsizeQ = false(num_chan(tag_num),1);
            
            C_hash{tag_num}.numwaves = nan(num_chan(tag_num),1);
            
            cur_inds = tsq_inds{tag_num}; % inds into _all for this tag
            
            for k=1:num_chan(tag_num)
                if ~isempty(cur_inds)
                    f = find(channel_all(cur_inds)==C_hash{tag_num}.channels(k));
                    C_hash{tag_num}.timerange(k,:) = timestamp_all(cur_inds(f([1 end])));
                    C_hash{tag_num}.tsqinds{k} = cur_inds(f);
                    C_hash{tag_num}.numwaves(k) = length(f);
                    cur_sizes = size_all(cur_inds(f));
                    C_hash{tag_num}.fixedsizeQ(k) = length(unique(cur_sizes))==1;
                    if C_hash{tag_num}.fixedsizeQ(k)
                        C_hash{tag_num}.size(k) = double(unique(cur_sizes));
                    else
                        C_hash{tag_num}.size(k) = nan;
                    end
                else
                    warning('Empty indices')
                end
                
            end
            
            
            %------ Scalars ("slow store")
            
        case EVTYPE.SCALAR
            
            C_hash{tag_num}.channels = C_chan_nums{tag_num};
            C_hash{tag_num}.timerange = nan(num_chan(tag_num),2);
            C_hash{tag_num}.data = cell(num_chan(tag_num),1);
            
            cur_inds = tsq_inds{tag_num}; % inds into _all for this tag
            
            for k=1:num_chan(tag_num)
                if ~isempty(cur_inds)
                    f = find(channel_all(cur_inds) == C_hash{tag_num}.channels(k));
                    C_hash{tag_num}.timerange(k,:) = timestamp_all(cur_inds(f([1 end])));
                    C_hash{tag_num}.data{k}.times = timestamp_all(cur_inds(f));
                    C_hash{tag_num}.data{k}.vals = nan(length(f),1); % NEED OFFSETS READ AS DOUBLES!
                    
                    for j=1:length(f)
                        fseek(FID_TSQ,HEADER_SIZE*(cur_inds(f(j))-1)+ 4+4+4+2+2+8,'bof');
                        C_hash{tag_num}.data{k}.vals(j) = fread(FID_TSQ,1,'double');
                    end
                    
                else
                    warning(['Empty indices in scalar hash (could be a button that was never pushed)'])
                end
            end
            
            %------ Mark (unused; only at file begining)
            
        case EVTYPE.MARK
            % Should happen only at beginning, and handled before this loop
            
            %------ Stron (ticks, TTL, DFORM is double (time?))
            
        case EVTYPE.STRON
            
            if C_chan_nums{tag_num}>0
                warning('Multiple STRON channels')
            end
            
            cur_inds = tsq_inds{tag_num}; % inds into _all for this tag
            
            if ~isempty(cur_inds)
                C_hash{tag_num}.times = timestamp_all(cur_inds);
                C_hash{tag_num}.timerange = [nan nan]; % TODO
            else
                warning('Empty indices')
            end
            
            
        otherwise
            warning(['Unknown hash type:' C_type(tag_num)]);
    end % switch TAGTYPE
end % for tag_num

clear channel_all timestamp_all


% catch
%
%   warning('Something amiss with TSQ load')
%
% end

%% Initial parsing done, close file

FID_TSQ = fclose(FID_TSQ);
if FID_TSQ
    warning('Problem closing data block header file')
end



%% Make structure/class
%
%
%

% List of all channels:
% DT.contents.hash.channels

CONTENTS = struct('tag',cellstr(tags_unique),'type',num2cell(C_type), ...
    'format',num2cell(C_format),'samprate',num2cell(C_samprate),'hash',cell(num_tags,1));

METADATA = struct('timezero',time_zero);

for tag_num = 1:num_tags
    CONTENTS(tag_num).hash = C_hash{tag_num};
end

DT = mk_class(TSQ_fullname,TEV_fullname,CONTENTS,METADATA);



end  % tdt_block()

function DT = mk_class(TSQ_fullname,TEV_fullname,CONTENTS,METADATA)

% This function is here mostly to enforce consistent class definition
% across different creation conditional paths

VERSTRING = '2.1';

DT_struct = struct('ver',VERSTRING,'tsq_file',TSQ_fullname,'tev_file',TEV_fullname,'contents',CONTENTS,'metadata',METADATA);
DT = class(DT_struct,'tdt_block_sam');

end % mk_class()
