function write_waveforms_ntt(filename, timestamps, wave_mat)
% Write waveform samples to an ntt file as a tetrode
% wave_mat should be a matrix of the electrode data, with dimensions:
% numWaveforms x 4 x waveformLength

HEADERLENGTH = 2^14;

%dwScNumber = uint32(0);
%dwCellNumber = uint32(0);
%dnParams = uint32(zeros(8,1));

RECFILL = uint32(zeros(10,1));

%--- Find thresh crossings

% cross_inds = find(diff(vseries>=thresh)==1)+1;

% wave_mat = nan(length(cross_inds),WAVELEN);

timestamps = uint64(round(timestamps*1e6));

% for k=1:length(cross_inds)
%   wave_mat(k,:) = vseries(cross_inds(k)+WAVEOFF+(0:WAVELEN-1));
% end

ADBV = 1/2^15;
% ADBV = 1;
wave_mat = int16(round(wave_mat/ADBV));

%figure(1000)
%imagesc(wave_mat)

% --- Write file

FID = fopen(filename,'w','ieee-le');
% try
if FID<1
    error('Error opening file for writing')
end

fwrite(FID,uint8(ones(HEADERLENGTH,1)),'char'); % Fill up header with junk
fseek(FID,0,'bof');
fwrite(FID,'######## Neuralynx Data File Header','char*1'); % MClust looks for this
fseek(FID,HEADERLENGTH,'bof');
for k=1:size(wave_mat,1)
    fwrite(FID,timestamps(k),'uint64');
    fwrite(FID,RECFILL,'uint32');
    cur_waves = squeeze(wave_mat(k,:, :));
    fwrite(FID,cur_waves(:),'int16'); % The tetrode channels will be fed into the file, interlaced
end

FID = fclose(FID);
if FID
    warning('Error closing file')
end

% catch e
%     fclose(FID);
%     rethrow(e);
% end