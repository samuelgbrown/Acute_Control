%6/19/2015
%Paul
%Increase input arguement numbers

%6/18/2015
%Paul
%Changed script into function


% 6/1/2015
% Extract TDT files for acute prep
% MP


%Lasr Channel:
%    1:Shutter non-zero float value at the time it was open and close. 0s falls in between
%    2:Laser Strength
%    3:Laser Duration
%    4:Photodiode
%
%Raw1: Neural activities

%Alternative:
%Lasr Channel: Fast Store!
%    1:Shutter non-zero float value at the time it was open and close. 0s falls in between
%    2:Photodiode
%    3:Laser Strength
%LSDS Channel: Slow Store!
%    1.Laser Strength
%    2.Laser Duration
%    3.Empty
%
%Raw1 Channel: Neural activities
%Buts Channel:
%     1. Record'5341_5'
%     2. Penetration
%     3. Site

%Info:Penetration number, Site number
%Fq = 24414.0625;

%Input the mouse folder directory and the folders names that contain the
%TDT files
function rawStructHead = TDTextraction(mouse_num,session_num)

TDT_mouse_dir = [pwd '\' num2str(mouse_num)]; %Mouse # Folder
TDT_trial_dir = {session_num}; %Folders that hold TDT files

%%
%Find the tsq files
for j = 1:length(TDT_trial_dir)
    TDT_dir_data{j}=[TDT_mouse_dir '\raw\' TDT_trial_dir{j}];
    structlisting = dir(TDT_dir_data{j});
    i = 1;
    ftsq = structlisting(i).name;
    check = isempty(strfind(ftsq,'tsq'));
    while check == 1
        i = i + 1;
        ftsq = structlisting(i).name;
        check = isempty(strfind(ftsq,'tsq'));
        if i > length(structlisting)
            disp('tsq file not found')
            return
        end
    end
    TDT_file_data{j} = ftsq;
end

%Get header information from tsq file and extract the data using getdata
%function
%If stores are unfamiliar, look in prep_header to find all stores found in
%the TDT file. Then extract them using the format data{z} =
%getdata(prep_header{z}, 'Storename', 'channel', channel#)
%%
for z=1:length(TDT_dir_data)
    clear prep_header Laser1 Laser2 Laser3 Spikes Pen Site
    prep_header = tdt_block([TDT_dir_data{z} filesep TDT_file_data{z}]);
    Laser1=getdata(prep_header,'Lasr','channel',1);
    Laser2=getdata(prep_header,'Lasr','channel',2);
    Laser3=getdata(prep_header,'Lasr','channel',3);
    disp('Lasr')
    Spikes=getdata(prep_header,'Raw1','channel',1);
    disp('Raw1')
    Pen=getdata(prep_header,'Buts','channel',2);
    Site=getdata(prep_header,'Buts','channel',3);
    TrailState =getdata(prep_header,'Buts','channel',4);
    SpikeThresh =getdata(prep_header,'Buts','channel',5);
    disp('Buts')
    %     SequenceInfo = getdata(prep_header,'Valu','channel',1);
    disp('converting to structure')
    
    name = ['m' TDT_trial_dir{z} '_raw'];
    strudata.(name).shutter = Laser1.vals;
    strudata.(name).diode = Laser2.vals;
    strudata.(name).strength = Laser3.vals;
    strudata.(name).spikes = Spikes.vals;
    strudata.(name).time = Spikes.times;
    strudata.(name).Pen = Pen.vals;
    strudata.(name).Site = Site.vals;
    strudata.(name).TrailState = TrailState;
    strudata.(name).SpikeThresh = SpikeThresh;
    %     strudata.(name).SequenceEnd = SequenceInfo;
    disp('saving files as .MAT')
    
    rawStructHead = strudata.(name);
    
    %Changing the directory and adding new folders for processed data and
    %saving mat files
    
    oldcd = cd(TDT_mouse_dir);
    dircheck = isdir([TDT_mouse_dir '\Processed Data']);
    if dircheck==0
        mkdir('Processed Data')
    end
    
    
    
    cd([TDT_mouse_dir '\Processed Data'])
    mkdir(TDT_trial_dir{z})
    cd(TDT_trial_dir{z})
    save(name,'-struct','strudata',name,'-v7.3')
    disp('Done')
    cd(oldcd)
    
end



