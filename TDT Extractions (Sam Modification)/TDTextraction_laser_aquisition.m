%6/30/2015
%MP
% Checks to see if Extracted laser files already exist
% Added a Catch in case 5 channels are in the Lasr Store
%All Laser recordings from TDT should be in a folder "Laser_Raw" directly
%under the mouse directory (ie with folders Raw and Processed)
% example C:\Users\palmieri\Desktop\Data\5335\Laser_Raw\LA_20150528_END
% where \Data is the current folder
% 5335 is the variable mousenum input into the function
% Laser_Raw is a directory that must be added which should contain the
% folder from TDT

function TDTextraction_laser_aquisition(mousenum)

TDT_laser_dir = fullfile(pwd,num2str(mousenum),'Laser_Raw'); %Mouse # Folder

files = dir(TDT_laser_dir);
TDT_aqusition_dir = {files(3:end).name}; %Folders that hold TDT files

% Check if files have already been extracted
existcheck = cellfun(@(x) exist(fullfile(pwd,num2str(mousenum),'Processed Laser',[x '.mat']),'file'),TDT_aqusition_dir);
nextract = find(existcheck == 0,1);
if isempty(nextract)
    disp('Files already extracted')
    return
end
TDT_aqusition_dir = TDT_aqusition_dir(nextract);

%Find the tsq files
for j = 1:length(TDT_aqusition_dir)
    TDT_dir_data{j}=fullfile(TDT_laser_dir,TDT_aqusition_dir{j});
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
for z=1:length(TDT_dir_data)
    clear prep_header Laser1 Laser2 Laser3 Spikes Pen Site
    name = TDT_aqusition_dir{z};
    prep_header = tdt_block([TDT_dir_data{z} filesep TDT_file_data{z}]);
    Laser1=getdata(prep_header,'Lasr','channel',1);%Shutter
    Laser2=getdata(prep_header,'Lasr','channel',2);%strength
    Laser3=getdata(prep_header,'Lasr','channel',3);%diode
    Laser4=getdata(prep_header,'Lasr','channel',4);%light sensor
    try % if 5 channels exsist
        Laser5=getdata(prep_header,'Lasr','channel',5);
        strudata.(name).shutter = Laser1.vals;
        strudata.(name).strength = Laser2.vals;
        strudata.(name).duration = Laser3.vals;
        strudata.(name).diode = Laser4.vals;
        strudata.(name).lightsensor = Laser5.vals;
        disp('5 channels')
    catch
        strudata.(name).shutter = Laser1.vals;
        strudata.(name).strength = Laser2.vals;
        strudata.(name).diode = Laser3.vals;
        strudata.(name).lightsensor = Laser4.vals;
    end
    disp('Lasr')
    

    
    disp('saving files as .MAT')
    
    
    
    %Changing the directory and adding new folders for processed data and
    %saving mat files
    
    
    dircheck = isdir(fullfile(pwd,num2str(mousenum),'Processed Laser'));
    if dircheck==0
        mkdir(fullfile(pwd,num2str(mousenum),'Processed Laser'))
    end
    
    
    ProcLaserLoc = fullfile(pwd,num2str(mousenum),'Processed Laser',TDT_aqusition_dir{z});
    save([ProcLaserLoc '.mat'],'-struct','strudata',name,'-v7.3')
    
end




