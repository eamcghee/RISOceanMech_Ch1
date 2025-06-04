%% FETCH SAC FILES FROM IRIS WITH MATLAB

% Download sac files in format compatible with template_extract_e.m
% Save to: '/Users/em/PROJECTS/PEAKTEM/DATA/DR01_HHE_--_2015_065.sac'
% File naming convention: DR01_HHE_--_2015_065.sac
% Naming of dates: Ex.) 065 in name above is 3/16/2015
% Re-name files for use in next code of workflow.

% USER INPUTS:
% 1) Change dates & year in line of code 'for day_s = 28:30' & 'year = '2015-';'
% 2) Change to your folders:
%    folderPath = ['/Volumes/CATALOGDR03/DATA/'];
%    IrisFetch.SACfiles('XH','RS14','--','HHZ',data_starttime_day,data_endtime_day,'/Volumes/CATALOGDR00/DATA/RS14');

% NOTE: RUN SECTIONS SEPARATELY

%% RETRIEVE SAC FILES

clear
close all

addpath('/Volumes/CATALOGDR00/DATA/RS14')
addpath('/Users/em/PROJECTS/PEAKTEM/irisFetch-matlab-2.0.12');
folderPath = ['/Volumes/CATALOGDR03/DATA/'];
javaaddpath('IRIS-WS-2.0.19.jar'); %Java addpath for enabling IRIS

% Metadata for naming file
netcode ='XH'; % set manually in fetch line to reduce fetch time
stacode ='RS14'; % set manually in fetch line to reduce fetch time
%chancode ='HHZ'; % set manually in fetch line to reduce fetch time
%chancode ='HHN'; % set manually in fetch line to reduce fetch time
%chancode ='HHE'; % set manually in fetch line to reduce fetch time
dashcode ='--'; % set manually in fetch line to reduce fetch time

% manually set year, mon, day_s, then change for loop day_s range

% NEXT: %start DR03 2016 Nov download...check these and
% Jun-Dec of 2015, then rename and run

for day_s = 28:30 
    day_e = day_s + 1;
    day_s = num2str(day_s);
    day_e = num2str(day_e);
    disp(['Processing day start: ' num2str(day_s)]);
    disp(['Processing day end: ' num2str(day_e)]);
    year = '2015-'; % DRO1: 2016 315 was last day! for DR02, start 2015 07/31
    mon = '01-'; %10 next for DR02, 2015.
    time = '00:00:00';
    if numel(day_s) < 2
        % If day_s has less than two digits, add a '0' to the formatted string
        data_starttime_day = sprintf('%s%s%s%s', year, mon, '0', day_s, ' ', time);
        else
        % If day_s has two or more digits, proceed without adding '0'
        data_starttime_day = sprintf('%s%s%s%s', year, mon, day_s, ' ', time);
    end
     
    if numel(day_e) < 2
    % If day_e has less than two digits, add a '0' to the formatted string
    data_endtime_day = sprintf('%s%s%s%s', year, mon, '0', day_e, ' ', time);
    else
    % If day_s has two or more digits, proceed without adding '0'
    data_endtime_day = sprintf('%s%s%s%s', year, mon, day_e, ' ', time);
    end

    %data_starttime_day = sprintf('%s%s%s%s', year, mon, '0', day_s, ' ', time);
    data_endtime_day = sprintf('%s%s%s%s', year, mon, '0', day_e, ' ', time);
    disp(['Processing day start: ' num2str(data_starttime_day)]);
    
    %data_starttime_day = '2015-12-16 00:00:00'; %First day of sequence yyyy-mm-dd hr:mi:se
    %data_endtime_day = '2015-12-17 00:00:00'; %First day of sequence

    data_starttime_num = datenum(data_starttime_day); %First day datenum
    data_endtime_num = datenum(data_endtime_day); %First day datenum
    
    ndays = data_endtime_num - data_starttime_num;
    spd = 1;%segments per day
    ovlp = 1; %percent of overlap between days (to not miss events at the boundaries of the data segments
    
    datetimeObj = datetime(data_starttime_day, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    
    % Extract date only and convert back to string
    data_starttime_day_notime = datestr(datetimeObj, 'yyyy-mm-dd');
    
    % Extract year only
    data_starttime_year = datestr(datetimeObj, 'yyyy');
    
    % Convert the string to a datetime object
    datetimeObj = datetime(data_starttime_day_notime, 'InputFormat', 'yyyy-MM-dd');
    
    % Calculate the Julian date (number of days elapsed in the year (doy))
    jdate = datenum(datetimeObj) - datenum(datetime(datetimeObj.Year, 1, 1)) + 1;
    
    % Convert to a string with leading zeros (so 65 becomes 065, for ex.)
    jdatee = sprintf('%03d', jdate);
    
    % Display the result
    disp(['Original date without time: ', data_starttime_day_notime]);
    disp(['Julian date: ', num2str(jdate)]);
    
    data_starttime_num=datenum(data_starttime_day);%First day datenum
    nsegments=round((ndays+1)*spd*(1/(1-(ovlp/100))));
    
    javaaddpath('IRIS-WS-2.0.19.jar'); %Java addpath for enabling IRIS Web Services Calls
    
    %'/Users/em/PROJECTS/PEAKTEM/DATA/DR01_--_2015_065.sac'
    %wf_filename produces: DR01_--_2015_075.mat
    wf_filename=[stacode,'_',dashcode,'_',data_starttime_year,'_',jdatee,'.sac'];
    fullpath = fullfile(['/Volumes/CATALOGDR00/DATA/RS14/',wf_filename]);  % Get the full path of the file
    
    % reminders, double check...
    % netcode='XH';
    % stacode='DR01';
    % chancode='HHZ';
    % %chancode='HHN';
    % %chancode='HHE';
    % dashcode='--';
    
    %ex: irisFetch.SACfiles('IU','ANMO','*','*Z','2010-02-27 06:30:00','2010-02-27 07:30:00','/example/directory/')
    irisFetch.SACfiles('XH','RS14','--','HHZ',data_starttime_day,data_endtime_day,'/Volumes/CATALOGDR00/DATA/RS14');
    irisFetch.SACfiles('XH','RS14','--','HHE',data_starttime_day,data_endtime_day,'/Volumes/CATALOGDR00/DATA/RS14');
    irisFetch.SACfiles('XH','RS14','--','HHN',data_starttime_day,data_endtime_day,'/Volumes/CATALOGDR00/DATA/RS14');
end

%% RENAME HHE
% RUN FROM CATALOGDR00/DATA FOLDER

for jdatee = 28:30 %for any single-double-triple, change originalFileNameTemplate1 leading 0's
    jdatee = num2str(jdatee);
% FILE NAME CHANGE FOR HHE

% Original filename template
originalFileNameTemplate1 = 'XH.RS14..HHE.2015.%03d.00.00.00.SAC'; %ex. 00%03d for 1 dig day; 0%03d for 2 dig; %03d for 3 dig..

% Create the filename with the current jdatee value
originalFileName1 = strrep(originalFileNameTemplate1, '%03d', jdatee);

disp(['Original filename: ', originalFileName1]);

% Extracting relevant information using regular expressions
expression = 'XH\.(\w+)\.\.(\w+)\.(\d{4})\.(\d{3})\.(\d{2})\.(\d{2})\.(\d{2})\.SAC';
tokens = regexp(originalFileName1, expression, 'tokens');

% Rearrange extracted information
if ~isempty(tokens)
    network = tokens{1}{2};
    station = tokens{1}{1};
    year = tokens{1}{3};
    doy = tokens{1}{4};

    % Forming the new file name
    newFileName1 = sprintf('%s_%s_--_%s_%s.sac', station, network, year, doy);
    disp(['New filename: ',newFileName1]);

    % Rename the file
    %movefile(originalFileName1, newFileName1);

    %disp('File renamed successfully.');
%else

folderPath = '/Volumes/CATALOGDR00/DATA/RS14/';  % Specify the path to the folder
oldFilePath = fullfile(folderPath, originalFileName1);  % Construct the full path to the old file
newFilePath = fullfile(folderPath, newFileName1);  % Construct the full path to the new file
movefile(oldFilePath, newFilePath);  
    %disp('Failed to extract information from the filename.');
end
end
%% RENAME HHN

% FILE NAME CHANGE FOR HHN
for jdatee = 28:30 %for any single-double-triple, change orig...Template1 leading 0's
    jdatee = num2str(jdatee);

% Original filename template
originalFileNameTemplate2 = 'XH.RS14..HHN.2015.%03d.00.00.00.SAC';

% Create the filename with the current jdatee value
originalFileName2 = strrep(originalFileNameTemplate2, '%03d', jdatee);

disp(['Filename: ', originalFileName2]);

% Extracting relevant information using regular expressions
expression = 'XH\.(\w+)\.\.(\w+)\.(\d{4})\.(\d{3})\.(\d{2})\.(\d{2})\.(\d{2})\.SAC';
tokens = regexp(originalFileName2 ...
    , expression, 'tokens');

% Rearrange extracted information
if ~isempty(tokens)
    network = tokens{1}{2};
    station = tokens{1}{1};
    year = tokens{1}{3};
    doy = tokens{1}{4};

    % Forming the new file name
    newFileName2 = sprintf('%s_%s_--_%s_%s.sac', station, network, year, doy);
    disp(['New filename: ',newFileName2]);
    % Rename the file
    %movefile(originalFileName1, newFileName1);

    %disp('File renamed successfully.');
%else

folderPath = '/Volumes/CATALOGDR00/DATA/RS14/';  % Specify the path to the folder
oldFilePath = fullfile(folderPath, originalFileName2);  % Construct the full path to the old file
newFilePath = fullfile(folderPath, newFileName2);  % Construct the full path to the new file
movefile(oldFilePath, newFilePath);  

    %disp('Failed to extract information from the filename.');
end
end
%% RENAME HHZ

% FILE NAME CHANGE FOR HHZ
for jdatee = 28:30 %for any single-double-triple, change orig...Template1 leading 0's
    jdatee = num2str(jdatee);
% Original filename template
originalFileNameTemplate3 = 'XH.RS14..HHZ.2015.%03d.00.00.00.SAC';

% Create the filename with the current jdatee value
originalFileName3 = strrep(originalFileNameTemplate3, '%03d', jdatee);

disp(['Filename: ', originalFileName3]);

% Extracting relevant information using regular expressions
expression = 'XH\.(\w+)\.\.(\w+)\.(\d{4})\.(\d{3})\.(\d{2})\.(\d{2})\.(\d{2})\.SAC';
tokens = regexp(originalFileName3 ...
    , expression, 'tokens');

% Rearrange extracted information
if ~isempty(tokens)
    network = tokens{1}{2};
    station = tokens{1}{1};
    year = tokens{1}{3};
    doy = tokens{1}{4};

    % Forming the new file name
    newFileName3 = sprintf('%s_%s_--_%s_%s.sac', station, network, year, doy);
    disp(['New filename: ',newFileName3]);
    % Rename the file
    %movefile(originalFileName1, newFileName1);

    %disp('File renamed successfully.');
%else

folderPath = '/Volumes/CATALOGDR00/DATA/RS14/';  % Specify the path to the folder
oldFilePath3 = fullfile(folderPath, originalFileName3);  % Construct the full path to the old file
newFilePath3 = fullfile(folderPath, newFileName3);  % Construct the full path to the new file
movefile(oldFilePath3, newFilePath3);  

    %disp('Failed to extract information from the filename.');
end
end

%% NOTES ON METADATA FOR DR 01-03

% metadata from IRIS for DR01-03: 
% DR01 IRISDMC 2014-11-20	2016-11-10	Ross Ice Shelf - DR01	-77.767097	178.345703	2.0
% DR02 IRISDMC 2014-11-21	2016-12-11	Ross Ice Shelf - DR02	-77.824303	-178.424896	0.0
% DR03 IRISDMC 2014-11-22	2016-12-12	Ross Ice Shelf - DR03	-78.263	-175.115997	7.0
% Website: http://ds.iris.edu/mda/XH/?starttime=2014-01-01T00%3A00%3A00&endtime=2017-12-31T23%3A59%3A59
% You can click on each of the stations. Here's more info for DR01, Instruments/Channels: 
% Nanometrics Trillium 120 Sec PH Response/Quanterra
% BHE 20.0Hz IRISDMC | BHN 20.0Hz IRISDMC | BHZ 20.0Hz IRISDMC
% HHE 200.0Hz IRISDMC | HHN 200.0Hz IRISDMC | HHZ 200.0Hz IRISDMC
% LHE 1.0Hz IRISDMC | LHN 1.0Hz IRISDMC | LHZ 1.0Hz IRISDMC