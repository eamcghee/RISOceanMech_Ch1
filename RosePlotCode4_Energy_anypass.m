%% Reshapes the data into a matrix where each column represents 
% one hour [Compute the desired statistic for each hour (e.g., mean)]

% caviats for chapter: 
% DR01 2016 023, corrupt, created "empty" arrays
% DR02 2016 008, had to be padded, only 4 hours
% DR02 2016 043, corrupt, created "empty" arrays

%% data 1, 2, 3 are velocity data (processed from counts to nominal velocity in template_extract_and_match_e2_loop.m
clear all 
close all
clc

% === Start timer ===
startTime = datetime('now');
disp(['Start time: ', datestr(startTime)])
tic 

addpath('/Volumes/CATALOGDR04/DR01/array/')
addpath('/Volumes/CATALOGDR03/DR01_MAT/') % DR01
%addpath('/Volumes/CATALOGDR02/MAT/') % DR02
%addpath('/Volumes/CATALOGDR00/') % DR03
%addpath('/Volumes/CATALOGDR03/') % DR02

sta = 'DR03';
year = '2016';

% Initialize the array to hold all hourly values across days
E_hr_all = [];
E_hr_all_v = [];

for icount = 091:091 % Loop through each day. 2015: 001:090; 2016: 001:091
    jjday = sprintf('%3.3d', icount); % Format the day count with leading zeros
    jjday
    
    %DR01 2015
    %load(['/Volumes/CATALOGDR03/DR01_MAT/', sta, '_', year, '_', jjday, '_matfile.mat'], 'data1', 'data2', 'data3', 'month', 'day');

    %DR02 2015
    %load(['/Volumes/CATALOGDR02/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'data1', 'data2', 'data3', 'month', 'day');
 
    %DR03 2015
    %load(['/Volumes/CATALOGDR00/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'data1', 'data2', 'data3', 'month', 'day');

    %DR01 2016, 01-22: %DR01_2016_023_matfile.mat. File might be corrupt.
    %DR01 2016, 24-75: 
    %load(['/Volumes/CATALOGDR03/DR01_MAT/', sta, '_', year, '_', jjday, '_matfile.mat'], 'data1', 'data2', 'data3'); %D, 'month', 'day');

    %DR02 2016, 01-19:
    %load(['/Volumes/CATALOGDR00/DR02_2016_001TO019_MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'data1', 'data2', 'data3', 'month', 'day');
    %DR02 2016, 20-42: DR02_2016_043_matfile.mat. File might be corrupt; DR02_2016_043_matfile.mat. File might be corrupt.
    %DR02 2016, 44-90: 
    %load(['/Volumes/CATALOGDR03/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'data1', 'data2', 'data3'); %, 'month', 'day');
    
    %DR03 2016, 01-90:
    %load(['/Volumes/CATALOGDR00/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'data1', 'data2', 'data3'); %, 'month', 'day');

    %DR01,02,03, 2016 HHZ day 091 only
    load(['/Volumes/CATALOGDR04/DATA/DR01/array/',sta,'_',year,'_',jjday,'_matfile.mat'], 'data1', 'data2', 'data3'); %, 'month', 'day');
    %% bandpass

    %fs=200; 
    fs = round(1/0.004999999888241291);
    % HHZ sachdr1.delta = 0.004999999888241291 (200 Hz sampling rate)
    % LHZ sachdr1.delta = 1 (1 Hz sampling rate)

    %fmin=0.12; 
    %fmax=1;

    fmin=0.07; 
    fmax=0.12;

    %fmin=0.03; 
    %fmax=0.06;
    
    %swell
    %fmin=0.03; 
    %fmax=0.12;

    %DF (double frequency)
    %fmin=0.14; 
    %fmax=0.2;

    % butterworth filter
    [bs,as]=butter(2,([fmin,fmax]/fs)*2);

    sdata1=filtfilt(bs,as,data1); % Apply Butterworth bandpass filter to data2 (Z component)
    sdata2=filtfilt(bs,as,data2); % Apply Butterworth bandpass filter to data2 (N component)
    sdata3=filtfilt(bs,as,data3); % Apply Butterworth bandpass filter to data3  (E component)

    % Pad sdata1 if needed
    target_length = 17280001;
    if length(sdata1) < target_length
        sdata1 = [sdata1; zeros(target_length - length(sdata1), 1)];
    end
    
    % Pad sdata2 if needed
    if length(sdata2) < target_length
        sdata2 = [sdata2; zeros(target_length - length(sdata2), 1)];
    end
    
    % Pad sdata3 if needed
    if length(sdata3) < target_length
        sdata3 = [sdata3; zeros(target_length - length(sdata3), 1)];
    end
    
    % get energy
    detectseries_mp=sqrt(sdata1.^2+sdata2.^2+sdata3.^2); % Energy of swell-band acceleration data (ZNE, 3-comp)
    detectseries_mp_v=sqrt(sdata1.^2); % Energy of swell-band acceleration data (Z, vertical-comp)

    % convert to hourly
    % Ensure detectseries_s has exactly 17,280,000 samples
    if numel(detectseries_mp) < 17280000
        detectseries_mp(end+1:17280000, 1) = 0; % Add zeros to make up the difference
    end

    % Sampling rate
    %fs = 200;
    samples_per_hour = 3600 * fs; % Total number of samples per hour

    % Ensure detectseries_s length is divisible by the number of hours
    detectseries_mp = detectseries_mp(1:(24 * samples_per_hour)); %size 17280000x1 double
    detectseries_mp_v = detectseries_mp_v(1:(24 * samples_per_hour)); %size 17280000x1 double, but normally 17279999x1 double

    try
        detectseries_mp = detectseries_mp(1:(24 * samples_per_hour));
        detectseries_mp_v = detectseries_mp_v(1:(24 * samples_per_hour));
    catch
        % because sometimes the array is ...999 instead of ...1000
        if length(detectseries_mp_v) == (24 * samples_per_hour - 1)
            detectseries_mp_v(end+1) = 0;  % Append a zero
            detectseries_mp_v = detectseries_mp_v(1:(24 * samples_per_hour));
        else
            rethrow(lasterror);  % Re-throw the error if it's not the expected case
        end
    end

    % Reshape the data into a matrix where each column represents one hour
    detectseries_reshaped = reshape(detectseries_mp, samples_per_hour, 24);
    detectseries_reshaped_v = reshape(detectseries_mp_v, samples_per_hour, 24);

    % Compute the desired statistic for each hour (e.g., mean)
    E_hr = mean(detectseries_reshaped, 1); % 1x24 array for hourly means
    E_hr_v = mean(detectseries_reshaped_v, 1); % 1x24 array for hourly means

    % Append the hourly values for this day to E_hr_all
    E_hr_all = [E_hr_all; E_hr];
    E_hr_all_v = [E_hr_all_v; E_hr_v];
end

% Assume E_hr_all is a 74x24 matrix
E_hr_all_flat = reshape(E_hr_all', [], 1); % Transpose and reshape into a column vector
E_hr_all_flat_v = reshape(E_hr_all_v', [], 1); % Transpose and reshape into a column vector

hourly_data = E_hr_all_flat;
hourly_data_v = E_hr_all_flat_v;

% Reshape the data from 74x24 to a single column (continuous hourly data)
%hourly_data = reshape(A, [], 1);  % Transpose and reshape to concatenate by columns

% Create a datetime array starting from 2015-01-01T00:00:00Z with an hourly increment
%start_time = datetime(2016, 1, 1, 0, 0, 0, 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
start_time = datetime(2016, 3, 31, 0, 0, 0, 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
num_hours = numel(hourly_data);  % Total number of hours
time_vector = start_time + hours(0:num_hours-1)';  % Generate timestamps

% Combine the data and time vector into one array/table
A3= table(time_vector, hourly_data, hourly_data_v, 'VariableNames', {'Time_UTC', 'Swellbandhi_hourly_avg_v', 'Swellbandhi_vert_hourly_avg_v'}); 
%A3= table(time_vector, hourly_data, hourly_data_v, 'VariableNames', {'Time_UTC', 'Swellbandlo_hourly_avg_v', 'Swellbandlo_vert_hourly_avg_v'}); 
%A4= table(time_vector, hourly_data_v, 'VariableNames', {'Time_UTC', 'midband_vert_hourly_avg'});

% Save to CSV
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR01_001to090_2016_pt07topt12Hz_3Comp_Vert_v.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR01_001to090_2016_pt03topt06Hz_3Comp_Vert_v.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR01_091to091_2016_pt03topt06Hz_3Comp_Vert_v.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR01_091to091_2016_pt07topt12Hz_3Comp_Vert_v.csv'

%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR02_001to090_2016_pt07topt12Hz_3Comp_Vert_v.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR02_001to090_2016_pt03topt06Hz_3Comp_Vert_v.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR02_091to091_2016_pt03topt06Hz_3Comp_Vert_v.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR02_091to091_2016_pt07topt12Hz_3Comp_Vert_v.csv'

%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR03_001to090_2016_pt07topt12Hz_3Comp_Vert_v.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR03_001to090_2016_pt03topt06Hz_3Comp_Vert_v.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR03_091to091_2016_pt03topt06Hz_3Comp_Vert_v.csv'
filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR03_091to091_2016_pt07topt12Hz_3Comp_Vert_v.csv'

writetable(A3, filename2);

% === End timer ===
endTime = datetime('now');
disp(['End time: ', datestr(endTime)])

% === Display total runtime ===
elapsedTime = toc;  % Total seconds elapsed

% Convert to hours, minutes, seconds format
hours = floor(elapsedTime/3600);
minutes = floor(mod(elapsedTime,3600)/60);
seconds = mod(elapsedTime,60);

fprintf('Total runtime: %d hours, %d minutes, %.2f seconds\n', hours, minutes, seconds);
%% acceleration data to hourly
% adata 1, 2, 3 are accleration data (processed from counts to nominal velocity, 
% then adata1=diff([data1(1);data1])*fs;, 
% then  adata1=data1.*tukeywin(length(data1),.005);
% in template_extract_and_match_e2_loop.m

clear all 
close all
clc
fprintf('velocity cleared, running acceleration now');

addpath('/Volumes/CATALOGDR04/DR01/array/')
addpath('/Volumes/CATALOGDR03/DR01_MAT/') % DR01
%addpath('/Volumes/CATALOGDR02/MAT/') % DR02
%addpath('/Volumes/CATALOGDR00/') % DR03
%addpath('/Volumes/CATALOGDR03/') % DR02

sta = 'DR03';
year = '2016';

% Initialize the array to hold all hourly values across days
E_hr_all = [];
E_hr_all_v = [];

for icount = 091:091% Loop through each day
    jjday = sprintf('%3.3d', icount); % Format the day count with leading zeros
    jjday
    
    %DR01 2015
    %load(['/Volumes/CATALOGDR03/DR01_MAT/', sta, '_', year, '_', jjday, '_matfile.mat'], 'adata1', 'adata2', 'adata3', 'month', 'day');

    %DR02 2015
    %load(['/Volumes/CATALOGDR02/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'adata1', 'adata2', 'adata3', 'month', 'day');
 
    %DR03 2015
    %load(['/Volumes/CATALOGDR00/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'adata1', 'adata2', 'adata3', 'month', 'day');

    %DR01 2016, 01-22: %DR01_2016_023_matfile.mat. File might be corrupt.
    %DR01 2016, 24-75: 
    %load(['/Volumes/CATALOGDR03/DR01_MAT/', sta, '_', year, '_', jjday, '_matfile.mat'], 'adata1', 'adata2', 'adata3'); %, 'month', 'day');

    %DR02 2016, 01-19:
    %load(['/Volumes/CATALOGDR00/DR02_2016_001TO019_MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'adata1', 'adata2', 'adata3', 'month', 'day');
    %DR02 2016, 20-42: DR02_2016_043_matfile.mat. File might be corrupt; DR02_2016_043_matfile.mat. File might be corrupt.
    %DR02 2016, 44-90: 
    %load(['/Volumes/CATALOGDR03/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'adata1', 'adata2', 'adata3'); %, 'month', 'day');
    
    %DR03 2016, 01-90:
    %load(['/Volumes/CATALOGDR00/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'adata1', 'adata2', 'adata3'); %, 'month', 'day');

    %DR01,02,03, 2016 HHZ day 091 only
    load(['/Volumes/CATALOGDR04/DATA/DR01/array/',sta,'_',year,'_',jjday,'_matfile.mat'], 'adata1', 'adata2', 'adata3'); %
    %% bandpass

    %fs=200; 
    fs = round(1/0.004999999888241291);
    % HHZ sachdr1.delta = 0.004999999888241291 (200 Hz sampling rate)
    % LHZ sachdr1.delta = 1 (1 Hz sampling rate)

    %fmin=0.12; 
    %fmax=1;

    fmin=0.07; 
    fmax=0.12;

    %fmin=0.03; 
    %fmax=0.06;
    
    %swell
    %fmin=0.03; 
    %fmax=0.12;

    %DF (double frequency)
    %fmin=0.14; 
    %fmax=0.2;

    % butterworth filter
    [bs,as]=butter(2,([fmin,fmax]/fs)*2);

    sdata1=filtfilt(bs,as,adata1); % Apply Butterworth bandpass filter to adata2 (Z component)
    sdata2=filtfilt(bs,as,adata2); % Apply Butterworth bandpass filter to adata2 (N component)
    sdata3=filtfilt(bs,as,adata3); % Apply Butterworth bandpass filter to adata3  (E component)

    % Pad sdata1 if needed
    target_length = 17280001;
    if length(sdata1) < target_length
        sdata1 = [sdata1; zeros(target_length - length(sdata1), 1)];
    end
    
    % Pad sdata2 if needed
    if length(sdata2) < target_length
        sdata2 = [sdata2; zeros(target_length - length(sdata2), 1)];
    end
    
    % Pad sdata3 if needed
    if length(sdata3) < target_length
        sdata3 = [sdata3; zeros(target_length - length(sdata3), 1)];
    end
    
    % get energy
    detectseries_mp=sqrt(sdata1.^2+sdata2.^2+sdata3.^2); % Energy of swell-band acceleration data (ZNE, 3-comp)
    detectseries_mp_v=sqrt(sdata1.^2); % Energy of swell-band acceleration data (Z, vertical-comp)

    % convert to hourly
    % Ensure detectseries_s has exactly 17,280,000 samples
    if numel(detectseries_mp) < 17280000
        detectseries_mp(end+1:17280000, 1) = 0; % Add zeros to make up the difference
    end

    % Sampling rate
    %fs = 200;
    samples_per_hour = 3600 * fs; % Total number of samples per hour

    % Ensure detectseries_s length is divisible by the number of hours
    %detectseries_mp = detectseries_mp(1:(24 * samples_per_hour)); %size 17280000x1 double
    %detectseries_mp_v = detectseries_mp_v(1:(24 * samples_per_hour)); %size 17280000x1 double, but normally 17279999x1 double

    try
        detectseries_mp = detectseries_mp(1:(24 * samples_per_hour));
        detectseries_mp_v = detectseries_mp_v(1:(24 * samples_per_hour));
    catch
        % because sometimes the array is ...999 instead of ...1000
        if length(detectseries_mp_v) == (24 * samples_per_hour - 1)
            detectseries_mp_v(end+1) = 0;  % Append a zero
            detectseries_mp_v = detectseries_mp_v(1:(24 * samples_per_hour));
        else
            rethrow(lasterror);  % Re-throw the error if it's not the expected case
        end
    end

    % Reshape the data into a matrix where each column represents one hour
    detectseries_reshaped = reshape(detectseries_mp, samples_per_hour, 24);
    detectseries_reshaped_v = reshape(detectseries_mp_v, samples_per_hour, 24);

    % Compute the desired statistic for each hour (e.g., mean)
    E_hr = mean(detectseries_reshaped, 1); % 1x24 array for hourly means
    E_hr_v = mean(detectseries_reshaped_v, 1); % 1x24 array for hourly means

    % Append the hourly values for this day to E_hr_all
    E_hr_all = [E_hr_all; E_hr];
    E_hr_all_v = [E_hr_all_v; E_hr_v];
end

% Assume E_hr_all is a 74x24 matrix
E_hr_all_flat = reshape(E_hr_all', [], 1); % Transpose and reshape into a column vector
E_hr_all_flat_v = reshape(E_hr_all_v', [], 1); % Transpose and reshape into a column vector

hourly_data = E_hr_all_flat;
hourly_data_v = E_hr_all_flat_v;

% Reshape the data from 74x24 to a single column (continuous hourly data)
%hourly_data = reshape(A, [], 1);  % Transpose and reshape to concatenate by columns

% Create a datetime array starting from 2015-01-01T00:00:00Z with an hourly increment
%start_time = datetime(2016, 1, 1, 0, 0, 0, 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
start_time = datetime(2016, 3, 31, 0, 0, 0, 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
num_hours = numel(hourly_data);  % Total number of hours
time_vector = start_time + hours(0:num_hours-1)';  % Generate timestamps

% Combine the data and time vector into one array/table
A3= table(time_vector, hourly_data, hourly_data_v, 'VariableNames', {'Time_UTC', 'Swellbandhi_hourly_avg_a', 'Swellbandhi_vert_hourly_avg_a'}); 
%A3= table(time_vector, hourly_data, hourly_data_v, 'VariableNames', {'Time_UTC', 'Swellbandlo_hourly_avg_a', 'Swellbandhlo_vert_hourly_avg_a'}); 
%A4= table(time_vector, hourly_data_v, 'VariableNames', {'Time_UTC', 'midband_vert_hourly_avg'});

% Save to CSV
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR01_001to090_2016_pt07topt12Hz_3Comp_Vert_a.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR01_001to090_2016_pt03topt06Hz_3Comp_Vert_a.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR01_091to091_2016_pt03topt06Hz_3Comp_Vert_a.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR01_091to091_2016_pt07topt12Hz_3Comp_Vert_a.csv'

%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR02_001to090_2016_pt07topt12Hz_3Comp_Vert_a.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR02_001to090_2016_pt03topt06Hz_3Comp_Vert_a.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR02_091to091_2016_pt03topt06Hz_3Comp_Vert_a.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR02_091to091_2016_pt07topt12Hz_3Comp_Vert_a.csv'

%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR03_001to090_2016_pt07topt12Hz_3Comp_Vert_a.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR03_001to090_2016_pt03topt06Hz_3Comp_Vert_a.csv'
%filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR03_091to091_2016_pt03topt06Hz_3Comp_Vert_a.csv'
filename2 = '/Users/em/PROJECTS/PEAKTEM/EnergyHrVA/Energy_all_comp_DR03_091to091_2016_pt07topt12Hz_3Comp_Vert_a.csv'
writetable(A3, filename2);

%% load only one mat file to see all arrays in there or 
% create "empty" arrays for adata and data, so the corrupt files still
% contain arrays but they are all equal to 1e-15. 

% clear all 
% close all
% clc
% 
% %addpath('/Volumes/CATALOGDR03/DR01_MAT/') % DR01
% %addpath('/Volumes/CATALOGDR02/MAT/') % DR02
% addpath('/Volumes/CATALOGDR03/MAT/') % DR02
% %addpath('/Volumes/CATALOGDR00/')          % DR03    
% 
% sta = 'DR02';
% year = '2016';
% 
% % Initialize the arrays to hold all hourly values across days
% %E_hr_all = [];
% %E_hr_all_v = [];
% 
% load(['/Volumes/CATALOGDR03/MAT/newDR02_2016_16_matfile.mat']); 

% for icount = 043:043   % Example loop (adjust for multiple days if needed)
%     jjday = sprintf('%3.3d', icount); % Format day with leading zeros (e.g., '001')
%     %jjdayy = '043';
%     % Load the full .mat file with all variables
%     %load(['/Volumes/CATALOGDR03/DR01_MAT/', sta, '_', year, '_', jjday, '_matfile.mat']);
% 
%     % DR01 2016, 01-22: %DR01_2016_023_matfile.mat. File might be corrupt.
%     % DR01 2016, 24-75: 
%     %load(['/Volumes/CATALOGDR03/DR01_MAT/', sta, '_', year, '_', jjday, '_matfile.mat']);
% 
%     %DR02 2016, 01-19:
%     %load(['/Volumes/CATALOGDR00/DR02_2016_001TO019_MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'data1', 'data2', 'data3','adata1', 'adata2', 'adata3', 'month', 'day');
%     %DR02 2016, 20-42: DR02_2016_043_matfile.mat. File might be corrupt; DR02_2016_043_matfile.mat. 
%     %DR02 2016, 44-90: 
%     load(['/Volumes/CATALOGDR03/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'data1', 'data2', 'data3', 'adata1', 'adata2', 'adata3', 'month', 'day');
%     % data1 = ones(17280001, 1) * 1e-15;
%     % data2 = data1;
%     % data3 = data1;
%     % adata1 = data1;
%     % adata2 = data1;
%     % adata3 = data1;
% 
%     % just load and save specific arrays (to a different folder, smaller file)
%     %save(['/Volumes/CATALOGDR03/MAT/new/', sta, '_', year, '_', jjdayy, '_matfile.mat'], 'data1', 'data2', 'data3', 'adata1', 'adata2', 'adata3');
%     % Save new arrays into the same .mat file
% 
%     % replicate "empty" adata's as data1 empty and append, save
%     % data1 = adata1(:,:);
%     % data2 = adata2(:,:);
%     % data3 = adata3(:,:);
%     %save(['/Volumes/CATALOGDR03/DR01_MAT/', sta, '_', year, '_', jjday, '_matfile.mat'], 'data1', 'data2', 'data3', '-append');
%     %save(['/Volumes/CATALOGDR03/MAT/new/', sta, '_', year, '_', jjdayy, '_matfile.mat'], 'adata1', 'adata2', 'adata3', '-append');
%     %save(['/Volumes/CATALOGDR03/MAT/', sta, '_', year, '_', jjdayy, '_matfile.mat'], 'month', 'day', '-append');
% end


%% Notes

%meanhp = mean(detectseries_hp);
%meansb = mean(detectseries_s);

%Va = [sdata3_h, sdata2_h]; % Combine seismic data components (N & E) into a matrix (2xN), where:
                           % sdata3_h = North component of seismic data
                           % sdata2_h = East component of seismic data

%Va = Va' * Va; % Compute the covariance matrix (2x2) by multiplying the transpose of Va with itself.
               % This results in a 2x2 symmetric matrix representing the energy distribution.

%[Ua, Da] = eig(Va); % Compute eigenvalues and eigenvectors of Va.
                    % Ua (2x2) contains the eigenvectors (principal directions).
                    % Da (2x2 diagonal) contains the eigenvalues, which represent signal power along each eigenvector.
                    % The eigenvectors describe the dominant orientation of motion.

%aza(i+1) = atan2d(Ua(1,2), Ua(2,2)); % Compute azimuth angle from the second eigenvector.
                                     % atan2d(y, x) returns the angle in degrees from the positive x-axis.
                                     % Here, Ua(1,2) and Ua(2,2) define the orientation of the second eigenvector.
                                     % The resulting aza value represents the azimuth in degrees from east.
                                     % Example interpretation: 010 means 100° from North.
%%

    %DR01 2014, Nov 21-Dec 31 (jdays = 325-365):
    %load(['/Volumes/CATALOGDR03/DR01_MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'adata1', 'adata2', 'adata3', 'month', 'day');
    %DR02 2014, Nov 21-Dec 31 (jdays = 325-365):
    %load(['/Volumes/CATALOGDR02/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'adata1', 'adata2', 'adata3', 'month', 'day');
    %DR03 2014, Nov 21-Dec 31 (jdays = 327-365):
    %load(['/Volumes/CATALOGDR00/MAT/',sta,'_',year,'_',jjday,'_matfile.mat'], 'adata1', 'adata2', 'adata3', 'month', 'day');
