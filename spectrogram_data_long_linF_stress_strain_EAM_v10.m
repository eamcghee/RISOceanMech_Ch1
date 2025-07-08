clear all
close all
clc
%% merge a file with specific dates

% 001:090 = 01 Jan 2015 to 31 Mar 2015 (goes to 1 Apr 2015 00:00:00)
% 091:181 = 1 Apr 2015 to 30 Jun 2015 (goes to 1 Jul 2015 00:00:00)
% 182:273 = 1 Jul 2015 to 30 Sep 2015 (goes to 1 Oct 2015 00:00:00)
% 274:365 = 1 Oct 2015 to 31 Dec 2015 (goes to 1 Jan 2016 00:00:00)

% folderPath = '/Volumes/CATALOGDR00/DATA/Longfiles2/DR01';
% year = 2015;
% start_jday = 274;
% end_jday = 364;
% delta = 1.0;
% start_datetime = datetime(year,1,1) + days(start_jday - 1);
% 
% % Initialize
% all_data = [];
% 
% % Loop through each Julian day file
% for jday = start_jday:end_jday
%     jday_str = sprintf('%03d', jday);
%     fname = sprintf('DR01_LHZ_--_2015_%s.sac', jday_str);
%     fullfile_path = fullfile(folderPath, fname);
% 
%     if exist(fullfile_path, 'file')
%         try
%             [data, hdr] = rdsac(fullfile_path);  % or readsac if rdsac isn't available
%             all_data = [all_data; data];  % Append data
%         catch ME
%             warning('[ERROR READING FILE] %s — %s', fullfile_path, ME.message);
%             % Optional: pad with zeros if file is corrupt
%         end
%     else
%         fprintf('[MISSING FILE] %s\n', fullfile_path);
%         % Optional: pad with zeros
%     end
% end
% 
% % Write merged SAC
% outputFile = fullfile(folderPath, sprintf('DR01_LHZ_--_%d_%03dto%03d_merged.sac', year, start_jday, end_jday));
% mksac(outputFile, all_data, delta, start_datetime);
% 
% disp(['Merged SAC file written: ', outputFile]);

%% merge files with tukeywin and sample trim

clear all
close all
clc

folderPath = '/Volumes/CATALOGDR00/DATA/Longfiles2/DR03';
%year = 2014; start_jday = 325; end_jday = 365; %325:365 DR01 2014
%year = 2015; start_jday = 001; end_jday = 365; %001:365; %DR01 2015
year = 2016; start_jday = 001; end_jday = 314; %001:314; %DR01 2016
%start_jday = 325; end_jday = 365;
expected_n = 86400;  % Target sample count
delta = 1.0;         % Sample interval (1 Hz)
start_datetime = datetime(year,1,1) + days(start_jday - 1);

% Initialize
all_data = [];

% Loop through each Julian day
for jday = start_jday:end_jday
    jday_str = sprintf('%03d', jday);
    fname = sprintf('DR03_LHZ_--_%d_%s.sac', year, jday_str);
    fullfile_path = fullfile(folderPath, fname);

    if exist(fullfile_path, 'file')
        try
            [data, hdr] = rdsac(fullfile_path);  % or readsac
            n = length(data);

            % Pad or trim to 86400 samples
            if n > expected_n
                data = data(1:expected_n);
                fprintf('Trimmed %s from %d to %d samples.\n', fname, n, expected_n);
            elseif n < expected_n
                data = [data; zeros(expected_n - n, 1)];
                fprintf('Padded %s from %d to %d samples.\n', fname, n, expected_n);
            end

            % Apply Tukey window
            data = data .* tukeywin(expected_n, 0.005);

            % Append to cumulative array
            all_data = [all_data; data];

        catch ME
            warning('[ERROR READING FILE] %s — %s', fullfile_path, ME.message);
            all_data = [all_data; zeros(expected_n, 1)];  % Pad if unreadable
        end
    else
        fprintf('[MISSING FILE] %s\n', fullfile_path);
        all_data = [all_data; zeros(expected_n, 1)];  % Pad missing day
    end
end

% Write merged SAC
outputFile = fullfile(folderPath, sprintf('DR03_LHZ_--_%d_%03dto%03d_merged_tt.sac', year, start_jday, end_jday));
mksac(outputFile, all_data, delta, start_datetime);

disp(['✅ Merged SAC file written: ', outputFile]);


%% Manually Plot velocity to check data in merged sac file

% clear all
% close all
% clc
% 
% % === GLOBAL FONT SETTINGS ===
% set(groot, 'DefaultAxesFontName', 'Helvetica');
% set(groot, 'DefaultTextFontName', 'Helvetica');
% set(groot, 'DefaultAxesFontSize', 18);
% set(groot, 'DefaultAxesFontWeight', 'bold');
% set(groot, 'DefaultAxesLineWidth', 2);
% set(groot, 'DefaultAxesTickLength', [0.01 0.01]);
% set(groot, 'DefaultAxesTickDir', 'in');  % or 'out' if preferred
% 
% % _001to090 = 01 Jan 2015 to 31 Mar 2015 (goes to 1 Apr 2015 00:00:00)
% % _091to181 = 1 Apr 2015 to 30 Jun 2015 (goes to 1 Jul 2015 00:00:00)
% % _182to273 = 1 Jul 2015 to 30 Sep 2015 (goes to 1 Oct 2015 00:00:00)
% % _274to365 = 1 Oct 2015 to 31 Dec 2015 (goes to 1 Jan 2016 00:00:00)
% 
% % xlim([datetime(2015,1,1), datetime(2015,3,31,23,59,59)]);
% % xtick_dates = datetime(2015,1,1):days(1):datetime(2015,4,1);
% % xtick_dates = [xtick_dates, datetime(2015,4,1)];
% 
% % === Load the merged SAC file ===
% sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_001to365_merged_tt.sac'; %1 Jan 2015 to 31 Dec 2015
% %sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_001to090_merged.sac'; %1 Jan 2015 to 31 Mar 2015
% %sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_091to181_merged.sac'; %1 Apr 2015 to 30 Jun 2015
% %sacfile = '/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_182to273_merged.sac'; %1 Jul 2015 to 30 Sep 2015 
% %sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_274to365_merged.sac'; %1 Oct 2015 to 31 Dec 2015 
% %sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_274to364_merged.sac'; %1 Oct 2015 to 30 Dec 2015 
% 
% %read-in sac file & assign data and header (there are no headers here)
% [data, hdr] = rdsac(sacfile);
% 
% % Preprocess: remove trend, use sensitivity to convert to velocity
% data=detrend(data);
% sensitivity = 5.038830e8;
% data=data/sensitivity; % 'data' output is in velocity
% 
% % === Build time vector ===
% dt = 1; % sample interval in seconds
% n = length(data);
% t = (0:n-1) * dt; % time in seconds
% 
% % Define start datetime based on file (Day 091 = April 1, 2015)
% %start_datetime = datetime(2015, 1, 1) ;  % Jan 1
% %%start_datetime = datetime(2015, 1, 1)  + days(90);  % April 1 = day 91
% %%start_datetime = datetime(2015, 1, 1)  + days(181);  % Jul 1 = day 182
% %%start_datetime = datetime(2015, 1, 1)  + days(273);  % Oct 1 = day 274
% %time_vec = start_datetime + seconds(t); % Build datetime vector from seconds
% 
% % === Generate datetime vector for table ===
% start_datetime = datetime(2015, 1, 1, 0, 0, 0);  % Start at midnight Jan 1, 2015
% time_vec = start_datetime + seconds(0:n-1);     % One timestamp per sample
% 
% % === Combine into table ===
% T = table(time_vec(:), data(:), ...
%     'VariableNames', {'Time_UTC', 'Velocity_mps'});
% % writetable(T, 'DR01_2015_Velocity_with_Time.csv');
% 
% datalength = length(data);
% 
% % === Define the quarter to plot ===
% quarter_start = datetime(2015,1,1,0,0,0);
% quarter_end   = datetime(2015,3,31,23,59,59);
% file_base_name = 'Velocity_gram_DR01_2015_1Janto31Mar_M';
% 
% % === Index into the table using datetime logic ===
% quarter_idx = T.Time_UTC >= quarter_start & T.Time_UTC <= quarter_end;
% T_quarter = T(quarter_idx, :);
% 
% % === Plot ===
% H = figure(2);
% %plot(time_vec, data,'k-','linewidth',2);
% plot(T_quarter.Time_UTC, T_quarter.Velocity_mps, 'k-', 'LineWidth', 2);
% xlabel('Time (UTC)');
% ylabel('Amplitude (m/s)');
% title('DR01 LHZ Velocity');
% grid off;
% %ylim([-6e-4 6e-4]);
% %ylim([-6e-3 6e-3]);
% %ylim([-1e4 1e4]);
% %ylim([-1e6 1e6]);
% %ylim([-1e7 1e7]);
% 
% % TEST
% %xlim([datetime(2015,1,1), datetime(2015,3,31,23,59,59)]);
% %xtick_dates = datetime(2015,1,1):days(15):datetime(2015,3,31);
% %xtick_dates = [xtick_dates, datetime(2015,4,1)];
% 
% % JAN-MAR
% % xlim([datetime(2015,1,1), datetime(2015,3,31,23,59,59)]);
% % xtick_dates = datetime(2015,1,1):days(15):datetime(2015,3,31);
% % xtick_dates = [xtick_dates, datetime(2015,3,31)];
% % file_base_name = sprintf('Velocity_gram_DR01_2015_1Janto31Mar');
% 
% % APR-JUN
% % xlim([datetime(2015,4,1), datetime(2015,6,30,23,59,59)]);
% % xtick_dates = datetime(2015,4,1):days(15):datetime(2015,6,30);
% % xtick_dates = [xtick_dates, datetime(2015,7,1)];
% % file_base_name = sprintf('Velocity_gram_DR01_2015_1Aprto30Jun');
% 
% % JUL-SEP
% % xlim([datetime(2015,7,1), datetime(2015,9,30,23,59,59)]);
% % xtick_dates = datetime(2015,7,1):days(15):datetime(2015,9,30);
% % xtick_dates = [xtick_dates, datetime(2015,10,1)];
% % file_base_name = sprintf('Velocity_gram_DR01_2015_1Julto30Sep');
% 
% % OCT-DEC
% % xlim([datetime(2015,10,1), datetime(2015,12,31,23,59,59)]);
% % xtick_dates = datetime(2015,10,1):days(15):datetime(2015,12,31);
% % xtick_dates = [xtick_dates, datetime(2016,1,1)];
% % file_base_name = sprintf('Velocity_gram_DR01_2015_1Octto31Dec');
% 
% %full year merged
% % JAN-DEC
% % xlim([datetime(2015,1,1), datetime(2015,12,31,23,59,59)]);
% % xtick_dates = datetime(2015,1,1):days(30):datetime(2015,12,31);
% % xtick_dates = [xtick_dates, datetime(2016,1,1)];
% % file_base_name = sprintf('Velocity_gram_DR01_2015_1Janto31Dec');
% 
% % xticks(xtick_dates);
% % xticklabels(datestr(xtick_dates, 'dd-mmm'));
% % xtickangle(0);
% 
% % Set x-axis ticks global, modular
% xtick_dates = quarter_start:days(15):quarter_end;
% xtick_dates = [xtick_dates, quarter_end + seconds(1)];  % Ensure last tick
% xticks(xtick_dates);
% xticklabels(datestr(xtick_dates, 'dd-mmm'));
% xtickangle(0);
% 
% set(H, 'Units', 'inches');
% %set(H, 'Position', [1, 1, 30, 4]);
% set(H, 'Position', [1, 1, 20, 3]);
% 
% 
%     % === SAVE AFTER EACH FIGURE ===
%     folder_path = '/Volumes/CATALOGDR01/2025FIGS/'; 
%     file_ext = '.png';
%     file_ext = '.fig';
%     dpi = 600;
%     %fullFilePath = fullfile(folder_path, [file_base_name file_ext]);
%     %fullFilePath2 = fullfile(folder_path, [file_base_name2 file_ext]);
%     %set(H, 'Units', 'inches');
%     %set(H, 'Position', [1, 1, 30, 2]);  %[1, 1, 5, 2]
% 
%     % Save Fig 1 as PNG
%     pngFilePath = fullfile(folder_path, [file_base_name, '.png']);
%     %print(H, pngFilePath, '-dpng', ['-r' num2str(dpi)]);
% 
%     % Save Fig 1 as MATLAB FIG
%     figFilePath = fullfile(folder_path, [file_base_name, '.fig']);
%     %saveas(H, figFilePath); 

%% Plot velocity gram QUARTERS loop

clear all
close all
clc

% === GLOBAL FONT SETTINGS ===
set(groot, 'DefaultAxesFontName', 'Helvetica');
set(groot, 'DefaultTextFontName', 'Helvetica');
set(groot, 'DefaultAxesFontSize', 18);
set(groot, 'DefaultAxesFontWeight', 'bold');
set(groot, 'DefaultAxesLineWidth', 2);
set(groot, 'DefaultAxesTickLength', [0.01 0.01]);
set(groot, 'DefaultAxesTickDir', 'in');  % or 'out' if preferred

% === Set station & year ===
station = 'DR01';
year = '2014';

% === Load the merged SAC file ===
sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2014_325to365_merged_tt.sac'; %21 Nov 2014 to 31 Dec 2014
%sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_001to365_merged_tt.sac'; %1 Jan 2015 to 31 Dec 2015
%sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2016_001to314_merged_tt.sac'; %1 Jan 2016 to 10 Nov 2016

%read-in sac file & assign data and header (there are no headers here)
[data, hdr] = rdsac(sacfile);

% Preprocess: remove trend, use sensitivity to convert to velocity
data=detrend(data);
sensitivity = 5.038830e8;
data=data/sensitivity; % 'data' output is in velocity

% === Build time vector ===
dt = 1; % sample interval in seconds
n = length(data);
t = (0:n-1) * dt; % time in seconds

% === Generate datetime vector for table ===
start_datetime = datetime(2014, 11, 21, 0, 0, 0);  % Start Nov 21, 2014 00:00:00
%start_datetime = datetime(2015, 1, 1, 0, 0, 0);  % Start Jan 1, 2015  00:00:00
%start_datetime = datetime(2016, 1, 1, 0, 0, 0);  % Start Jan 1, 2016 00:00:00
time_vec = start_datetime + seconds(0:n-1);     % One timestamp per sample

% === Combine into table ===
T = table(time_vec(:), data(:), ...
    'VariableNames', {'Time_UTC', 'Velocity_mps'});
% writetable(T, 'DR01_2015_Velocity_with_Time.csv');

datalength = length(data);

% === Define quarters for 2014 ===
 quarters = {
     datetime(2014,11,21,0,0,0),  datetime(2014,12,31,23,59,59), '21Nov to 31Dec 2014';
 };

% === Define quarters for 2015 ===
% quarters = {
%     datetime(2015,1,1,0,0,0),  datetime(2015,3,31,23,59,59), '1Jan to 31Mar 2015';
%     datetime(2015,4,1,0,0,0),  datetime(2015,6,30,23,59,59), '1Apr to 30Jun 2015';
%     datetime(2015,7,1,0,0,0),  datetime(2015,9,30,23,59,59), '1Jul to 30Sep 2015';
%     datetime(2015,10,1,0,0,0), datetime(2015,12,31,23,59,59), '1Oct to 31Dec 2015';
% };

% === Define quarters for 2016 ===
% quarters = {
%     datetime(2016,1,1,0,0,0),  datetime(2016,3,31,23,59,59), '1Jan to 31Mar 2016';
%     datetime(2016,4,1,0,0,0),  datetime(2016,6,30,23,59,59), '1Apr to 30Jun 2016';
%     datetime(2016,7,1,0,0,0),  datetime(2016,9,30,23,59,59), '1Jul to 30Sep 2016';
%     datetime(2016,10,1,0,0,0), datetime(2016,11,10,23,59,59), '1Oct to 10Nov 2016';
% };

% === Loop through each quarter ===
for i = 1:size(quarters,1)
    % Extract start, end, and label
    quarter_start = quarters{i,1};
    quarter_end = quarters{i,2};
    label = quarters{i,3};
    
    % Index into T for current quarter
    quarter_idx = T.Time_UTC >= quarter_start & T.Time_UTC <= quarter_end;
    T_quarter = T(quarter_idx, :);
    
    % Plot
    H = figure(i); clf;
    plot(T_quarter.Time_UTC, T_quarter.Velocity_mps, 'k-', 'LineWidth', 2);
    xlabel('Time (UTC)');
    ylabel('Velocity (m/s)');
    title(sprintf('%s LHZ Velocity (%s)', station, label));
    grid off;
    ylim([-4e-3 4e-3]);

    % X-axis limits and ticks
    xlim([quarter_start, quarter_end + seconds(1)]);
    %xtick_dates = quarter_start:days(15):quarter_end;
    %xtick_dates = [xtick_dates, quarter_end + seconds(1)];
    %xticks(xtick_dates);
    %xticklabels(datestr(xtick_dates, 'dd-mmm'));
    xtickangle(0);

    % Set custom ticks and labels
    xtick_dates = [ ...
        datetime(2014,11,21), ...
        datetime(2014,11,27), ...
        datetime(2014,12,4), ...
        datetime(2014,12,11), ...
        datetime(2014,12,18), ...
        datetime(2014,12,25), ...
        datetime(2015,1,1)];
    xticks(xtick_dates);
    xticklabels(datestr(xtick_dates, 'mmm dd'));
    
    % Figure size and save name
    set(H, 'Units', 'inches');
    set(H, 'Position', [1, 1, 20, 3]);

    % Save
    %file_base_name = sprintf('Velocity_gram_DR01_2015_%s_M', label);
    %saveas(H, ['/Volumes/CATALOGDR01/2025FIGS', file_base_name, '.png']);
    %saveas(H, ['/Volumes/CATALOGDR01/2025FIGS', file_base_name, '.fig']);
    
    %Velocity_gram_DR01_2015_21Nov to 31Dec 2014_M

    save_folder = '/Volumes/CATALOGDR01/2025FIGS/';
    %file_base_name = sprintf('Velocity_gram_', station, '_', year, '_%s_M', label);
    file_base_name = sprintf('Velocity_gram_%s_%s_%s_M', station, year, label);

    % Construct full file names using fullfile (OS-safe)
    png_file = fullfile(save_folder, [file_base_name, '.png']);
    fig_file = fullfile(save_folder, [file_base_name, '.fig']);
    
    % Save the files
    saveas(H, png_file);
    saveas(H, fig_file);

    fprintf('✅ Saved %s.png and .fig\n', file_base_name);

end

%% Spectrogram in Pressure, Strain, 2 Figures

clear all
close all
clc

% === GLOBAL FONT SETTINGS ===
set(groot, 'DefaultAxesFontName', 'Helvetica');
set(groot, 'DefaultTextFontName', 'Helvetica');
set(groot, 'DefaultAxesFontSize', 18);
set(groot, 'DefaultAxesFontWeight', 'bold');
set(groot, 'DefaultAxesLineWidth', 2);
set(groot, 'DefaultAxesTickLength', [0.01 0.01]);
set(groot, 'DefaultAxesTickDir', 'in');  % or 'out' if preferred

% =========================================================================
% SPECTROGRAMS WITH STRAIN (SWELL-BAND ONLY)
% =========================================================================
% Purpose:
% This script processes 1 Hz seismic data (long merged SAC files) to estimate 
% swell and infragravity band strain using the Uz -> σ_xx transfer function 
% from Lipovsky (2018), Eq. (6). It can be modified for extensional stress 
% and strain calculations using Eq. (11) from the same reference.
%
% Method:
% 1. Loads and preprocesses seismic SAC data.
% 2. Applies spectral analysis (spectrogram).
% 3. Converts spectral power density to stress/strain using transfer functions.
% 4. Visualizes results in frequency and time domains.
% =========================================================================

% Set component
comp = 'LHZ';

% === Load the merged SAC file === all dates are 1 mmm 0000Z to & including 31 mmm 2300Z
%mergedfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR03/DR03_LHZ_--_2014_325to365_merged_tt.sac'; %21 Nov 2014 to 31 Dec 2014
%mergedfile='/Volumes/CATALOGDR00/DATA/Longfiles2/DR03/DR03_LHZ_--_2015_001to365_merged_tt.sac'; %1 Jan 2015 to 31 Dec 2015 
mergedfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR03/DR03_LHZ_--_2016_001to314_merged_tt.sac'; %1 Jan 2016 to 10 Nov 2016

%read-in sac file & assign data and header (there are no headers here)
[data, hdr] = rdsac(mergedfile);

% ===== Preprocessing ======
% remove trend, use sensitivity to convert to velocity
data=detrend(data);
sensitivity = 5.038830e8;
data=data/sensitivity; % 'data' output is in velocity

% Option: decimate
D=1; %Decimation Factor (not using rn)
% adata=decimate(adata,D); % adata is for (data conv to) acceleration
% ddata=decimate(data,D); % ddata is decimated data 

% === Build time vector ===
dt = 1; % sample interval in seconds
n = length(data);
t = (0:n-1) * dt; % time in seconds

% === Generate datetime vector for table ===
%start_datetime = datetime(2014, 11, 21, 0, 0, 0);  % Start at midnight Jan 1, 2015
%start_datetime = datetime(2015, 1, 1, 0, 0, 0);  % Start at midnight Jan 1, 2015
start_datetime = datetime(2016, 1, 1, 0, 0, 0);  % Start at midnight Jan 1, 2016
time_vec = start_datetime + seconds(0:n-1);     % One timestamp per sample

% === Combine into table ===
T = table(time_vec(:), data(:), ...
    'VariableNames', {'Time_UTC', 'Velocity_mps'});
% writetable(T, 'DR01_2015_Velocity_with_Time.csv');
% '01-Jan-2015 00:00:00' is the format for Time_UTC, column 1 of T

% === Define quarters for 2014 ===
%data_q1 = [datetime(2014,11,21,0,0,0), datetime(2014,12,31,23,59,59)];

% === Define quarters for 2015 ===
% data_q1 = [datetime(2015,1,1,0,0,0), datetime(2015,3,31,23,59,59)]; %90 days
% data_q2 = [datetime(2015,4,1,0,0,0), datetime(2015,6,30,23,59,59)]; %91 days
% data_q3 = [datetime(2015,7,1,0,0,0), datetime(2015,9,30,23,59,59)]; %92 days
% data_q4 = [datetime(2015,10,1,0,0,0), datetime(2015,12,31,23,59,59)]; %92 days

% === Define quarters for 2016 ===
data_q1 = [datetime(2016,1,1,0,0,0),  datetime(2016,3,31,23,59,59)]; %
data_q2 = [datetime(2016,4,1,0,0,0),  datetime(2016,6,30,23,59,59)]; %
data_q3 = [datetime(2016,7,1,0,0,0),  datetime(2016,9,30,23,59,59)]; %
data_q4 = [datetime(2016,10,1,0,0,0), datetime(2016,11,10,23,59,59)]; %

% === SELECT QUARTER TO ANALYZE ===
selected_q = data_q3;  % <--- change this to q1, q2, q3, or q4

% Index into table T to get the selected quarter
T_q = T(T.Time_UTC >= selected_q(1) & T.Time_UTC <= selected_q(2), :);
data = T_q.Velocity_mps; % Extract the data and time vector for spectrogram
time_vec = T_q.Time_UTC; 
datacheck = length(data)/86400;
Td = datenum(T_q.Time_UTC);
xl = [datenum(selected_q(1)), datenum(selected_q(2) + seconds(1))];
%xl = datenum(selected_q); % x-axis limits for plots, matching selected quarter
%for example, xl = [735870 735961]  % i.e., '01-Jul-2015' to '30-Sep-2015'
% === Generate dynamic xticks for the selected quarter ===
xtick_dates = linspace(xl(1), xl(2), 7);  % 7 ticks spaced evenly across the quarter

% Td vector with manually set dates and no table indexing
%Td = datenum(selected_q(1)) + (0:length(data)-1) * delta / 86400;

%  === Quarter labels  ===
quarter_start = selected_q(1);
quarter_end = selected_q(2);
quarter_label = sprintf('%sTo%s', datestr(quarter_start, 'ddmmmyyyy'), datestr(quarter_end, 'ddmmmyyyy'));

% Update output filenames
%2014
% year = '2014';
% file_base_name = sprintf('Strain_PaGram_DR03_2014_%s_M', quarter_label);
% file_base_name2 = sprintf('Strain_DR03_2014_%s_M', quarter_label);

%2015
% year = '2015';
% file_base_name = sprintf('Strain_PaGram_DR03_2015_%s_M', quarter_label);
% file_base_name2 = sprintf('Strain_DR03_2015_%s_M', quarter_label);

%2016
year = '2016';
file_base_name = sprintf('Strain_PaGram_DR03_2016_%s_M', quarter_label);
file_base_name2 = sprintf('Strain_DR03_2016_%s_M', quarter_label);

% === Variables manual input to make main code run ===
nstas = 1; % number of stations = 1
station = 'DR03';
Fs = 1;
delta = 1; 

% =========================================================================
% MAIN PROCESSING LOOP
% =========================================================================

for i=1:nstas
    H=figure(i);

% =========================================================================
% SPECTROGRAM COMPUTATION. Option 1: Chen 2019; Option 2: Lipovsky 2018
% =========================================================================

    % Spectrogram window length
    NW=2048;
    freqs=[];

    Ax=[110 130 1000 1200];
    set(H,'Position',Ax);
    [S,F,T,PSD]=spectrogram(data,hanning(NW),round(0.99*NW),freqs,Fs);
    %[S,F,T,PSD]=spectrogram(ddata,hanning(NW),round(0.99*NW),freqs,Fs);
    %[S,F,T,PSD]=spectrogram(detrend(adata),1e4,[],[],1);
    %F=F/sachdr.delta/D; % Normalize frequency axis
    F=F/delta/D; % Normalize frequency axis
    [~,C]=size(PSD); 

% =======
% Convert spectrogram to stress consistent with Lipovsky (2018)
% =======

    %instrument lower corner corrections
    f_ins=1/120; % Instrument response correction frequency
    T_ins=1-(1./(1+(F/f_ins).^2)); % Instrument correction?
    T_ins(1)=inf; % Zero frequency fix
    
    % Ice and ocean properties
    h = 220; % Ice thickness (m)
    H = 489; % Ocean cavity thickness (m)
    rhoi = 916; % Ice density (kg/m³)
    rhow = 1028; % Water density (kg/m³)
    g = 9.8; % Gravity (m/s²)
    nu = 0.33; % Poisson's ratio
    E = 8.7e9; % Elastic modulus (Pa)
    Ep = E/(1-nu^2); % Effective elastic modulus
    D = Ep/12*h^3; % Flexural rigidity

    % Wavenumber and phase speed computation (Lipovsky, Eq. 4)
    kmax = log10(2/h);
    k = logspace(-10,kmax,100);
    w =  sqrt( (D*k.^5/rhow + g*k) ./ (k*h*rhoi/rhow + coth(k*H)) );
    c = w./k;
    
    % Equation 6 in Lipovsky (2018)
    % Compute transfer function from Uz to stress (Lipovsky, Eq. 6)
    Z = -1j*Ep*h*k./c;
    f_xfer = w/(2*pi); % Convert to frequency domain
    
    Zi = interp1(f_xfer,Z,F,'spline'); % Interpolated transfer function from vertical velocity to extensional strain 
    T_stress=abs(Zi)./(T_ins); %Amplitude transfer function from vertical velocity to extensional strain, with instrument correction
    Tmat_stress=repmat(T_stress,1,C); 
    PSD_s=PSD.*Tmat_stress.^2; % Convert PSD to stress PSD
 
    % Smooth the PSD spectrograms a bit (quasi-Gaussian)
    PSDsm=imgaussfilt(PSD_s,[3,11]);

% =======
% Figure 1: Stress, Strain, Spectrogram
% =======

    H=figure(1);    
    % Subplots 3-4: Spectrogram with PSD (rel to 1 Pa^2/Hz) 
    h=subplot(3,1,1:2);
    imagesc(Td,F,10*log10(PSDsm+eps));
    xlim(xl)
    %datetick('x',3)
    set(gca,'xtick',[])
    axis xy
    %bookfonts
    ylabel('F (Hz)')
    colormap(jet)
    ax=get(h,'Position');
    %s=mean(std(10*log10(PSDsm)));
    %m=median(median(10*log10(PSDsm)));
    %caxis([m-3*s m+4*s])
    %caxis(PSDlim)
    hc=colorbar;
    set(get(hc,'label'),'string','PSD (dB rel. 1 Pa^2/Hz)');
    set(h,'Position',ax);
    %axis tight
    ylim([0 0.15]) %undecimated
    %ylim([0 0.12]) %for decimated D=4
    caxis([50 100])
    
    % tick_locations=[datenum(2015,1,1:365),datenum(2016,1,1:365)];
    % set(gca,'XMinorTick','on','YMinorTick','on')
    % set(gca,'XTick',tick_locations)
    % datetick('x',6,'keeplimits', 'keepticks')
    % xtickangle(90)
    title([station, ' ', year, ' Pressure Spectrogram, 0.03–0.12 Hz Bandpass'])

    %axis tight
    set(h,'Position',ax)
    set(gcf,'units','normalized','position',[0.1300 0.1100 0.8 0.7])
    
    Td_s=Td;
    N=numel(data);
    if N/2 ~= round(N/2)
        N=N-1;
        data=data(1:end-1);
        Td_s=Td(1:end-1);
    end

% =========================================================================
% STRESS TIME SERIES COMPUTATION- mid-spectrogram.
% =========================================================================
    
    T_ins_xfer=1-(1./(1+(f_xfer/f_ins).^2));
    f_spec = Fs*(0:(N/2))/N;
    nyquist = 0.5;
    
    %Swell band
    [b,a]=butter(2,[0.03/nyquist, 0.12/nyquist],'bandpass');
    data_swell_spec=fft(filtfilt(b,a,data));
    
    %IG band
    %[b,a]=butter(2,[0.005/nyquist, 0.02/nyquist],'bandpass');
    %data_ig_spec=fft(detrend(filtfilt(b,a,data)));

    %Highpass band
    %[b,a]=butter(2,[0.05/nyquist, 0.49/nyquist],'bandpass');
    %data_hp_spec=fft(detrend(filtfilt(b,a,data)));
    
    Z_spec = interp1(f_xfer,Z./T_ins_xfer,f_spec)'; % Interpolate transfer function
    Z_spec(1) = 0;
    
    Tf_stress = zeros(size(data_swell_spec));
    Tf_stress(2:(N/2+1) ) = Z_spec(2:end);              % + frequencies
    Tf_stress(N/2 + 2:end) = fliplr( Z_spec(3:end));   % - frequencies

    data_stress_swell = real(ifft(data_swell_spec.* Tf_stress));
    %data_stress_ig = real(ifft(data_ig_spec.*Tf_stress));
    %data_stress_hp = real(ifft(data_hp_spec.*Tf_stress));

% =======Return to figure 1 
% Figure 1: Stress, Strain, Spectrogram,Swell-band Stress-strain:
% ======= 
    % Subplot 2: Swell-band stress (left y-axis) & strain (right y-axis) 
    subplot(3,1,3)
    
    % Swell-band stress (left y-axis)
    yyaxis left
    % plot(Td_s,data_stress_swell,'b-','linewidth',0.1);
    % hold on   
    plot(Td_s,movmax(abs(data_stress_swell),3600),'k-','linewidth',2);
    hold on
    %plot(1200*ones(size(Td_s)),'k--','linewidth',2);
    %plot(3800*ones(size(Td_s)),'k--','linewidth',2);
    hold off
    ylim(1e4*[-.1 1])
    ylabel('\sigma_x (Pa)')
    
    % Swell-band strain (right y-axis) 
    yyaxis right
    ylim(1e4/Ep*[-.1 1])
    ylabel('\epsilon_x')
    %xlim(xl)
    %bookfonts
    %title([sachdr.kstnm(1:4),' (0.03 - 0.12 Hz)'])
    %set(gca,'xtick',[])
    title([station, ' ', year, ' Vertical Velocity Stress-Strain, 0.03–0.12 Hz Bandpass'])
    ax=gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';

    xlim(xl)
    set(gca, 'XTick', xtick_dates)
    datetick('x','mmm dd','keeplimits','keepticks')
    xtickangle(0)
    xlabel('Date')
    % 
    % xlim(xl)
    % xtick_locations = datenum(2015,1,1):15:datenum(2015,4,1);
    % %set(gca, 'XTick', xtick_locations)
    % set(gca, 'XTick', xtick_dates)
    % datetick('x','mmm dd','keeplimits'); %,'keepticks')
    % %datetick('x', 'mmm dd', 'keeplimits', 'keepticks')
    % xtickangle(0)
    

% =======
% Figure 2 of only Swell-band Stress-strain: 
% =======

G=figure(2);
    % Subplot 2: Swell-band stress (left y-axis) & strain (right y-axis) 
    subplot(1,1,1)
    
    % Swell-band stress (left y-axis)
    yyaxis left
    % plot(Td_s,data_stress_swell,'b-','linewidth',0.1);
    % hold on   
    plot(Td_s,movmax(abs(data_stress_swell),3600),'k-','linewidth',2);
    hold on
    %plot(1200*ones(size(Td_s)),'k--','linewidth',2);
    %plot(3800*ones(size(Td_s)),'k--','linewidth',2);
    hold off
    ylim(1e4*[-.1 1])
    ylabel('\sigma_x (Pa)')
    
    % Swell-band strain (right y-axis) 
    yyaxis right
    ylim(1e4/Ep*[-.1 1])
    ylabel('\epsilon_x')
    %xlim(xl)
    %bookfonts
    %title([sachdr.kstnm(1:4),' (0.03 - 0.12 Hz)'])
    %set(gca,'xtick',[])
    title([station, ' ', year, ' Vertical Velocity Stress-Strain, 0.03–0.12 Hz Bandpass'])
    %title([station(1), ' ', year, ' Vertical Velocity Stress-STrain, 0.03-0.12 Hz Bandpass']) % for no SAC header
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    set(gcf,'units','inches','position',[1 1 20 3])  % Also 5x2 inches

    xlim(xl)
    set(gca, 'XTick', xtick_dates)
    datetick('x','mmm dd','keeplimits','keepticks')
    xtickangle(0)
    xlabel('Date')
end
%%
% save figs
    % === SAVE AFTER EACH FIGURE ===
    folder_path = '/Volumes/CATALOGDR01/2025FIGS/'; 
    file_ext = '.png';
    file_ext = '.fig';
    dpi = 600;
    %fullFilePath = fullfile(folder_path, [file_base_name file_ext]);
    %fullFilePath2 = fullfile(folder_path, [file_base_name2 file_ext]);
    %set(H, 'Units', 'inches');
    %set(H, 'Position', [1, 1, 30, 2]);  %[1, 1, 5, 2]

    % Save Fig 1 as PNG
    pngFilePath = fullfile(folder_path, [file_base_name, '.png']);
    print(H, pngFilePath, '-dpng', ['-r' num2str(dpi)]);
    
    % Save Fig 1 as MATLAB FIG
    figFilePath = fullfile(folder_path, [file_base_name, '.fig']);
    saveas(H, figFilePath); 

    % Save Fig 1 as PNG
    pngFilePath2 = fullfile(folder_path, [file_base_name2, '.png']);
    print(G, pngFilePath2, '-dpng', ['-r' num2str(dpi)]);
    
    % Save Fig 1 as MATLAB FIG
    figFilePath2 = fullfile(folder_path, [file_base_name2, '.fig']);
    saveas(G, figFilePath2); 

%% NEXT:
% spectrogram_data_long_linF_E_2025_v2.m

%% unused/removed code -->

























% % ~~
% % IG-band Stress-strain:     
% 
%     % Subplot 3: IG-band stress (left y-axis) & strain (right y-axis) 
%     subplot(5,1,3)
% 
%     % detrending infragravity band stress- not sure why.
%     data_stress_ig=detrend(data_stress_ig-movmean(data_stress_ig,3600));
% 
%     % IG-band stress (left y-axis)
%     yyaxis left
%     %plot(Td_s,detrend(data_stress_ig),'b-','linewidth',0.1);
%     %hold on
%     plot(Td_s,movmax(abs(data_stress_ig),3600),'k-','linewidth',2);
%     hold on
%     %plot(1200*ones(size(Td_s)),'k--','linewidth',2);
%     %plot(3800*ones(size(Td_s)),'k--','linewidth',2);
%     hold off
%     ylim(1e4*[-.1 1])
%     ylabel('\sigma_x (Pa)')
% 
%     % IG-band strain (right y-axis) 
%     yyaxis right
%     ylim(1e4/Ep*[-1 1])
%     ylabel('\epsilon_x')
%     xlim(xl)
%     bookfonts
%     title([sachdr.kstnm(1:4),' (0.005 - 0.02 Hz)'])
%     set(gca,'xtick',[])
%     ax=gca;
%     ax.YAxis(1).Color = 'k';
%     ax.YAxis(2).Color = 'k';
%     set(gca,'xtick',[])

% ~~
% % Hp-band Stress-strain: 
% 
%     % Subplot 1: Hp-band stress (left y-axis) & strain (right y-axis) 
%     subplot(5,1,1)
% 
%     % does the hp band need to be detrended too? swell-band wasn't.
%     data_stress_hp=detrend(data_stress_hp-movmean(data_stress_hp,3600));
% 
%     % Hp-band stress(left y-axis) 
%     yyaxis left
%     %plot(Td_s,detrend(data_stress_ig),'b-','linewidth',0.1);
%     %hold on
%     plot(Td_s,movmax(abs(data_stress_hp),3600),'k-','linewidth',2);
%     hold on
%     %plot(1200*ones(size(Td_s)),'k--','linewidth',2);
%     %plot(3800*ones(size(Td_s)),'k--','linewidth',2);
%     hold off
%     ylim(1e4*[-.1 1])
%     ylabel('\sigma_x (Pa)')
% 
%     % Hp-band strain (right y-axis) 
%     yyaxis right
%     ylim(1e4/Ep*[-1 1])
%     ylabel('\epsilon_x')
%     xlim(xl)
%     bookfonts
%     title([sachdr.kstnm(1:4),' (0.05 - 0.5 Hz)'])
%     set(gca,'xtick',[])
%     ax=gca;
%     ax.YAxis(1).Color = 'k';
%     ax.YAxis(2).Color = 'k';
%     set(gca,'xtick',[])

% =======
% Figure 2: Stress, Flexural mode impedance component (Zxxz) vs Frequency (Hz)
    % Equation 6

    % figure(2)
    % Fplot=logspace(-3,-1,1000);
    % Zinterp=interp1(F,abs(Zi),Fplot,'spline');
    % semilogx(Fplot,Zinterp,'linewidth',2)
    % xlabel('F (Hz)')
    % ylabel('Z^F_{XXZ} (Pa/m/s)')
    % grid on
    % bookfonts
    % xlim([1e-3 .1])

% =======
% Figure 3: Wave phase velocity (m/s) vs Frequency (Hz)     

    % figure(3)
    % wplot=2*pi*Fplot;
    % cinterp=interp1(w,c,wplot,'spline');
    % semilogx(wplot/(2*pi),cinterp,'linewidth',2)
    % xlabel('F (Hz)')
    % ylabel('c (m/s)')
    % bookfonts
    % xlim([.001 .1])
    % grid on
    % 

%end

% added by EAM
% total_season_stress_swell = trapz(abs(data_stress_swell));
% total_season_stress_hp = trapz(abs(data_stress_hp));
% total_season_stress_ig = trapz(abs(data_stress_ig));
% 
% med_season_stress_swell = median(abs(data_stress_swell));
% med_season_stress_hp= median(abs(data_stress_hp));
% med_season_stress_ig = median(abs(data_stress_ig));
% 
% mean_season_stress_swell = mean(abs(data_stress_swell));
% mean_season_stress_hp = mean(abs(data_stress_hp));
% mean_season_stress_ig = mean(abs(data_stress_ig));
% 
% std_stress_swell = std(abs(data_stress_swell));
% std_stress_hp = std(abs(data_stress_hp));
% std_stress_ig = std(abs(data_stress_ig));
% 
% max_stress_swell = max(abs(data_stress_swell));
% max_stress_hp = max(abs(data_stress_hp));
% max_stress_ig = max(abs(data_stress_ig));
% 
% min_stress_swell = min(abs(data_stress_swell));
% min_stress_hp = min(abs(data_stress_hp));
% min_stress_ig = min(abs(data_stress_ig));
% 
% % table 
% stress_table = table(...
%     total_season_stress_swell, total_season_stress_hp, total_season_stress_ig, ...
%     med_season_stress_swell, med_season_stress_hp, med_season_stress_ig, ...
%     mean_season_stress_swell, mean_season_stress_hp, mean_season_stress_ig, ...
%     std_stress_swell, std_stress_hp, std_stress_ig, ...
%     max_stress_swell, max_stress_hp, max_stress_ig, ...
%     min_stress_swell, min_stress_hp, min_stress_ig, ...
%     'VariableNames', {...
%         'Total_Stress_Swell', 'Total_Stress_HP', 'Total_Stress_IG', ...
%         'Median_Stress_Swell', 'Median_Stress_HP', 'Median_Stress_IG', ...
%         'Mean_Stress_Swell', 'Mean_Stress_HP', 'Mean_Stress_IG', ...
%         'Std_Stress_Swell', 'Std_Stress_HP', 'Std_Stress_IG', ...
%         'Max_Stress_Swell', 'Max_Stress_HP', 'Max_Stress_IG', ...
%         'Min_Stress_Swell', 'Min_Stress_HP', 'Min_Stress_IG'});


% Save the table to a CSV file
%writetable(stress_table, 'stress_metrics_DR01_2015.csv');

%% Notes

%set up for 1 Hz data long (merged) sac files
%estimates swell and infragravity band strain using the Lipovsky (2018)
%Uz -> sigma_xx transfer function; equation (6)
%could be modified for extensional stress and strain as well using equation
%(11) in the same reference given a velocity time series
%

% trapz: Q = trapz(Y) computes the approximate integral of Y via the 
% trapezoidal method with unit spacing. The size of Y determines the 
% dimension to integrate along:
% If Y is a vector, then trapz(Y) is the approximate integral of Y.


%% Stress Histogram
% 
% num_bins = 100; 
% 
% figure;
% histogram(abs(data_stress_swell), num_bins, 'Normalization', 'pdf', 'FaceColor', [0.3010 0.7450 0.9330], 'EdgeColor', 'b');
% hold on;
% %[f, xi] = ksdensity(data_stress_swell); 
% y_limits = ylim; % Get y-axis limits for vertical lines
% xlim([0 2000])
% xlabel('Stress Value');
% ylabel('Probability Density');
% title('Distribution of Stress Values, DR01 1/1-3/16 2015');
% set(gca, 'FontSize', 12);
% grid on;
% hold off;

