%% Merged file, Seismogram and spectrogram with the longfile 1 Hz data
clear all
close all
clc
%%
%addpath('/Users/em/PROJECTS/PEAKTEM/Strain/Elisa_Stress_Strain/DATA')

% === GLOBAL FONT SETTINGS ===
set(groot, 'DefaultAxesFontName', 'Helvetica');
set(groot, 'DefaultTextFontName', 'Helvetica');
set(groot, 'DefaultAxesFontSize', 18);
set(groot, 'DefaultAxesFontWeight', 'bold');
set(groot, 'DefaultAxesLineWidth', 2);
set(groot, 'DefaultAxesTickLength', [0.01 0.01]);
set(groot, 'DefaultAxesTickDir', 'in');  % or 'out' if preferred

% Set component
comp = 'LHZ';

% === Load the merged SAC file === all dates are 1 mmm 0000Z to & including 31 mmm 2300Z
%mergedfile='/Volumes/CATALOGDR00/DATA/Longfiles2/DR02/DR02_LHZ_--_2014_325to365_merged_tt.sac'; %21 Nov 2014 to 31 Dec 2014
%mergedfile='/Volumes/CATALOGDR00/DATA/Longfiles2/DR02/DR02_LHZ_--_2015_001to365_merged_tt.sac'; %1 Jan 2015 to 31 Dec 2015 
mergedfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR02/DR02_LHZ_--_2016_001to314_merged_tt.sac'; %1 Jan 2016 to 10 Nov 2016

%read-in sac file & assign data and header (there are no headers here)
[data, hdr] = rdsac(mergedfile);

% ===== Preprocessing ======
% remove trend, use sensitivity to convert to velocity
data=detrend(data);
sensitivity = 5.038830e8;
data=data/sensitivity; % 'data' output is in velocity
adata=diff([0;data]);  %convert to acceleration (1 s/s)

% Option: decimate
D=1; %Decimation Factor (not using rn)
% adata=decimate(adata,D); % adata is for (data conv to) acceleration
% ddata=decimate(data,D); % ddata is decimated data 

% === Build time vector ===
dt = 1; % sample interval in seconds
n = length(data);
na = length(adata);
%n = length(adata);
t = (0:n-1) * dt; % time in seconds

% === Generate datetime vector for table ===
%start_datetime = datetime(2014, 11, 21, 0, 0, 0);  % Start at midnight Jan 1, 2015
%start_datetime = datetime(2015, 1, 1, 0, 0, 0);  % Start at midnight Jan 1, 2015
start_datetime = datetime(2016, 1, 1, 0, 0, 0);  % Start at midnight Jan 1, 2016
time_vec = start_datetime + seconds(0:n-1);     % One timestamp per sample

% === Velocity Combine into table ===
%T = table(time_vec(:), data(:), ...
%    'VariableNames', {'Time_UTC', 'Velocity_mps'});
% writetable(T, 'DR01_2015_Velocity_with_Time.csv');
% '01-Jan-2015 00:00:00' is the format for Time_UTC, column 1 of T

% === Acceleration Combine into table ===
% T = table(time_vec(:), adata(:), ...
%     'VariableNames', {'Time_UTC', 'Accel_mpss'});

% === Acceleration and Velocity Combine into table ===
T = table(time_vec(:), data(:), adata(:), ...
    'VariableNames', {'Time_UTC', 'Velocity_mps', 'Accel_mpss'});

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
selected_q = data_q1;  % <--- change this to q1, q2, q3, or q4

% Index into table T to get the selected quarter
T_q = T(T.Time_UTC >= selected_q(1) & T.Time_UTC <= selected_q(2), :);
%data = T_q.Velocity_mps; % data, Velocity Extract the data and time vector for spectrogram
adata = T_q.Accel_mpss; % adata, Accel, Extract the data and time vector for spectrogram
data = T_q.Velocity_mps; % data, vel, Extract the data and time vector for velo plot
time_vec = T_q.Time_UTC; %
%datacheck = length(data)/86400;
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
% file_base_name = sprintf('SpectrogramSwell1Hz_DR02_2014_%s_M', quarter_label);
% file_base_name2 = sprintf('Velocity1Hz_DR02_2014_%s_M', quarter_label);

%2015
% year = '2015';
% file_base_name = sprintf('SpectrogramSwell1Hz_DR02_2015_%s_M', quarter_label);
% file_base_name2 = sprintf('Velocity1Hz_DR02_2015_%s_M', quarter_label);

%2016
year = '2016';
file_base_name = sprintf('SpectrogramSwell1Hz_DR02_2016_%s_M', quarter_label);
file_base_name2 = sprintf('Velocity1Hz_DR02_2016_%s_M', quarter_label);

% === Variables manual input to make main code run ===
nstas = 1; % number of stations = 1
station = 'DR02';
Fs = 1;
delta = 1; 
ny = Fs/2; % nyquist is 0.5 Hz

% =========================================================================
% MAIN SPECTROGRAM PROCESSING LOOP
% =========================================================================

for i=1:nstas
    
    H=figure(i);

    low_cutoff = 0.03;  % Lower cutoff frequency in Hz
    high_cutoff = 0.12;  % Upper cutoff frequency in Hz
    Wn = [low_cutoff high_cutoff] / ny; % Normalize frequencies by Nyquist frequency (Fs/2)
    [b, a] = butter(2, Wn, 'bandpass');  % 2nd-order Butterworth bandpass filter
    %NW=4096;
    %NW=2048;
    NW=1024;
    %NW=512; %nov-dec 2014

    freqs=[];

    [S,F,T,PSD]=spectrogram(adata,hanning(NW),round(0.95*NW),freqs,1);
    %[S,F,T,PSD]=spectrogram(data,hanning(NW),round(0.95*NW),freqs,1);
    %[S,F,T,PSD]=spectrogram(detrend(fdata),hanning(NW),round(0.95*NW),freqs,1);
    %[S,F,T,PSD]=spectrogram(detrend(adata),1e4,[],[],1);
    %F=F/sachdr.delta/D;
    F=F/delta/D;

    % === MEAN, SMOOTHING =======================

    % medPSD=mean(PSD,2);
    % [~,C]=size(PSD); %build a matrix for pre-event spectral normalization
    % PSDnorm=repmat(medPSD,1,C);
    % PSDsm=imgaussfilt(PSD,[3,11]); % Smoothed PSD, PSDsm
    % PSDdiff=PSD./PSDnorm; % Normalized PSD (rel med), PSDdiff
    % PSDdiffsm=PSDdiff; % Smoothed normalized PSD (rel med), PSDdiffsm
    % %PSDdiffsm=imgaussfilt(PSDdiff,[0.5,5]);
    % PSDdiffsm=imgaussfilt(PSDdiff,[1.5,11]);

    % === NO MEAN, NO SMOOTHING, FOR COMPARING ACROSS THE YEAR ========

    PSDsm=imgaussfilt(PSD,[3,11]); % Smoothed PSD, PSDsm

    % === MAIN PLOT =========================
    naa = length(adata);
    %PSDlim=[-160 -80];
    
    h=subplot(2,1,1:2);
    imagesc(Td,F,10*log10(PSDsm)); %STANDARDIZED FOR X-SEASON COMPARE
    %imagesc(Td,F,10*log10(PSDdiffsm));
    %imagesc(Td,F,10*log10(PSDdiff));
    %xlim(xl)
    axis xy
    ylabel('F (Hz)'); %, 'FontSize', 20) % Increase Y-axis label size
    colormap(jet)
 
    ax=get(h,'Position');
    %s=mean(std(10*log10(PSDdiffsm)));
    %m=median(median(10*log10(PSDdiffsm+eps)));
    %s=mean(std(10*log10(PSDdiff)));
    %m=median(median(10*log10(PSDdiff+eps)));
    %caxis([m-3*s m+3*s])

    % y-axis stuff
    ylim([0.03 0.12]) % swell band high limit
    yticks([0.04 0.06 0.08 0.10 0.12])  % include 0.12 as a tick

    % c-axis stuff
    %caxis([-35 15]) %for PSD rel med
    caxis([-160 -50]) % for non-rel med
    hc=colorbar;
    %set(get(hc, 'label'), 'string', 'PSD (dB rel. med. 1 (m/s)^{2})/Hz)') % velocity power
    %set(get(hc, 'label'), 'string', 'PSD (dB rel. med. 1 (m/s^{2})^{2})/Hz)') % acceleration power
    set(get(hc, 'label'), 'string', 'PSD (dB rel. 1 (m/s^{2})^{2})/Hz)') % acceleration power
    set(hc); %, 'FontSize', 20) % Increase colorbar tick size
    %set(h,'Position',ax);

    % figure title stuff and size
    %title([station, ' ', year, ' Velocity Spectrogram, 0.03–0.12 Hz Bandpass'])
    title([station, ' ', year, ' Acceleration Spectrogram, 0.03–0.12 Hz Bandpass'])
    set(h,'Position',ax)
    %set(gcf,'units','normalized','position',[0.1300 0.1100 0.8 0.7])
    %set(gcf,'units','inches','position',[1, 1, 18, 4])  % Also 5x2 inches
    %set(gcf,'units','inches','position',[1, 1, 6.15, 0.65])
    set(gcf,'units','inches','position',[1, 1, 20, 4])

    % x-axis stuff
    xlim(xl)
    set(gca, 'XTick', xtick_dates)
    datetick('x','mmm dd','keeplimits','keepticks')
    xtickangle(0)
    xlabel('Date')

    G = figure(3)
    plot(Td, data, 'k', LineWidth=2);
    xlabel('Time (UTC)');
    %ylabel('Amplitude (counts)');
    ylabel('Velocity (m/s)');
    title('Merged 1 Hz SAC Data: DR01 LHZ');
    grid off;
    ylim([-4e-3 4e-3])
    axis xy
    title([station, ' ', year, ' Velocity, 0.03–0.12 Hz Bandpass'])
    set(h,'Position',ax)
    %set(gcf,'units','normalized','position',[0.1300 0.1100 0.8 0.7])
    %set(gcf,'units','inches','position',[1, 1, 18, 4])  % Also 5x2 inches
    %set(gcf,'units','inches','position',[1, 1, 6.15, 0.65])
    set(gcf,'units','inches','position',[1, 1, 20, 3])

    % x-axis stuff
    xlim(xl)
    set(gca, 'XTick', xtick_dates)
    datetick('x','mmm dd','keeplimits','keepticks')
    xtickangle(0)
    xlabel('Date')

        % === SAVE AFTER EACH FIGURE ===
    folder_path = '/Volumes/CATALOGDR01/2025FIGS/'; 
    %file_base_name = sprintf('DR01_2016_1Janto1Apr_1024_sm_Swell_%02d', i);
    %file_base_name = sprintf('DR01_2016_1Janto1Apr_1024_sm_Swell2');
    file_ext = '.png';
    file_ext = '.fig';
    dpi = 600;
    fullFilePath = fullfile(folder_path, [file_base_name file_ext]);

    %set(H, 'Units', 'inches');
    %set(H, 'Position', [1, 1, 30, 4]);  %[1, 1, 5, 2]

    % Save as PNG
    pngFilePath = fullfile(folder_path, [file_base_name, '.png']);
    pngFilePath2 = fullfile(folder_path, [file_base_name2, '.png']);
    %print(H, pngFilePath, '-dpng', ['-r' num2str(dpi)]);
    %print(G, pngFilePath2, '-dpng', ['-r' num2str(dpi)]);
    
    % Save as MATLAB FIG
    %figFilePath = fullfile(folder_path, [file_base_name, '.fig']);
    %figFilePath2 = fullfile(folder_path, [file_base_name2, '.fig']);
    %saveas(H, figFilePath);  % Alternatively: savefig(H, figFilePath);
    %saveas(G, figFilePath2); 

end


%% old way of saving:
%folder_path = '/Volumes/CATALOGDR01/200HzSpectro/';
%matfilesname = 'DR01_2016_1Janto1Apr_1024_sm_Swell.png';
%fullFilePath = fullfile(folder_path, matfilesname);
%saveas(H, fullFilePath);


