clear all
close all
clc

% spectrogram data long linF stress strain EAM v2 SLOPE.m
% slpe 1: 4.236781186270549e-07
% slpe 2: 4.236781186270552e-07

%% 1: figure of selected slope points in PSD max from spectrogram
% =========================================================================
% Spectrogram, selected section of time series, and resulting slope
% =========================================================================
% ~ 3 min to run full  script

clear
fclose('all');

startTime = datetime('now'); % or use datestr(now) for a string version
disp(['Script started at: ', datestr(startTime)])
tic; % Start timing

comp = 'LHZ';

% Set the path to your data folder
%data_folder = '/Volumes/CATALOGDR00/Longfiles';system(['/bin/ls ', data_folder, '/DR01-2015-0-xx_', comp, '.sac > list.all']);
data_folder = '/Volumes/CATALOGDR00/Longfiles';system(['/bin/ls ', data_folder, '/DR01-2016-0-xx_', comp, '.sac > list.all']);
!cat list.all | wc -l > list.num

% Set date range and component
%xl=[datenum('01-01-2015'),datenum('03-31-2015')]; % 03-16 0000Z (incl 3-15)
%comp='LHZ';

% List all SAC files to process
%eval(['!/bin/ls DATA/DR01-2015-0-xx_',comp,'.sac > list.all'])
%!cat list.all | wc -l > list.num

% Open file lists
fid1 = fopen('list.num');
fid2 = fopen('list.all');

% Read file names into an array
i=1;
while 1
           tline = fgetl(fid2);
           
            if ~ischar(tline), break, end
            if i==1
            stnames=tline;
            else
                stnames=char(stnames,tline);
            end
            i=i+1;
        end
        fclose(fid2);

% Read the number of stations/files        
nstas=fscanf(fid1,'%d');
fclose(fid1);

% =========================================================================
% LOAD, CONVERT TO VELOCITY, DECIMATE OPTION, DATE-TIME 
% =========================================================================

for i=1:nstas

    % === GLOBAL FONT SETTINGS ===
    set(groot, 'DefaultAxesFontName', 'Helvetica');
    set(groot, 'DefaultTextFontName', 'Helvetica');
    set(groot, 'DefaultAxesFontSize', 18);
    set(groot, 'DefaultAxesFontWeight', 'bold');
    set(groot, 'DefaultAxesLineWidth', 2);
    set(groot, 'DefaultAxesTickLength', [0.01 0.01]);
    set(groot, 'DefaultAxesTickDir', 'in');  % or 'out' if preferred

    H=figure(i);
    %set(H, 'Units', 'inches');         % Set units to inches
    %set(H, 'Position', [1, 1, 10, 4]);   % [x, y, width, height] in inches

    % Load SAC file data and header    
    [sachdr,data]=load_sac(stnames(i,:));

    % Preprocess: remove trend, use sensitivity to convert to velocity
    data=detrend(data);
    sensitivity = 5.038830e8;
    data=data/sensitivity; % 'data' output is in velocity

    % Option: convert to displacement (1 s/s)
    %data=cumsum(data);

    % Option: convert to acceleration (1 s/s)
    %adata=diff([0;data]);

    % Option: decimate
    D=1; %Decimation Factor (not using rn)
    % adata=decimate(adata,D); % adata is for (data conv to) acceleration
    % ddata=decimate(data,D); % ddata is decimated data 

    % Sampling rate (samples/s)
    Fs = 1;
    delta = 1; % no headers in this data

    %~~~~~~~
    % STORM DATES 
    % 2015
    %~~~~~~~~
   
    % tStart = datenum(2015, 1, 18, 06, 0, 0); % Storm1: Jan 18, 10:00Z
    % tEnd   = datenum(2015, 1, 19, 06, 0, 0); % Storm1: Jan 19, 06:00Z
    % xl = [datenum('01-17-2015'), datenum('01-20-2015')]; % Storm 2
    % file_base_name = sprintf('STORM_DR01_2015_18Janto21Jan_METHOD_');
    % %Storm1: 0.04-0.08
    % fMin = 0.04; fMax = 0.080;
    
    % tStart = datenum(2015, 1, 28, 02, 0, 0); % Storm2: Jan 28, 00:00Z  
    % tEnd   = datenum(2015, 1, 28, 18, 0, 0); % Storm2: Jan 29, 00:00Z  
    % xl = [datenum('01-27-2015'), datenum('01-30-2015')]; % Storm 1
    % file_base_name = sprintf('STORM_DR01_2015_27Janto31Jan_METHOD');
    % %Storm2: 0.04-0.08
    % fMin = 0.04; fMax = 0.080;
    
    % Didn't use
    % tStart = datenum(2015, 2, 6, 00, 0, 0); % Storm3: Didn't use
    % tEnd   = datenum(2015, 2, 7, 15, 0, 0);
    % xl = [datenum('02-5-2015'), datenum('02-8-2015')]; 
    % file_base_name = sprintf('STORM_DR01_2015_5Febto9Feb_METHOD_');

    % tStart = datenum(2015, 2, 19, 02, 0, 0); % Storm3 primary
    % tEnd   = datenum(2015, 2, 19, 18, 0, 0);  
    % xl = [datenum('02-18-2015'), datenum('02-21-2015')]; % Storm 3 alt
    % file_base_name = sprintf('STORM_DR01_2015_19to20Feb_METHOD_');
    % %Storm3: 0.04-0.08
    % fMin = 0.04; fMax = 0.070;

    % tStart = datenum(2015, 3, 5, 7, 0, 0);  % Storm4: March 5, 12:00Z
    % tEnd   = datenum(2015, 3, 6, 6, 0, 0);  % Storm4: March 6, 12:00Z
    % xl = [datenum('03-4-2015'), datenum('03-7-2015')]; % Storm 4
    % file_base_name = sprintf('STORM_DR01_2015_3Marto8Mar_METHOD_');
    % %Storm4: 0.04-0.08
    % fMin = 0.04; fMax = 0.080;

    %~~~~~~~
    % STORM DATES 
    % 2016
    %~~~~~~~~

    % tStart = datenum(2016, 1, 11, 0, 0, 0); % Storm1
    % tEnd   = datenum(2016, 1, 12, 18, 0, 0); % Storm1: 
    % xl = [datenum('01-10-2016'), datenum('01-13-2016')]; % Storm 1
    % file_base_name = sprintf('STORM_DR01_2016_xxx_METHOD_');
    % % %Storm1: 0.04-0.08
    % fMin = 0.04; fMax = 0.080;
    
    % tStart = datenum(2016, 1, 23, 02, 0, 0); % Storm2  
    % tEnd   = datenum(2016, 1, 23, 16, 0, 0); % Storm2
    % xl = [datenum('01-22-2016'), datenum('01-25-2016')]; % Storm 2
    % file_base_name = sprintf('STORM_DR01_2016_23to25Jan_METHOD_');
    % %Storm2: 0.045-0.08
    % fMin = 0.045; fMax = 0.060;
    
    % tStart = datenum(2016, 2, 6, 14, 0, 0); % Storm3
    % tEnd   = datenum(2016, 2, 7, 05, 0, 0);  % Storm3
    % xl = [datenum('02-5-2016'), datenum('02-8-2016')]; % Storm 3
    % file_base_name = sprintf('STORM_DR01_2016_xxx_METHOD_');
    % %Storm3: 0.04-0.08
    % fMin = 0.04; fMax = 0.080;
    
    tStart = datenum(2016, 2, 28, 9, 0, 0);  % Storm4
    tEnd   = datenum(2016, 2, 28, 21, 0, 0);  % Storm4
    xl = [datenum('02-27-2016'), datenum('03-01-2016')]; % Storm 4
    file_base_name = sprintf('STORM_DR01_2016_xxx_METHOD_');
    % %Storm4: 0.068-0.1
    fMin = 0.06; fMax = 0.088;

    % for full 2016 spectrogram 
    % tStart = datenum(2016, 1, 01, 0, 0, 0);  % 2016 full
    % tEnd = datenum(2016, 4, 01, 0, 0, 0);  % 2016 full
    % xl_full = [datenum('01-01-2016'), datenum('04-01-2016')]; % 2016 full

    % did not use
    % tStart = datenum(2016, 2, 22, 12, 0, 0); 
    % tEnd   = datenum(2016, 2, 23, 12, 0, 0);  
    % xl = [datenum('02-21-2016'), datenum('02-24-2016')]; 
    % file_base_name = sprintf('STORM_DR01_2016_xxx_METHOD_');
    % did not use

% =======
% ~~~ TRIM DATA, CREATE TIME VECTOR
% ========
    % Time vector conversion from SAC header
    [month,day]=monthday(sachdr.nzyear,sachdr.nzjday);
    tzmin=datenum(sachdr.nzyear,month,day,sachdr.nzhour,sachdr.nzmin,sachdr.sec);
    tzmax=tzmin+length(data)*sachdr.delta/(24*60*60);
    Td=linspace(tzmin,tzmax,length(data));

    % Trim data and Td to xl dates
    trimMask = (Td >= xl(1)) & (Td <= xl(2));
    
    if ~any(trimMask)
        warning('No data in selected tStart–tEnd window.');
        continue; % or break, depending on how you want to handle it
    end
    
    data = data(trimMask);
    Td = Td(trimMask);
    
    % Recalculate time bounds for trimmed data (optional but clean)
    tzmin = Td(1);
    tzmax = Td(end);

% =======
% ~~~ BAND PASS DATA
% ========    

    %=== Band-pass filter between 0 or 0.03 and 0.12 Hz ===
    f_low = 0.03;   % lower cutoff frequency in Hz
    f_high = 0.12;  % upper cutoff frequency in Hz
    Fs = 1 / delta;  % actual sampling frequency from SAC header
    nyquist = Fs / 2;

    % Design 4th-order Butterworth bandpass filter
    [b, a] = butter(4, [f_low, f_high] / nyquist, 'bandpass');

    % Apply zero-phase filtering (filtfilt avoids phase distortion)
    data = filtfilt(b, a, data);

% =========================================================================
% SPECTROGRAM COMPUTATION PSD MAX - SLOPE DATES TO DETERMINE STORM DISTANCE
% =========================================================================

    % Spectrogram window length
    NW=2048; 
    freqs=[];

    Ax=[110 130 1000 1200];
    set(H,'Position',Ax);
    [S,F,T,PSD]=spectrogram(data,hanning(NW),round(0.99*NW),freqs,Fs);
    %[S,F,T,PSD]=spectrogram(adata,hanning(NW),round(0.99*NW),freqs,Fs);
    %[S,F,T,PSD]=spectrogram(detrend(adata),1e4,[],[],1);
    F=F/sachdr.delta/D; % Normalize frequency axis
    [~,C]=size(PSD); 

% =======
% ~~~ ALTER PSD W FILTER & MEDIAN
% ========
    % alter PSD first
    
    medPSD=mean(PSD,2);

    [~,C]=size(PSD); %build a matrix for pre-event spectral normalization
    
    PSDnorm=repmat(medPSD,1,C);

    % Smooth the PSD spectrograms a bit (quasi-Gaussian) Smoothed PSD, PSDsm
    %PSDsm=imgaussfilt(PSD,[3,11]);
    %PSDsm=imgaussfilt(PSD,[1.5,11]); 
    PSDsm=imgaussfilt(PSD,[0.5,5]); 
    %PSDsm = imgaussfilt(10*log10(PSD ./ PSDnorm), [3, 11]);  % dB before smoothing

    % Normalized PSD (rel med), PSDdiff
    PSDdiff=PSD./PSDnorm;
    
    % Smoothed normalized PSD (rel med), PSDdiffsm
    PSDdiffsm=PSDdiff;
    %PSDdiffsm=imgaussfilt(PSDdiff,[0.5,5]);
    %PSDdiffsm=imgaussfilt(PSDdiff,[1.5,11]);
    %PSDdiffsm=imgaussfilt(PSDdiff,[3,11]);
    %~~ new end

% ========
% ~~~ CALCULATE SLOPE GIVEN TIMES AND PSD
% ========

    % First, convert PSD peak times and frequencies as before
    %PSD_slope = abs(PSD).^2;  
    %PSD_slope = PSD;
    PSD_slope = PSDsm; 
    %PSD_slope = PSDdiffsm;
    [maxVals, maxIdxs] = max(PSD_slope, [], 1);  % max across F for each T
    peakFreqs = F(maxIdxs); % frequency values of max    
    %peakFreqs_ofSlope = peakFreqs(1, ())

    specTimes = tzmin + T / (24*60*60); % time values in datenum
    
    specTimes = specTimes(:);    % force column vector
    peakFreqs = peakFreqs(:);    % force column vector

    mask = (specTimes >= tStart) & (specTimes <= tEnd) & ...
       (peakFreqs >= fMin) & (peakFreqs <= fMax);

    fitTimes = specTimes(mask);
    fitFreqs = peakFreqs(mask);

% ========
% ~~~ POLYFIT OF SELECTED (MASKED) STORM TIMES WITH MAX PSD VALUES
% ========

    if ~isempty(fitTimes)
        
        t0 = fitTimes(1);  
        fitTimesCentered = fitTimes - t0;

        p = polyfit(fitTimesCentered, fitFreqs, 1);
        fitLine = polyval(p, fitTimesCentered);  

        slope = p(1); % slope is in Hz/day here
        %slope_hz_per_sec = slope * (24 * 60 * 60);  % convert to Hz/s
        slope_hz_per_sec = slope * (1 / 86400); % 1 day = 24x60x60 = 86400s

    else
        disp('No data points in selected time/frequency range.');
    end

    fprintf('Fit times span: %.2f days (%.2f hours)\n', ...
        max(fitTimesCentered) - min(fitTimesCentered), ...
        (max(fitTimesCentered) - min(fitTimesCentered)) * 24);

% ========
% ~~~ FIGURE 1: SPECTROGRAM WITH PSD, SELECTED DATA POINTS FOR STORM DATES
% ========
    figure(1)
    %set(G, 'Units', 'inches');         % Set units to inches
    %set(G, 'Position', [1, 1, 10, 4]);   % [x, y, width, height] in inches
    clf

    % Subplot 1: Spectrogram with PSD (rel to 1 Pa^2/Hz) 
    h=subplot(4,1,1:2);
    imagesc(Td,F,10*log10(PSDsm+eps));
    %imagesc(Td,F,10*log10(PSDdiffsm+eps));
    xlim(xl)
    %xlim(xl_full)
    %datetick('x',3)
    set(gca,'xtick',[])
    axis xy
    %bookfonts
    ylabel('F (Hz)')
    colormap(jet)
    ax=get(h,'Position');
    s=mean(std(10*log10(PSDsm)));
    m=median(median(10*log10(PSDsm)));
    %caxis([-10 15])
    %caxis([-35 5])
    %caxis([-160 -40])
    %caxis([m-3*s m+4*s])
    %caxis(PSDlim)
    hc=colorbar;
    %set(get(hc,'label'),'string','PSD (dB rel. 1 Pa^2/Hz)');
    set(get(hc,'label'),'string','PSD (dB rel. 1 (m/s)^2 /Hz)');
    set(h,'Position',ax);
    %axis tight
    %ylim([0 0.12])
    %ylim([0.01 0.12])
    ylim([0.03 0.12])
    %ylim([0 0.5])
    %caxis([-35 15])
    title('DR01 storm spectrogram, selected max PSD data points for L2-norm fit'); %,'FontSize', 16);
   
    tick_locations=[datenum(2015,1,1:365),datenum(2016,1,1:365)];
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'XTick',tick_locations)
    datetick('x',6,'keeplimits', 'keepticks')
    xtickangle(0)
    
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

    % Subplot 2: Linear trend (fitFreqs vs fitTimes)
    subplot(4,1,3:4)
    hold on
    plot(specTimes, peakFreqs, '.', 'Color', 'black'); % All peaks in black
    plot(fitTimes, fitFreqs, 'ro', 'LineWidth', 2);   % Linear fit slope in red
    %datetick('x',6,'keeplimits','keepticks')
    %xtickangle(90)
    %ylim([0 0.12])
    ylim([0.03 0.12])
    xlim(xl)
    ylabel('Freq (Hz)')
    xlabel('Date')

    tick_locations=[datenum(2015,1,1:365),datenum(2016,1,1:365)];
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca,'XTick',tick_locations)
    datetick('x',6,'keeplimits', 'keepticks')
    xtickangle(0)
    %title(sprintf('Linear fit: slope = %.2e Hz/day', slope))
    %bookfonts
end


%% LINEAR MODEL SLOPE CALCULATION
% Figure 2, Slope
% Figure 3, histogram of residuals (Hz)

%~~~~~~~
% CALCULATE SLOPE GIVEN TIMES AND PSD
% using fitlm, an L2-NORM FIT
%~~~~~~~~

% === GLOBAL FONT SETTINGS ===
set(groot, 'DefaultAxesFontName', 'Helvetica');
set(groot, 'DefaultTextFontName', 'Helvetica');
set(groot, 'DefaultAxesFontSize', 18);
set(groot, 'DefaultAxesFontWeight', 'bold');
set(groot, 'DefaultAxesLineWidth', 2);
set(groot, 'DefaultAxesTickLength', [0.01 0.01]);
set(groot, 'DefaultAxesTickDir', 'in');  % or 'out' if preferred

% Center time to the first time point (limits eval to slope day)
t0 = fitTimes(1);
fitTimesCentered = fitTimes - t0;

tbl = table(fitTimesCentered, fitFreqs);
lm = fitlm(tbl, 'fitFreqs ~ fitTimesCentered'); % fitlm is an L2 fit
% fitFreqs is a 2479x1 double
% fitTimesCentered is a 2479x1 double
% slope of the line is delta f over delta t, which is 4.237E-07Hz/s

% === Regression statistics ===
Rsq = lm.Rsquared.Ordinary;     % R-squared
RMSE = lm.RMSE;                 % Root Mean Square Error
anovaTbl = anova(lm, 'summary');
Fval = anovaTbl.F(1);           % F-statistic
pval = anovaTbl.pValue(1);      % p-value
residuals = lm.Residuals.Raw;
std_resid = std(residuals);  % Standard deviation of residuals

%The fitlm function in MATLAB (from the Statistics and Machine Learning Toolbox) 
%performs an ordinary least squares (OLS) regression, which is an L2 norm method. 
%It minimizes the sum of squared residuals.

%L1-norm (least absolute deviations) minimizes and is more robust to outliers.
%L2-norm (least squares) minimizes and is more sensitive to outliers, 
% but gives the classic "best linear unbiased estimate" under Gaussian noise assumptions.

% Extract slope and convert
slope = lm.Coefficients.Estimate(2);  
slope_hz_per_sec_2 = slope / (24 * 60 * 60);  % Convert from Hz/day to Hz/s

% Predict values and get confidence intervals (95% default)
[y_pred, y_ci] = predict(lm, fitTimesCentered);

% Print output
fprintf('Slope (Hz/s) using fitlm: %.3e\n', slope_hz_per_sec_2);

fprintf('Fit times span: %.2f days (%.2f hours)\n', ...
        max(fitTimesCentered) - min(fitTimesCentered), ...
        (max(fitTimesCentered) - min(fitTimesCentered)) * 24);

fprintf('Slope (Hz/s) using polyfit: %.3e\n', slope_hz_per_sec); 
%fprintf('Pt 2 Slope (Hz/s): %.3e\n', slope_hz_per_sec_2);
fprintf('R-squared: %.4f\n', Rsq);
fprintf('F-statistic: %.2f (p = %.4g)\n', Fval, pval);
fprintf('RMSE: %.3e Hz\n', RMSE);

% ~~~ Plot slope with residuals shaded II
figure(2);
hold on;

% Plot data points
h1 = plot(fitTimes, fitFreqs, 'o', ...
          'MarkerEdgeColor', 'black', ...
          'MarkerFaceColor', 'black', ...
          'MarkerSize', 5); 

% Fitted line
y_fit = predict(lm, fitTimesCentered);
h2 = plot(fitTimes, y_fit, 'r-', 'LineWidth', 3); 

% Format and label
datetick('x', 'keeplimits');
xlabel('Time');
ylabel('Frequency (Hz)');
ylim([0.04 0.08]);

startStr = datestr(tStart, 'mm/dd');
endStr = datestr(tEnd, 'mm/dd');
yearStr = datestr(tStart, 'yyyy');
title({'L2-norm (least squares) Linear Fit', ...
       ['of Max PSD from ', startStr, '–', endStr, ' ', yearStr]});
set(gca, 'FontSize', 18);

% Correct legend with all handles
legend([h1, h2], ...
       {'Max PSD', ...
        sprintf('\\Deltaf/\\Deltat = %.3e Hz/s', slope_hz_per_sec_2)}, ...
       'Location', 'northwest');

% Create histogram of residuals
figure(3);
histogram(residuals, 30, 'FaceColor', [0.2 0.2 0.8], 'EdgeColor', 'black');
xlabel('Residual (Hz)');
ylabel('Count');
title('Histogram of Residuals');
grid on;
set(gca, 'FontSize', 14);


%% Distance to foci, travel time, and error calculations

% %Equation 1 and equation 2 formulas and variable outputs to save in table below

% Cathles et al., 2009
% Eqn 1: the distance to focus x of waves arriving at a particular
% observation site may be written: 
% x = (g/4pi)*[(df/dt)^-1]

% The L2-fit model provides a RMSE in terms of frequency (+/- on y-axis)
% over the selected duration of time on the x-axis (need to conv fitTimes to (s)). 
% this RMSE (Hz) is divided by the selected time duration (s) and used to
% find the distance error calculation by this equation: 

% [(df/dt)^2]*4pi*[RMSE in Hz / time duration (s)]

% ~~~~~for example, storm 1 2015


pi = 3.141592653589793; 
g = 9.81; 

s = slope_hz_per_sec_2; % this is slope, df/dt (Hz/s) from lmfit L2-fit above "slope_hz_per_sec_2"
s2 = s*s;

R_H = RMSE; % this is the RMSE (Hz) from L2-fit model above "RMSE" 

% this is the total time over which slope is calculated, 
delta_t_d = max(fitTimesCentered) - min(fitTimesCentered); %days
delta_t_s = (max(fitTimesCentered) - min(fitTimesCentered)) * (60 * 60 * 24); % sec

distance = (g/(4*pi))*((s)^-1); 

% Error in meters = [(df/dt)^2]*4pi*[RMSE in Hz / time duration (s)]
error_m = (g/((s^2)*4*pi))*(R_H/delta_t_s); % error in meters

error_dist_max = error_m + distance; % error_m +/- to distance
error_dist_min = distance - error_m; % error_m +/- to distance

error_perc = (error_m/distance)*100; 

% Eqn 2: the total time of transit deltat(f) of waves at frequency
% (selected, a given f) required to travel the distance x is given by: 
% delta t = [4pi(f_sel))/g]*x

range = fMax - fMin; % 0.05-0.03 = 0.02
f_int = range/3; % 0.02/3 = 0.006666666666667
f_int0 = fMin; 
f_int1 = fMin+f_int; 
f_int2 = f_int1+f_int; 
f_int3 = f_int2+f_int; %should equal fMax

% use equation two where f~ is equal to iterations of f_int0 through f_int3
% each of the outputs will be a travel time
% the f_int0 output of travel time will be the low end of the time
% f_int3 will be high end of time window where storm origin can be searched
% input is in x dist in meters, output is time in s
% convert s to min to hrs to days

travel_t0 = [(4*pi*f_int0)/g]*distance; % f_int0 associated with f_min;
travel_t1 = [(4*pi*f_int1)/g]*distance;
travel_t2 = [(4*pi*f_int2)/g]*distance;
travel_t3= [(4*pi*f_int3)/g]*distance; % f_int3 assc w f_max; 

travel_t0_sec = travel_t0;
travel_t0_min = (travel_t0_sec / 60); % seconds / 60 = min
travel_t0_hrs= travel_t0_min / 60; % min / 60 = hours
travel_t0_days = travel_t0_hrs / 24; % hrs / 24 = days

travel_t3_sec = travel_t3;
travel_t3_min = (travel_t3_sec / 60); % seconds / 60 = min
travel_t3_hrs= travel_t3_min / 60; % min / 60 = hours
travel_t3_days = travel_t3_hrs / 24; % hrs / 24 = days

%% Table summarizing slope and regression results ===

% Define output folder
outputFolder = '/Volumes/CATALOGDR01/2025FIGS/StormCatalog/Stats/';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Convert start/end datetime to formatted strings for filename
tStart_str = datestr(tStart, 'yyyy-mm-dd_HH-MM-SS');
tEnd_str   = datestr(tEnd,   'yyyy-mm-dd_HH-MM-SS');
tStart_HMS = datestr(tStart, 'HH:MM:SS');
tEnd_HMS   = datestr(tEnd,   'HH:MM:SS');

% Number of data points used
nPoints = height(tbl);

% % Force slope to display as full decimal string (up to 15 digits precision)
slope_decimal = str2double(sprintf('%.25f', slope_hz_per_sec_2));
error_m_str = str2double(sprintf('%.25f', error_m)); % error in meters
error_dist_max_str= str2double(sprintf('%.25f', error_dist_max)); % error_m +/- to distance
error_dist_min_str= str2double(sprintf('%.25f',error_dist_min)); % error_m +/- to distance

% Create result table with split numeric frequency columns
resultTbl = table(...
    fMin, ...
    fMax, ...
    {datestr(tStart, 'yyyy-mm-dd')}, ...
    {tStart_HMS}, ...
    {datestr(tEnd, 'yyyy-mm-dd')}, ...
    {tEnd_HMS}, ...
    slope_decimal, ...
    RMSE, ...
    Rsq, ...
    nPoints, ...
    travel_t0_sec, ...
    travel_t0_days, ...
    travel_t3_sec, ...
    travel_t3_days, ...
    distance, ...
    error_m_str, ...
    error_dist_max_str, ...
    error_dist_min_str, ...
    error_perc, ...
    'VariableNames', {...
        'fMin', ...
        'fMax', ...
        'StartDate', ...
        'StartTime', ...
        'EndDate', ...
        'EndTime', ...
        'Slope_Hz_per_s', ...
        'RMSE_Hz', ...
        'R_squared', ...
        'NumPointsUsed', ...
        'min_travel_time_s', ...
        'min_travel_time_d', ...
        'max_travel_time_s', ...
        'max_travel_time_d', ...
        'distance', ...
        'error_m', ...
        'error_dist_max', ...
        'error_dist_min', ...
        'error_perc'});

% Compose file path and save
csv_filename = sprintf('SlopeSummary2_%s_to_%s.csv', tStart_str, tEnd_str);
csv_path = fullfile(outputFolder, csv_filename);
writetable(resultTbl, csv_path);

fprintf('Results saved to: %s\n', csv_path);

elapsedTime = toc; % End timing
stopTime = datetime('now');
disp(['Script ended at: ', datestr(stopTime)])
disp(['Total elapsed time: ', num2str(elapsedTime, '%.2f'), ' seconds'])


%% END
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% removed plot code for slope figure

% create a table of all inputs to save as csv
% frequency range set, dates in hr-min-sec start and end, slope, rmse, r
% squared values, and total number of sampled values used in linear model. 
% save the csv by the tstart and tend values to a specified folder.

% === Create table summarizing slope and regression results ===

% % Convert start/end datetime to formatted strings (hr-min-sec)
% tStart_str = datestr(tStart, 'yyyy-mm-dd_HH-MM-SS');
% tEnd_str   = datestr(tEnd,   'yyyy-mm-dd_HH-MM-SS');
% 
% % Convert datetime to separate time-of-day strings for table
% tStart_HMS = datestr(tStart, 'HH:MM:SS');
% tEnd_HMS   = datestr(tEnd,   'HH:MM:SS');
% 
% % Number of data points used
% nPoints = height(tbl);  % or numel(fitTimes)
% 
% % Define frequency bounds
% fMin_str = sprintf('%.3f', fMin);
% fMax_str = sprintf('%.3f', fMax);
% fRange = [fMin_str, '–', fMax_str, ' Hz'];
% % FrequencyRange_Hz results in a cell: 0.045‚Äì0.080 Hz
% % instead, create two columns, noe for fMin and one for Mfax, values only
% 
% % Create result table
% resultTbl = table(...
%     {fRange}, ...
%     {datestr(tStart, 'yyyy-mm-dd')}, ...
%     {tStart_HMS}, ...
%     {datestr(tEnd, 'yyyy-mm-dd')}, ...
%     {tEnd_HMS}, ...
%     slope_hz_per_sec_2, ... %Slope_Hz_per_s change format of 4.23678118627055e-07 to all decimal
%     RMSE, ...
%     Rsq, ...
%     nPoints, ...
%     'VariableNames', {...
%         'FrequencyRange_Hz', ...
%         'StartDate', ...
%         'StartTime', ...
%         'EndDate', ...
%         'EndTime', ...
%         'Slope_Hz_per_s', ...
%         'RMSE_Hz', ...
%         'R_squared', ...
%         'NumPointsUsed'});
% 
% 
% % Define output folder
% outputFolder = '/Volumes/CATALOGDR01/2025FIGS/StormCatalog/Stats/';
% 
% csv_filename = sprintf('SlopeSummary_%s_to_%s.csv', tStart_str, tEnd_str);
% csv_path = fullfile(outputFolder, csv_filename);
% 
% % Save table as CSV
% writetable(resultTbl, csv_path);
% 
% fprintf('Results saved to: %s\n', csv_path);
    % % Acceleration Time vector conversion from SAC header
    % [month,day]=monthday(sachdr.nzyear,sachdr.nzjday);
    % tzmin=datenum(sachdr.nzyear,month,day,sachdr.nzhour,sachdr.nzmin,sachdr.sec);
    % tzmax=tzmin+length(adata)*sachdr.delta/(24*60*60);
    % Td=linspace(tzmin,tzmax,length(adata));

    % to bandpass acceleration data to swell-band / did not use    
    % Td_s=Td;
    % N=numel(adata);
    % if N/2 ~= round(N/2)
    %     N=N-1;
    %     adata=adata(1:end-1);
    %     Td_s=Td(1:end-1);
    % end
    % 
    % f_spec = Fs*(0:(N/2))/N;
    % nyquist = 0.5;
    % 
    % %Swell band
    % [b,a]=butter(2,[0.03/nyquist, 0.12/nyquist],'bandpass');
    % data_swell_spec=fft(filtfilt(b,a,adata));

%~~~~~~~ Plot slope 
% figure(2);
% plot(fitTimes, fitFreqs, 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black', 'MarkerSize', 5); 
% %plot(fitTimes, fitFreqs, '.', 'Color', 'black');
% hold on;
% plot(fitTimes, predict(lm, fitTimesCentered), 'r-','LineWidth', 2);  % Fitted line
% datetick('x', 'keeplimits');
% xlabel('Time'); %, 'FontSize', 16);
% ylabel('Frequency (Hz)'); %, 'FontSize', 16);
% %ylim([0 0.12])
% ylim([0.04 0.075])
% startStr = datestr(tStart, 'mm/dd');
% endStr = datestr(tEnd, 'mm/dd');
% yearStr = datestr(tStart, 'yyyy');
% title({'L2-norm (least squares) Linear Fit', ...
%        ['of Max PSD from ', startStr, '–', endStr, ' ', yearStr]});
% set(gca, 'FontSize', 18);
% legend('Max PSD', ...
%        sprintf('\\Deltaf/\\Deltat = %.3e Hz/s', slope_hz_per_sec), ...
%        'Location', 'northwest'); %, 'FontSize', 16);

% ~~~~~~~~Plot slope with residuals shaded

% figure(2);
% plot(fitTimes, fitFreqs, 'o', ...
%      'MarkerEdgeColor', 'black', ...
%      'MarkerFaceColor', 'black', ...
%      'MarkerSize', 5); 
% hold on;
% 
% % Fitted line
% y_fit = predict(lm, fitTimesCentered);
% plot(fitTimes, y_fit, 'r-', 'LineWidth', 2); 
% 
% % Residual shaded region (±1 std dev)
% y_upper = y_fit + std_resid;
% y_lower = y_fit - std_resid;
% x_vals = [fitTimes; flipud(fitTimes)];
% y_vals = [y_upper; flipud(y_lower)];
% fill(x_vals, y_vals, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% % Create shaded area between upper and lower bounds
% x_vals = [fitTimes; flipud(fitTimes)];
% y_vals = [y_upper; flipud(y_lower)];
% 
% h_fill = fill(x_vals, y_vals, 'r', ...
%               'FaceAlpha', 0.04, ...
%               'EdgeColor', 'none');  % Capture handle for legend
% 
% % Format and label
% datetick('x', 'keeplimits');
% xlabel('Time');
% ylabel('Frequency (Hz)');
% ylim([0.04 0.075]);
% 
% title({'L2-norm (least squares) Linear Fit', ...
%        ['of Max PSD from ', startStr, '–', endStr, ' ', yearStr]});
% set(gca, 'FontSize', 18);
% 
% legend('Max PSD', ...
%        sprintf('\\Deltaf/\\Deltat = %.3e Hz/s', slope_hz_per_sec), ...
%        'Location', 'northwest');

% ~~~~ plot slope with varying residuals blue highlight
% blue shaded 95% confidence band that varies with fitTimes

% figure(2);
% hold on;
% 
% % Plot data points
% h1 = plot(fitTimes, fitFreqs, 'o', ...
%           'MarkerEdgeColor', 'black', ...
%           'MarkerFaceColor', 'black', ...
%           'MarkerSize', 5); 
% 
% % Fitted line
% y_fit = predict(lm, fitTimesCentered);
% h2 = plot(fitTimes, y_fit, 'r-', 'LineWidth', 3); 
% 
% % Residual shaded region (±1 std dev, red, more transparent)
% y_upper_std = y_fit + std_resid;
% y_lower_std = y_fit - std_resid;
% x_vals = [fitTimes; flipud(fitTimes)];
% y_vals_std = [y_upper_std; flipud(y_lower_std)];
% h_fill_std = fill(x_vals, y_vals_std, 'r', ...
%               'FaceAlpha', 0.08, ...
%               'EdgeColor', 'none');
% 
% % 95% Confidence Interval band (blue, varying width)
% % y_ci was obtained earlier by [y_pred, y_ci] = predict(lm, fitTimesCentered);
% x_vals = [fitTimes; flipud(fitTimes)];
% y_vals_ci = [y_ci(:,1); flipud(y_ci(:,2))];
% h_fill_ci = fill(x_vals, y_vals_ci, 'b', ...
%                  'FaceAlpha', 0.15, ...
%                  'EdgeColor', 'none');
% 
% % Bring plot elements to front
% uistack(h2, 'top');
% uistack(h1, 'top');
% 
% % Format and label
% datetick('x', 'keeplimits');
% xlabel('Time');
% ylabel('Frequency (Hz)');
% ylim([0.04 0.08]);
% 
% title({'L2-norm (least squares) Linear Fit', ...
%        ['of Max PSD from ', startStr, '–', endStr, ' ', yearStr]});
% set(gca, 'FontSize', 18);
% 
% % Legend: Add handles for the confidence band
% legend([h1, h2, h_fill_std, h_fill_ci], ...
%        {'Max PSD', ...
%         sprintf('\\Deltaf/\\Deltat = %.3e Hz/s', slope_hz_per_sec), ...
%         '±1σ Residual Band', ...
%         '95% Confidence Band'}, ...
%        'Location', 'northwest');



















%% demonstrate L1 vs L2 fit

% === FIGURE 4: Compare L1 (LAR) and L2 Fits ===

% % Convert time to days from t0 (for numeric stability)
% fitTimesDays = fitTimesCentered * 24;
% 
% % L1-norm fit (using 'LAR' for least absolute residuals)
% [fitL1Obj, gofL1] = fit(fitTimesDays, fitFreqs, 'poly1', 'Robust', 'LAR');

% % for figure(4), try this update to fix polyfit badly conditioned error
%[fitL1Obj, gofL1] = fit(fitTimesDays, fitFreqs, 'poly1', 'Robust', 'LAR');

% % L2-norm fit using polyfit
% fitL2Obj = polyfit(fitTimesDays, fitFreqs, 1);  % L2 slope/intercept
% 
% % Smooth x-axis vector in days
% x_plot_days = linspace(min(fitTimesDays), max(fitTimesDays), 100);
% 
% % Predict y-values
% yL2 = polyval(fitL2Obj, x_plot_days);
% yL1 = feval(fitL1Obj, x_plot_days);
% 
% % Convert x values back to datenum for plotting
% x_plot = t0 + x_plot_days / 24;
% 
% % Plot
% figure(4); clf;
% plot(fitTimes, fitFreqs, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); hold on;
% plot(x_plot, yL2, 'b-', 'LineWidth', 2);  % L2 line
% plot(x_plot, yL1, 'r--', 'LineWidth', 2); % L1 line
% datetick('x', 'keeplimits');
% xlabel('Time', 'FontSize', 16);
% ylabel('Peak Frequency (Hz)', 'FontSize', 16);
% title('L1 vs L2 Fit of Peak Frequency (Storm 1)', 'FontSize', 16);
% legend('Data', 'L2: Least Squares (OLS)', 'L1: Least Absolute Deviation (LAR)', ...
%        'Location', 'northwest', 'FontSize', 14);
% set(gca, 'FontSize', 16);
% box on


%%
% Notes on lm in Toolbox:
    % b_0 is the intercept (constant term).
    % b_1 is the slope (coefficient for the predictor variable x).
    % for ex: lm = fitlm(x, y);  % Fits a linear model y = b0 + b1*x
    % we're fitting fitTimesCentered = b1*fitFreqs + b0 (in y = mx+ b format)
    % if you type in: lm.Coefficients
    % This will display a table with the following columns:
    % Name: The name of each coefficient (e.g., "Intercept", "x" for a single predictor).
    % Estimate: The estimated value for each coefficient.
    % SE: Standard error of the coefficient estimate.
    % tStat: t-statistic for the hypothesis test (Estimate / SE).
    % pValue: The p-value for testing if the coefficient is significantly different from zero.
    % If you access lm.Coefficients.Estimate, it will give you a vector of these values. 
    % For example:
    % intercept = lm.Coefficients.Estimate(1);  % Intercept (b_0)
    % slope = lm.Coefficients.Estimate(2);      % Slope (b_1)

%% 
% figure(2);
% plot(fitTimes, fitFreqs, 'co');        % Original data
% hold on;
% plot(fitTimes, fitLine, 'b-', 'LineWidth', 2);  % Fitted line
% datetick('x', 'keeplimits');           % Format time axis
% xlabel('Time');
% ylabel('Peak Frequency (Hz)');
% title('Linear Fit of Peak Frequency vs Time');
% %legend('Data', 'Linear Fit', upperright); % change to upper right

    % %~~ new start
    %     % MEAN, SMOOTHING
    % 
    % medPSD=mean(PSD,2);
    % 
    % [~,C]=size(PSD); %build a matrix for pre-event spectral normalization
    % 
    % PSDnorm=repmat(medPSD,1,C);
    % 
    % % Smooth the PSD spectrograms a bit (quasi-Gaussian) Smoothed PSD, PSDsm
    % PSDsm=imgaussfilt(PSD,[3,11]);
    % %PSDsm=imgaussfilt(PSD,[1.5,11]); 
    % %PSDsm=imgaussfilt(PSD,[0.5,5]); 
    % %PSDsm = imgaussfilt(10*log10(PSD ./ PSDnorm), [3, 11]);  % dB before smoothing
    % 
    % % Normalized PSD (rel med), PSDdiff
    % PSDdiff=PSD./PSDnorm;
    % 
    % % Smoothed normalized PSD (rel med), PSDdiffsm
    % PSDdiffsm=PSDdiff;
    % %PSDdiffsm=imgaussfilt(PSDdiff,[0.5,5]);
    % %PSDdiffsm=imgaussfilt(PSDdiff,[1.5,11]);
    % PSDdiffsm=imgaussfilt(PSDdiff,[3,11]);
    % %~~ new end
