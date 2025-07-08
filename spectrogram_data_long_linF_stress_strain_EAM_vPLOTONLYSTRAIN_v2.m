clear all
close all
clc
fclose('all');
%%
% =========================================================================
% SPECTROGRAMS WITH STRAIN (SWELL-BAND AND IG-BAND)
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

%% Load and plot merged sac file

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
% % === Load the merged SAC file ===
% sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_001to090_merged.sac'; %1 Jan 2015 to 31 Mar 2015
% %sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_091to181_merged.sac'; %1 Apr 2015 to 30 Jun 2015
% %sacfile = '/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_182to273_merged.sac'; %1 Jul 2015 to 30 Sep 2015 
% %sacfile ='/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_274to365_merged.sac'; %1 Oct 2015 to 31 Dec 2015 
% 
% %read-in sac file & assign data and header (there are no headers here)
% [data, hdr] = rdsac(sacfile);
% 
% % === Build time vector ===
% dt = 1; % sample interval in seconds
% n = length(data);
% t = (0:n-1) * dt; % time in seconds
% 
% % Convert seconds to hours
% %t_hours = t / 3600;
% %t_days = t_hours/24;
% 
% % Define start datetime based on file (Day 091 = April 1, 2015)
% start_datetime = datetime(2015, 1, 1) ;  % Jan 1
% %start_datetime = datetime(2015, 1, 1)  + days(90);  % April 1 = day 91
% %start_datetime = datetime(2015, 1, 1)  + days(181);  % Jul 1 = day 182
% %start_datetime = datetime(2015, 1, 1)  + days(273);  % Oct 1 = day 274
% time_vec = start_datetime + seconds(t); % Build datetime vector from seconds
% 
% % === Plot ===
% H = figure(3);
% plot(time_vec, data, 'k');
% xlabel('Time (UTC)');
% ylabel('Amplitude (counts)');
% title('Merged 1 Hz SAC Data: DR01 LHZ');
% grid off;
% ylim([-1e6 1e6]);
% %ylim([-1e7 1e7]);
% 
% % JAN-MAR
% xlim([datetime(2015,1,1), datetime(2015,3,31,23,59,59)]);
% xtick_dates = datetime(2015,1,1):days(15):datetime(2015,3,31);
% xtick_dates = [xtick_dates, datetime(2015,3,31)];
% file_base_name = sprintf('Counts_DR01_2015_1Janto31Mar');
% 
% % APR-JUN
% %xlim([datetime(2015,4,1), datetime(2015,6,30,23,59,59)]);
% %xtick_dates = datetime(2015,4,1):days(15):datetime(2015,6,15);
% %xtick_dates = [xtick_dates, datetime(2015,6,30)];
% %file_base_name = sprintf('Counts_DR01_2015_1Aprto30Jun');
% 
% % JUL-SEP
% % xlim([datetime(2015,7,1), datetime(2015,9,30,23,59,59)]);
% % xtick_dates = datetime(2015,7,1):days(15):datetime(2015,9,15);
% % xtick_dates = [xtick_dates, datetime(2015,9,30)];
% % file_base_name = sprintf('Counts_DR01_2015_1Julto30Sep');
% 
% % OCT-DEC
% %xlim([datetime(2015,10,1), datetime(2015,12,31,23,59,59)]);
% %xtick_dates = datetime(2015,10,1):days(15):datetime(2015,12,15);
% %xtick_dates = [xtick_dates, datetime(2015,12,31)];
% %file_base_name = sprintf('Counts_DR01_2015_1Decto31Dec');
% 
% xticks(xtick_dates);
% xticklabels(datestr(xtick_dates, 'dd-mmm'));
% 
% set(H, 'Units', 'inches');
% set(H, 'Position', [1, 1, 30, 4]);
% 
%     % === SAVE AFTER EACH FIGURE ===
%     folder_path = '/Volumes/CATALOGDR01/2025FIGS/'; 
%     file_ext = '.fig';
%     dpi = 600;
%     fullFilePath = fullfile(folder_path, [file_base_name file_ext]);
% 
%     % Save as PNG
%     pngFilePath = fullfile(folder_path, [file_base_name, '.png']);
%     %print(H, pngFilePath, '-dpng', ['-r' num2str(dpi)]);
% 
%     % Save as MATLAB FIG
%     figFilePath = fullfile(folder_path, [file_base_name, '.fig']);
%     %saveas(H, figFilePath);  % Alternatively: savefig(H, figFilePath);
% 
%     %is 001to090 only Start date: 01-Jan-2015, End date:   11-Feb-2015
%     %01:53:53?

%%

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

% --- 1. Set integration window 
date_int_start = datenum(2015,1,1); 
date_int_end   = datenum(2015,4,1);
year = 2015;
year = num2str(year);
file_base_name = sprintf('DR01_2015_1Janto1Apr_Strain_10000_MS');
%file_base_name = sprintf('DR01_2015_1Aprto30Jun_Strain_10000');

% Set date range and component
xl = [datenum('01-01-2015'), datenum('04-01-2015')]; 
%xl = [datenum('04-01-2014'), datenum('07-01-2015')];

% === Manually define header info ===
% station = 'DR01';
% year = 2015;
% month = 1;
% day = 1;
% hour = 0;
% minute = 0;
% second = 0;
%delta = 1;       % sample interval in seconds (1 Hz → 1 sec)
%Fs = 1 / delta;  % sampling frequency

comp = 'LHZ';

% Set the path to your folder containing many individual SAC files
data_folder = '/Volumes/CATALOGDR00/Longfiles';

% === SYSTEM COMMAND TO BUILD LIST OF FILES ===
% This line executes a system-level `ls` command to list SAC files matching a pattern,
% and writes that list to a file called 'list.all'.
% Pattern explanation: 'DR01-2015-1-xx_LHZ.sac' likely refers to one SAC file per day/hour/etc.
system(['/bin/ls ', data_folder, '/DR01-2015-1-xx_', comp, '.sac > list.all']);

% Set the path to your data folder
%data_folder = '/Volumes/CATALOGDR00/Longfiles';system(['/bin/ls ', data_folder, '/DR01-2015-1-xx_', comp, '.sac > list.all']);
%data_folder = '/Volumes/CATALOGDR00/Longfiles';system(['/bin/ls ', data_folder, '/DR01-2014-3-xx_', comp, '.sac > list.all']); %Oct-Dec
%data_folder = '/Volumes/CATALOGDR00/DATA/Longfiles2/DR01/DR01_LHZ_--_2015_001to090_merged.sac'

%At this point, list.all is a text file with full paths (or relative paths) to each SAC file matching the pattern. For example:

% === COUNT THE NUMBER OF FILES ===
% This line uses `cat` and `wc -l` to count how many lines (files) are in list.all.
!cat list.all | wc -l > list.num

% Open file lists
fid1 = fopen('list.num'); % File containing number of files
fid2 = fopen('list.all'); %File containing list of SAC file names

% === READ FILE NAMES INTO AN ARRAY ===
i=1;
while 1
           tline = fgetl(fid2);  % Read one line (i.e., one SAC file path)
           
            if ~ischar(tline), break, end  % End of file reached
            if i==1
            stnames=tline; % First line
            else
                stnames=char(stnames,tline); % Append additional lines
            end
            i=i+1;
        end
        fclose(fid2);

        % === OPEN FILE LISTS ===
fid1 = fopen('list.num');  % File containing number of files
fid2 = fopen('list.all');  % File containing list of SAC file names

    %Now stnames is a character array where each row is a SAC file path like:
    %stnames(1,:) = '/Volumes/CATALOGDR00/Longfiles/DR01-2015-1-01_LHZ.sac'
    %stnames(2,:) = '/Volumes/CATALOGDR00/Longfiles/DR01-2015-1-02_LHZ.sac'

% === READ NUMBER OF FILES FROM list.num ===
nstas = fscanf(fid1, '%d');  % E.g., 90 files
fclose(fid1);

% =========================================================================
% MAIN PROCESSING LOOP
% =========================================================================

for i=1:nstas
    H=figure(i);

    % Load SAC file data and header    
    [sachdr,data]=load_sac(stnames(i,:));

    % For data with no SAC header
    %[~, data] = load_sac(stnames(i,:));  % Ignore header

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

    % % Time vector conversion from SAC header
    [month,day]=monthday(sachdr.nzyear,sachdr.nzjday);
    tzmin=datenum(sachdr.nzyear,month,day,sachdr.nzhour,sachdr.nzmin,sachdr.sec);
    tzmax=tzmin+length(data)*sachdr.delta/(24*60*60);
    Td=linspace(tzmin,tzmax,length(data));

    % Compute time vector manually because no SAC header
    % tzmin = datenum(year, month, day, hour, minute, second);
    % tzmax = tzmin + (length(data) - 1) * delta / (24*60*60);  % end time
    % Td = linspace(tzmin, tzmax, length(data));

    % % Set start datetime manually
    % tzmin = datenum(2015, 1, 1, 0, 0, 0);
    % tzmax = datenum(2015, 3, 31, 23, 0, 0);
    % 
    % % Calculate delta from known duration
    % actual_duration_sec = (tzmax - tzmin) * 86400;
    % delta = actual_duration_sec / (length(data) - 1);
    % Fs = 1 / delta;
    % Td = tzmin + (0:length(data)-1) * delta / 86400;  % Create accurate time vector

% =========================================================================
% SPECTROGRAM COMPUTATION. Option 1: Chen 2019; Option 2: Lipovsky 2018
% =========================================================================

    % Spectrogram window length
    NW=2048;
    freqs=[];

    Ax=[110 130 1000 1200];
    set(H,'Position',Ax);
    [S,F,T,PSD]=spectrogram(data,hanning(NW),round(0.99*NW),freqs,Fs);
    %[S,F,T,PSD]=spectrogram(detrend(adata),1e4,[],[],1);
    F=F/sachdr.delta/D; % Normalize frequency axis
    %F=F/delta/D; %for no SAC header
    [~,C]=size(PSD); 

% =======
% Option 1: Chen 2019
% Convert spectrogram to strain consistent with Chen et al. (2019)

    %     p=2.5;
    %     fc=.017;
    %     A=1;
    %     T_GSV=A./(1+(F/fc).^(2*p));
    %     k=sqrt(F.^1.7/(.02^1.7))*1e-3;
    %     c=2*pi*F./k;
    %     %transfer function for displacement
    %     T_GSV_k2=T_GSV.*k.^2;
    %     %transfer function for acceleration
    %     %to strain
    %     T_GSV_k2_s=T_GSV_k2./((2*pi*F).^2);
    %     Tmat_s=repmat(T_GSV_k2_s,1,C)*1e6;
    %     %to strain rate
    %     T_GSV_k2_sr=T_GSV_k2./(2*pi*F);
    %     Tmat_sr=repmat(T_GSV_k2_sr,1,C)*1e6;
    %     %instrument lower corner corrections
    %     f_ins=1/120;
    %     T_ins=1-(1./(1+(F/f_ins).^2));
    %     %zero frequency fix
    %     T_ins(1)=inf;
    %     T_ins_mat=repmat(T_ins,1,C);
    %     PSD_s=PSD.*Tmat_s.^2./(T_ins_mat.^2);
    %     PSD_sr=PSD.*Tmat_sr.^2./(T_ins_mat.^2);

% =======
% Option 2: Lipovsky 2018
% Convert spectrogram to stress consistent with Lipovsky (2018)
    
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
    figure(1); clf;
    
    % Spectrogram: 
    h=subplot(1,1,1);
    imagesc(Td,F,10*log10(PSDsm+eps));
    xlim(xl)
    axis xy
    %bookfonts
    ylabel('Frequency (Hz)')
    colormap(jet)
    hc=colorbar;
    set(get(hc,'label'),'string','PSD (dB rel. 1 Pa^2/Hz)');
    
    ylim([0 0.15])
    caxis([50 100])
    tick_locations = xl(1):15:xl(2); % Tick every 15 days
    set(gca,'XTick',tick_locations)
    datetick('x','mmm dd','keeplimits','keepticks')
    xtickangle(0)
    xlabel('Date')
    title([sachdr.kstnm(1:4), ' Spectrogram (0–0.15 Hz)'])
    %title([station(1), 'Spectrogram (0-0.15 Hz)'])
    set(gcf,'units','inches','position',[1 1 30 2])  % 5x2 inches
  
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
% Figure 1: Stress, Strain, Spectrogram

% ~~
% Swell-band Stress-strain: 
    %figure(2); clf;
    H = figure(2); clf;
    
    subplot(1,1,1)
    
    % Swell-band stress (left y-axis)
    yyaxis left
    plot(Td_s, movmax(abs(data_stress_swell), 3600), 'k-', 'LineWidth', 2); 
    ylabel('\sigma_x (Pa)')
    %ylim(1e4*[-.2 2])
    ylim([-1000 10000])

    % Swell-band strain (right y-axis)
    % Note: Where stress data is divided by (Ep) Young's modulus (in Pascals)
    % to find strain (unitless)
    % movmax(..., 3600) applies a moving maximum over a 1-hour window
    % abs(stress) ensures the stress is treated in magnitude only (ignoring sign).

    yyaxis right
    %plot(Td_s, movmax(abs(data_stress_swell), 3600)/Ep, 'k-'); %, 'linewidth', 2);
    ylabel('\epsilon_x')
    %ylim(1e4/Ep*[-.2 2])
    ylim([-1000/Ep 10000/Ep])
    
    xlim(xl)
    xtick_locations = xl(1):15:xl(2);
    set(gca,'XTick',xtick_locations)
    datetick('x','mmm dd','keeplimits','keepticks')
    xtickangle(0)
    xlabel('Date')
    
    title([sachdr.kstnm(1:4), ' ', year, ' Vertical Velocity Stress-Strain, 0.03–0.12 Hz Bandpass'])
    %title([station(1), ' ', year, ' Vertical Velocity Stress-STrain, 0.03-0.12 Hz Bandpass']) % for no SAC header
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    set(gcf,'units','inches','position',[1 1 20 3])  % Also 5x2 inches

    % === SAVE AFTER EACH FIGURE ===
    folder_path = '/Volumes/CATALOGDR01/2025FIGS/'; 
    file_ext = '.png';
    file_ext = '.fig';
    dpi = 600;
    fullFilePath = fullfile(folder_path, [file_base_name file_ext]);

    %set(H, 'Units', 'inches');
    %set(H, 'Position', [1, 1, 30, 2]);  %[1, 1, 5, 2]

    % Save as PNG
    pngFilePath = fullfile(folder_path, [file_base_name, '.png']);
    %print(H, pngFilePath, '-dpng', ['-r' num2str(dpi)]);
    
    % Save as MATLAB FIG
    figFilePath = fullfile(folder_path, [file_base_name, '.fig']);
    %saveas(H, figFilePath);  % Alternatively: savefig(H, figFilePath);

end
%%
    % =========================================================================
    % INTEGRATE SWELL-BAND STRESS OVER SPECIFIED DATE RANGE
    % =========================================================================
    
    % --- 1. Set integration window (set at top of code)
     
    % --- 2. Select stress values in that date range
    ind = (Td_s >= date_int_start) & (Td_s <= date_int_end);
    t_sub = Td_s(ind);
    stress_sub = movmax(abs(data_stress_swell), 3600);  % Already smoothed
    stress_sub = stress_sub(ind);

    % --- 3. Integrate using trapz
    dt_days = mean(diff(t_sub));              % timestep in days
    dt_sec = dt_days * 86400;                 % convert to seconds
    stress_integral = trapz(stress_sub) * dt_sec;  % Units: Pa·s
    
    % --- 4. 
    strain_integral = stress_integral/Ep; 

    % --- 5. Display result
    fprintf('%s Integrated stress from %s to %s: %.2e Pa·s\n', ...
        sachdr.kstnm(1:4), datestr(date_int_start), datestr(date_int_end), stress_integral);

    %fprintf('%s Integrated stress from %s to %s: %.2e Pa·s\n', ...
    %    station(1), datestr(date_int_start), datestr(date_int_end), stress_integral);

    fprintf('%s Integrated strain from %s to %s: %.2e', ...
        sachdr.kstnm(1:4), datestr(date_int_start), datestr(date_int_end), strain_integral);

    %fprintf('%s Integrated strain from %s to %s: %.2e', ...
    %    station(1), datestr(date_int_start), datestr(date_int_end), strain_integral);


    % Note: Where stress data is divided by (Ep) Young's modulus (in Pascals)
    % to find strain (unitless)
    % movmax(..., 3600) applies a moving maximum over a 1-hour window
    % abs(stress) ensures the stress is treated in magnitude only (ignoring sign).

%%
% =========================================================================
% INTEGRATE SWELL-BAND STRESS OVER SPECIFIED DATE RANGE (ABOVE 1.2 kPa ONLY)
% =========================================================================

% --- 1. Set integration window 

% --- 2. Select stress values in that date range
ind = (Td_s >= date_int_start) & (Td_s <= date_int_end);
t_sub = Td_s(ind);
stress_sub = movmax(abs(data_stress_swell), 3600);  % Already smoothed
stress_sub = stress_sub(ind);

% --- 3. Apply threshold mask: include only values >= 1.2 kPa (1200 Pa)
threshold = 1200;  % in Pa
mask = stress_sub >= threshold;
t_sub_thresh = t_sub(mask);
stress_sub_thresh = stress_sub(mask);

% --- 4. Integrate using trapz (only over values >= 1.2 kPa)
dt_days = mean(diff(t_sub));          % Average timestep in days
dt_sec = dt_days * 86400;             % Convert to seconds
stress_integral = trapz(stress_sub_thresh) * dt_sec;  % Units: Pa·s

% --- 5. Calculate strain integral
strain_integral = stress_integral / Ep;

% --- 6. Display result
fprintf('%s Integrated stress (≥1.2 kPa) from %s to %s: %.2e Pa·s\n', ...
    sachdr.kstnm(1:4), datestr(date_int_start), datestr(date_int_end), stress_integral);

%fprintf('%s Integrated stress (≥1.2 kPa) from %s to %s: %.2e Pa·s\n', ...
%    station(1), datestr(date_int_start), datestr(date_int_end), stress_integral);

fprintf('%s Integrated strain (≥1.2 kPa) from %s to %s: %.2e\n', ...
    sachdr.kstnm(1:4), datestr(date_int_start), datestr(date_int_end), strain_integral);

%fprintf('%s Integrated strain (≥1.2 kPa) from %s to %s: %.2e\n', ...
%    station(1), datestr(date_int_start), datestr(date_int_end), strain_integral);
    
    %% NOTES

    % The movmax(...,3600) smooths over 1 hour of data; your stress is in Pascals.
    
    % The integration gives Pascals × seconds, which is consistent with stress applied over time.
    
    % Compute the average stress:
    %stress_mean = mean(stress_sub);
