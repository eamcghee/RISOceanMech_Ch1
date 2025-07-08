clear all
close all
clc

% download images from url with coords
% url: https://worldview.earthdata.nasa.gov/?v=-2842718.1915111137,-2897383.0069994843,3167436.8267968767,55397.30320002348&p=antarctic&l=Reference_Labels_15m(hidden),Reference_Features_15m(hidden),Coastlines_15m,AMSRU2_Sea_Ice_Concentration_12km(min=24),VIIRS_NOAA21_CorrectedReflectance_TrueColor(hidden),VIIRS_NOAA20_CorrectedReflectance_TrueColor(hidden),VIIRS_SNPP_CorrectedReflectance_TrueColor(hidden),MODIS_Aqua_CorrectedReflectance_TrueColor(hidden),MODIS_Terra_CorrectedReflectance_TrueColor(hidden)&lg=false&t=2015-02-25-T15%3A34%3A23Z
% coords: -77.0801, 106.4898 (top right); -57.8866, -146.8847 (bottom left)
% jpg, 1km.  

% Define image folder and get list of files
%image_folder = '/Users/em/PROJECTS/SeaIce/IMGS/2014_Nov1toJan1_v2';
%image_folder = '/Users/em/PROJECTS/SeaIce/IMGS/2015_Jan1toApr1_v2';
%image_folder = '/Users/em/PROJECTS/SeaIce/IMGS/2015_Apr1toJul1_v2';
%image_folder = '/Users/em/PROJECTS/SeaIce/IMGS/2015_Jul1toOct1_v2';
%image_folder = '/Users/em/PROJECTS/SeaIce/IMGS/2015_Oct1toJan1_v2';
%image_folder = '/Users/em/PROJECTS/SeaIce/IMGS/2016_Jan1toApr1_v2';
%image_folder = '/Users/em/PROJECTS/SeaIce/IMGS/2016_Apr1toJul1_v2';
%image_folder = '/Users/em/PROJECTS/SeaIce/IMGS/2016_Jul1toOct1_v2';
image_folder = '/Users/em/PROJECTS/SeaIce/IMGS/2016_Oct1toJan1_v2';
image_files = dir(fullfile(image_folder, '*.jpg'));

%% interactively select the trapezoid to mask other-than Channel
% go all the way up to 69.3678 deg S

% % Load one image to define the trapezoid on
% img = imread('/Users/em/PROJECTS/SeaIce/IMGS/2015_Jan1toApr1_v2/snapshot-2015-03-01T00_00_00Z.jpg');  % Replace with your image path
% imshow(img);
% title('Click 4 corners of the trapezoid (clockwise or counterclockwise)');
% 
% % Use ginput to get 4 (x,y) points from mouse clicks
% [x_trap, y_trap] = ginput(4);
% 
% % Draw the trapezoid on the image
% hold on;
% plot([x_trap; x_trap(1)], [y_trap; y_trap(1)], 'r-', 'LineWidth', 2);
% title('Trapezoid defined');

%% manually define trapezoid (or run above and input x & y trap coords)

% snapshot-2015-02-25T00_00_00Z.png
% Define image folder and get list of files
%image_folder = '/Users/em/PROJECTS/SeaIce/IMGS/2014_Nov1toJan1_v2';
image_files = dir(fullfile(image_folder, '*.jpg'));  % Or '*.jpg'

%x_trap = [100, 400, 450, 50];  % X coordinates of corners 
%y_trap = [50, 50, 300, 300];   % Y coordinates of corners

x_trap = [1559.45550715619, 1723.64219041693, 849.032358431861, 640.641568139390];%East Corridor v2
y_trap = [798.542003733665, 893.265090230242, 1928.90416925949, 1543.69695084007];%East Corridor v2

% Define trapezoid coordinates (in pixels)
%x_trap = [554.329001367989, 95.5383036935704, 254.465800273598, 674.274281805746]; %East Corridor 
%y_trap = [818.376880984952, 1241.18399452804, 1433.09644322845, 914.333105335157]; %East Corridor

img = imread(fullfile(image_folder, image_files(1).name));  % change index as needed
imshow(img);
hold on;
plot([x_trap, x_trap(1)], [y_trap, y_trap(1)], 'r-', 'LineWidth', 2);
title('Trapezoid Overlay');
%% Given input image coordinates, use x_trap and y_trap values to determine the Degree Decimal coordinates of the trapezoid. 

% -83.8731, 126.3745 (top right); -64.5801, -160.1703 (bottom left)

% take coordinates and create a shapefile to import into overview figures
% so that the corridors match exactly. 

%%
% Loop over images
for k = 1:length(image_files)
    % Read image
    image_path = fullfile(image_folder, image_files(k).name);
    img = imread(image_path);

    % Convert to grayscale if it's RGB
    if size(img, 3) == 3
        img_gray = rgb2gray(img);
    else
        img_gray = img;
    end

    % Create a binary mask using poly2mask
    [rows, cols] = size(img_gray);
    mask = poly2mask(x_trap, y_trap, rows, cols);

    % Apply mask to grayscale image
    region_pixels = img_gray(mask);

    % Count non-black and black pixels
    non_black_count = sum(region_pixels > 0);
    black_count = sum(region_pixels == 0);
    total_count = non_black_count + black_count;
    ratio_non_black = non_black_count / total_count;
    perc_non_black = ratio_non_black * 100;

    % Extract datetime string from filename
    % Example: snapshot-2015-02-25T00_00_00Z.png
    filename = image_files(k).name;
    tokens = regexp(filename, 'snapshot-(\d{4}-\d{2}-\d{2})T\d{2}_\d{2}_\d{2}Z', 'tokens');
    if ~isempty(tokens)
        date_str = tokens{1}{1};  % '2015-02-25'
        date_num = datetime(date_str, 'InputFormat', 'yyyy-MM-dd');
    else
        date_num = NaT;  % Not a time, if the format doesn't match
    end

    % Display result
    fprintf('%s: Non-black pixels in trapezoid = %d\n', filename, non_black_count);
    fprintf('%s: Percent sea ice in east corridor (by pixels): %.2f%%\n', filename, perc_non_black);

    % Store results
    results(k).filename = filename;
    results(k).date = date_num;
    results(k).nonblack = non_black_count;
    results(k).black_count = black_count;
    results(k).total_count = total_count;
    results(k).ratio_non_black = ratio_non_black;
    results(k).perc_non_black = perc_non_black;
end

% Convert results to table
T = struct2table(results);

% Save results to CSV
writetable(T, '/Users/em/PROJECTS/SeaIce/SeaIceCorridor_Apr25/EastCorr_nonblack_pixel_counts_multday3_2016_Oct1toJan1_v2.csv');

%% add Time_UTC column to csv

% Read the CSV file
T = readtable('/Users/em/PROJECTS/SeaIce/SeaIceCorridor_Apr25/EastCorr_nonblack_pixel_counts_multday3_2016_Oct1toJan1_v2.csv');

% Extract the datetime string using regular expressions
expr = 'snapshot-(\d{4}-\d{2}-\d{2}T\d{2}_\d{2}_\d{2}Z)';
tokens = regexp(T.filename, expr, 'tokens', 'once');

% Convert the tokens from cell arrays to character vectors
timestampStrs = cellfun(@(x) x{1}, tokens, 'UniformOutput', false);

% Replace underscores with colons to match ISO 8601 format
timestampStrs = strrep(timestampStrs, '_', ':');

% Store in a new column
T.Time_UTC = timestampStrs;

% (Optional) Save to a new CSV file
writetable(T, '/Users/em/PROJECTS/SeaIce/SeaIceCorridor_Apr25/EastCorr_nonblack_pixel_counts_multday_2016_Oct1toJan1_new_v2.csv');