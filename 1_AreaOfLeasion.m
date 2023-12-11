clear all; close all

% Define image paths
img_path1 = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin Screenshots/T1w_20181213_pos1.png'
img_path2 = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin Screenshots/T1w_20210609_pos1.png';
img_path3 = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin Screenshots/T1w_20230221_pos1.png';

% Load images
img1 = imread(img_path1);
img2 = imread(img_path2);
img3 = imread(img_path3);

% Store the images in a cell array
images = {img1, img2, img3};

years = [2018, 2021, 2023]


% Define the pixel size in millimeters
pixelSize = 0.1; % 0.1 mm per pixel


% Convert to grayscale if they're RGB
img_grays = cell(size(images));
for i = 1:length(images)
    if size(images{i}, 3) == 3
        img_grays{i} = rgb2gray(images{i});
    else
        img_grays{i} = images{i};
    end
end

%% Display the first image and draw the ROI
figure;
imshow(img_grays{1});
title('Draw the Region of Interest');
roi = drawfreehand();
mask = roi.createMask(); % Create a binary mask of the ROI

% re establish draw roi tool, execute in command window
%roi.Selected = true;


%% OR load previously saved mask

load('/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Figure_Images/Lesion/mask.mat')


%% Draw the width of the brain

% Display the image
figure;
imshow(img_grays{1});
title('Draw a horizontal line on the image');

% Use drawline to let user draw on the image
hLine = drawline('Color', 'r');

%% OR load previously saved width

load('/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Figure_Images/Lesion/width.mat')

%% Calculate Pixels per centimeter

% Get line position in pixels
linePosition = hLine.Position;

% Calculate the length of the line in pixels
lineLengthPixels = sqrt(sum((linePosition(2,:) - linePosition(1,:)).^2));

% Set the length of the line as 5 cm
knownLengthCm = 5;

% Calculate pixels per centimeter
pixelsPerCm = lineLengthPixels / knownLengthCm;

% Display the result
disp(['Pixels per centimeter: ' num2str(pixelsPerCm)]);

pixelsPerCm = 100
%% Make contrast adjusted images

% Manually note down the low and high values you've chosen for Image 1
contrastLowImage1 = 20;  % Replace with your chosen low value
contrastHighImage1 = 150; % Replace with your chosen high value

% Normalize the contrast limits to [0, 1] if needed for Image 1
if isinteger(img_grays{1})
    maxVal = double(intmax(class(img_grays{1})));
    contrastLowHighImage1 = [contrastLowImage1, contrastHighImage1] / maxVal;
elseif isa(img_grays{1}, 'double')
    % No normalization needed, values should already be in range [0,1]
end

% Apply the same contrast limits to all three images
adjusted_grays = cell(size(img_grays));
for idx = 1:3

    % Apply the same contrast values for Images 2 and 3
    adjusted_grays{idx} = imadjust(img_grays{idx}, contrastLowHighImage1, []);

end



%% Calculate area and visualize

thresholdValue = 0.4; % Adjust this value as needed

% Initialize an array to store the calculated areas
areas_cm2 = zeros(1, 3);

figure;

for idx = 1:3
    % Binarize the adjusted image based on the threshold
    binaryImg = ~imbinarize(adjusted_grays{idx}, thresholdValue);

    % Identify the regions within the ROI where the binary mask is true
    binaryImg = binaryImg & mask;  % Apply the mask to the binary image

    % Calculate the area by summing all the non-zero pixels and converting to square centimeters
    area_cm2 = sum(binaryImg(:)) / (pixelsPerCm^2);

    % Store the area in the 'areas_cm2' array
    areas_cm2(idx) = area_cm2;

    % Create an RGB version of the grayscale image for overlay purposes
    overlayImg = cat(3, adjusted_grays{idx}, adjusted_grays{idx}, adjusted_grays{idx});

    % Overlay the lesion area in red on the RGB image
    redRegion = binaryImg;  % This is already defined by combining mask & binaryImg
    overlayImg(:,:,1) = overlayImg(:,:,1) + uint8(redRegion * 255);  % Increase red channel intensity
    overlayImg(:,:,2) = overlayImg(:,:,2) - uint8(redRegion * 255);  % Decrease green channel intensity
    overlayImg(:,:,3) = overlayImg(:,:,3) - uint8(redRegion * 255);  % Decrease blue channel intensity

    % Display the original image to the left
    subplot(3, 2, 2 * idx - 1);
    imshow(img_grays{idx});
    title(['Year ' num2str(years(idx))]);

    % Display the image with the overlay and calculated area to the right
    subplot(3, 2, 2 * idx);
    imshow(overlayImg);
    title(['Lesion Cross Section Area: ' num2str(area_cm2, '%.2f') ' cm²']);
end

% Display the calculated areas in square centimeters
disp('Areas (in square centimeters) for each image:');
disp(areas_cm2);


%% bar

% Display the calculated areas in square cm
disp('Areas (in square centimeters) for each image:');
disp(areas_cm2);

% Plot the areas in square centimeters for each image
figure;
bar(1:3, areas_cm2);
xlabel('Scan Year');
ylabel('Area (cm)');
title('Area of Lesion Crossection Over Time');
xticklabels(years);
ylim([0 .5])

% Add text labels above the bars
for idx = 1:3
    text(idx, areas_cm2(idx), [num2str(areas_cm2(idx), '%.2f')], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

%% change

firstchange = (((areas_cm2(2))/(areas_cm2(1)))-1)*100

secondchange = (((areas_cm2(3))/(areas_cm2(2)))-1)*100


%% plot
% Display the calculated areas in square cemtimeters
disp('Areas (in square cemtimeters) for each image:');
disp(areas_cm2);

% Create a vector for the x-axis (years)
x = 1:3;

% Plot the areas in square cemtimeters as points with "o" markers
figure;
plot(x, areas_cm2, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Scan Year');
ylabel('Area (cm²)');
title('Area of Lesion Crossection Over Time');
xticks(x);
xticklabels(years);
% ylim([7500 8600])


% Add data labels above the "o" markers
for idx = 1:3
    text(x(idx), areas_cm2(idx), [num2str(areas_cm2(idx), '%.2f') ' cm²'], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end




%% Set the opasity on lesion area

roi_opacity = .5;
figure;
image_version = adjusted_grays %adjusted_grays OR img_grays

for idx = 1:3





    % Create a pure red image for the overlay
    redOverlay = cat(3, uint8(redRegion*255), zeros(size(binaryImg), 'uint8'), zeros(size(binaryImg), 'uint8'));

    RGBs = cat(3, image_version{idx}, image_version{idx}, image_version{idx});

    redChannel = RGBs(:,:,1);
    redChannel(redRegion) = uint8(double(redChannel(redRegion)) * (1-roi_opacity) + 255 * roi_opacity);

    greenChannel = RGBs(:,:,2);
    greenChannel(redRegion) = uint8(double(greenChannel(redRegion)) * (1-roi_opacity));

    blueChannel = RGBs(:,:,3);
    blueChannel(redRegion) = uint8(double(blueChannel(redRegion)) * (1-roi_opacity));

    overlayImg = cat(3, redChannel, greenChannel, blueChannel);

    % Display the blended image in the original figure
    figure(1);
    subplot(3, 1, idx);
    imshow(overlayImg);
    title(['Year ' num2str(years(idx)) ' - Lesion Cross Section Area: ' num2str(area_cm2, '%.2f') ' cm²']);
end

% Display the calculated areas in square cemtimeters
disp('Areas (in square cemtimeters) for each image:');
disp(areas_cm2);

%%
% Specify the desired filename and format
outputFilename = 'fuuuuuuuuuuuuuuuk';

% Set the desired resolution in DPI
dpiValue = 1200;

% Save the figure
print(gcf, '-dpng', ['-r', num2str(dpiValue)], outputFilename);

disp(['Figure saved to: ', outputFilename]);

% Specify the desired filename and format
outputFilename = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Figure_Images/Lesion/Lesion_Overlay.svg';

% Set the desired resolution in DPI
dpiValue = 1200; % This won't affect the SVG, but will be used for embedded raster images

% Save the figure
print(gcf, '-dsvg', ['-r', num2str(dpiValue)], outputFilename);

disp(['Figure saved to: ', outputFilename]);


%% Methods

% Methods:
% We conducted an analysis of brain scan images obtained from a monkey at
% three different time points: 2018, 2021, and 2023. The images were processed
% using MATLAB with the following methods:
%
% - Image Preprocessing: The images were loaded and converted to grayscale.
%   Contrast adjustment was applied to enhance visibility in Image 1. The chosen
%   contrast limits for Image 1 were a low value of 20 and a high value of 150.
%   These contrast limits were then applied to all three images.
%
% - Region of Interest (ROI) Selection: An ROI was manually selected using the
%   'drawfreehand' function in Image 1, and a binary mask of the ROI was created
%   to isolate the lesion area.
%
% - Area Calculation: A threshold value of 0.3 was applied to binarize the adjusted
%   images. For each image (2018, 2021, and 2023), the lesion area was calculated
%   in square cemtimeters by summing all non-zero pixels, considering each pixel
%   represented 0.1 cm². The calculated areas were stored for analysis.
%

%% Results and Interpretation:
% The analysis reveals percentage changes in the size of the lesion over the
% three time points (2018, 2021, and 2023). These percentage changes indicate
% the relative variations in the lesion's size:
%
% - From 2018 to 2021, the lesion's area increased by approximately 9.22%.
% - From 2021 to 2023, there was a slight decrease in the lesion's area by
%   approximately 0.27%.
%
% These changes, while noticeable, are relatively small and may be within the
% range of noise inherent to the fMRI scans. The minimal overall change and the
% small decrease between 2021 and 2023 suggest a relatively stable size and
% shape of the lesion over the observed time period. Further studies and
% evaluations may help clarify any significant trends or variations in the
% lesion's characteristics.
%


