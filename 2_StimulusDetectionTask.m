% STIMULUS DETECTION TASK:
% The task involves a monkey fixating on a visual cue and responding to stimuli randomly distributed throughout the visual feild.

% TargetLocations:          A two-column matrix where each row corresponds to a unique target location (represented by x and y coordinates). 
% allCondID:                For each trial, an index indicating the row in 'TargetLocations' that corresponds to the target location used.
% allBehaviorOutcomes:      For each trial, a binary value indicating if the trial was successful (1) or not (0).
% allBehaviorOutcomesSTR:   For each trial, a string describing the specific outcome.
% allFixationTime:          For each trial, the duration of fixation prior to the stimulus.
% allResponseTime:          For each trial, two times indicating when in ms after the stimulus the eye left the fixation window (1st row) and when it entered the target window (2nd row).

clear all; clc; close all

% Directory of data files 
SaveDir = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/StimulusWEyetrace/Figures';
%FileDir = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/StimulusWEyetrace/Stimulus';

FileDir = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/StimulusWEyetrace/AllStimulus';

cd(SaveDir);

% Get a list of all files in the directory
fileList = dir(fullfile(FileDir, '*.mat'));
numFiles = numel(fileList);

% Initialize variables to store the concatenated data
concatenated_allBehaviorOutcomes = [];
concatenated_allBehaviorOutcomesSTR = [];
concatenated_allCondID = [];
concatenated_allFixationTime = [];
concatenated_allResponseTime = [];
concatenated_TargetLocations = [];

% Iterate over the files
for i = 1:numFiles
    % Get the file name and full path
    fileName = fileList(i).name;
    filePath = fullfile(FileDir, fileName);
    
    % Load the data from the file
    data = load(filePath);
    
    % Concatenate the arrays
    concatenated_allBehaviorOutcomes = [concatenated_allBehaviorOutcomes, data.allBehaviorOutcomes];
    concatenated_allBehaviorOutcomesSTR = [concatenated_allBehaviorOutcomesSTR, data.allBehaviorOutcomesSTR];
    concatenated_allCondID = [concatenated_allCondID, data.allCondID];
    concatenated_allFixationTime = [concatenated_allFixationTime, data.allFixationTime];
    concatenated_allResponseTime = [concatenated_allResponseTime, data.allResponseTime];
    concatenated_TargetLocations = [concatenated_TargetLocations, data.TargetLocations];
end

% Specify the filenames and the names of the arrays you want to concatenate
filenames = {fileList.name};
array_names = {'allBehaviorOutcomes', 'allBehaviorOutcomesSTR', 'allCondID', 'allFixationTime', 'allResponseTime', 'TargetLocations'};



% Fixing variable names

% Removing concatinated from variable names=
allBehaviorOutcomes = concatenated_allBehaviorOutcomes;
allBehaviorOutcomesSTR = concatenated_allBehaviorOutcomesSTR;
allCondID = concatenated_allCondID;
allFixationTime = concatenated_allFixationTime;
allResponseTime = concatenated_allResponseTime;

% Adding TargetLocations
TargetLocations = data.TargetLocations;

% Clearing unnecessary variables
% clearvars -except allBehaviorOutcomes allBehaviorOutcomesSTR allCondID allFixationTime allResponseTime SaveDir TargetLocations;



%% %%%%%%%%%%%%%%%%%% %%
%%% Basic Statistics %%%
%%%% %%%%%%%%%%%%%% %%%%

%Percent Correct trials
sum(allBehaviorOutcomes)/numel(allBehaviorOutcomes)




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram of trial outcomes

% Find the unique Trial outcomes
unique_outcomes = unique(allBehaviorOutcomesSTR)

% Count the occurrences of each unique string
counts = histcounts(categorical(allBehaviorOutcomesSTR), categorical(unique_outcomes))

% Create a bar plot to visualize the counts
bar(counts)
xticklabels(unique_outcomes)
xlabel('Trial Outcomes')
ylabel('Counts')


%% histogram of Fixation Times
histogram(allFixationTime,'BinWidth', 10, 'EdgeColor', 'black')

xlabel('Time (ms)');
ylabel('Number of Trials');
title('Fixation Time Distribution');


%% Histogram of all Response times
histogram(allResponseTime,'BinWidth', 1, 'EdgeColor', 'black')

xlabel('Time (ms)');
ylabel('Number of Trials');
title('Response Time Distribution');

%% Histogram of coordinates for each trial
histogram(allCondID,'BinWidth', 1, 'EdgeColor', 'black')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Heat Maps  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Black and white Heat Map of Trial Counts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_targets = length(TargetLocations);
n_trials = length(allResponseTime);
c_n_trials = zeros(n_targets, 1);  %How many Correct of no Responce trials at each target


%Calculation the responce time for every coordinate
for i = 1:n_trials;
    cond_id = allCondID(i);
    if strcmp(allBehaviorOutcomesSTR(i), 'CorrectResponse') == 1 || strcmp(allBehaviorOutcomesSTR(i), 'NoResponse') == 1 || strcmp(allBehaviorOutcomesSTR(i), 'FlyingOverTime') == 1|| strcmp(allBehaviorOutcomesSTR(i), 'WrongChoice') == 1|| strcmp(allBehaviorOutcomesSTR(i), 'EyeCrossRespWindow') == 1 ; %
        c_n_trials(cond_id) = c_n_trials(cond_id) + 1;
    else;
    end;
end;


% Reshape trial count to a grid
x = unique(TargetLocations(:, 1));
y = unique(TargetLocations(:, 2));
[X, Y] = meshgrid(x, y);
Z = nan(length(y), length(x));
for i = 1:n_targets;
    ix = find(x == TargetLocations(i, 1));
    iy = find(y == TargetLocations(i, 2));
    Z(iy, ix) = c_n_trials(i);
end

% Flip the orientation of the y-axis
Z = flipud(Z);

%Define Grid of trial count 
cnGrid = Z






% Create a black and white colormap
colormap(gray);

% Create a heatmap of the trial count values
m = max(Z(:));
h = imagesc(x, y, Z, [0 m]); % specify the x, y, Z, and the color range

% Add this line to make each cell square
axis equal;

colorbar('Ticks', [0 m*.2 m*.4 m*.6 m*.8 m]);
set(get(colorbar, 'Label'), 'String', 'Trial Count');
xlabel('Visual Degrees');
ylabel('Visual Degrees');
title('Trial counts CorrectResponse and NoResponse Trials');
yticklabels(flipud(get(gca, 'yticklabels')));

% Set NaN values to be transparent
set(h, 'alphadata', ~isnan(Z));



% Save the figure as SVG in the specified folder with the name '1'
saveas(gcf, '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Figure_Images/StimulusDetection/1.svg', 'svg');


%% Trial count basic stats
nanmean(cnGrid(:))
min(cnGrid(:))
max(cnGrid(:))
nanstd(cnGrid(:))

%% input coordinate to trouble shoot  

% Input the target coordinate
targetX = 0
targetY = -1

% Find the index of the target location in TargetLocations
targetIndex = find(TargetLocations(:, 1) == targetX & TargetLocations(:, 2) == targetY);

% Validate if targetIndex is not empty
if isempty(targetIndex)
    fprintf('No trials found for the target location (%d, %d)\n', targetX, targetY);
else
    % Find the trials that used this target location
    targetTrials = find(allCondID == targetIndex);
    
    % Get the count of trials and allBehaviorOutcomesSTR values for these trials
    trialCount = length(targetTrials);
    trialOutcomes = allBehaviorOutcomesSTR(targetTrials);
    
    % Display the count and allBehaviorOutcomesSTR values
    fprintf('Count of trials for location (%d, %d): %d\n', targetX, targetY, trialCount);
    fprintf('Behavior outcomes for these trials are:\n');
    disp(trialOutcomes);
end





%% Heat map of accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Count the number of correct and incorrect trials for each coordinate
n_targets = length(TargetLocations);
n_trials = length(allBehaviorOutcomes);
n_correct = zeros(n_targets, 1);
n_incorrect = zeros(n_targets, 1);
for i = 1:n_trials
    cond_id = allCondID(i);
    if strcmp(allBehaviorOutcomesSTR(i), 'CorrectResponse') == 1
        n_correct(cond_id) = n_correct(cond_id) + 1;
    elseif strcmp(allBehaviorOutcomesSTR(i), 'NoResponse') == 1
        n_incorrect(cond_id) = n_incorrect(cond_id) + 1;
    end
end


% Compute the accuracy for each coordinate
accuracy = n_correct ./ (n_correct + n_incorrect);

% Reshape accuracy to a grid
x = unique(TargetLocations(:, 1));
y = unique(TargetLocations(:, 2));
[X, Y] = meshgrid(x, y);
Z = nan(length(y), length(x));
for i = 1:n_targets
    ix = find(x == TargetLocations(i, 1));
    iy = find(y == TargetLocations(i, 2));
    Z(iy, ix) = accuracy(i);
end

% Flip the orientation of the y-axis
Z = flipud(Z);

%Define Grid of Accuracy 
accGrid = Z;

% non-smoothed heatmap of the accuracy values
% Create a colormap with an extra color for NaN values
figure;
map = [gray; jet(64)];
colormap(map);




% Create a heatmap of the accuracy values

h = imagesc(x, y, Z, [0 1]); % specify the x, y, Z and the color range
colormap(jet);


% Add this line to make each cell square
axis equal;

cbar = colorbar('Ticks', [0 0.2 0.4 0.6 0.8 1]);
set(get(cbar, 'Label'), 'String', 'Probability');
xlabel('Visual Degrees');
ylabel('Visual Degrees');
title('Accuracy CorrectResponse vs NoResponse');
% Set the y-axis tick labels
yticklabels(flipud(get(gca, 'yticklabels')));

% Set NaN values to be transparent
set(h, 'alphadata', ~isnan(Z));

% Remove NaN values from the colorbar
%set(cbar, 'ylim', [0 1]);




% ... (previous plotting code)

% Save the figure as SVG in the specified folder with the name '1'
saveas(gcf, '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Figure_Images/StimulusDetection/2.svg', 'svg');



%% Accuracy basic stats
nanmean(accGrid(:))
min(accGrid(:))
max(accGrid(:))
nanstd(accGrid(:))


%% Heat Map of mean RT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_targets = length(TargetLocations);
n_trials = length(allResponseTime);
r_trials = zeros(n_targets, 1);
c_n_trials = zeros(n_targets, 1);
allResponseTime2 = allResponseTime(2,:);
%allResponseTime2 = allResponseTime(1,:);  % when fixation breaks instead of 


%Calculation the responce time for every coordinate
for i = 1:n_trials;
    cond_id = allCondID(i);
    if strcmp(allBehaviorOutcomesSTR(i), 'CorrectResponse') == 1 ;
    %if strcmp(allBehaviorOutcomesSTR(i), 'NoResponse') == 1;

        r_trials(cond_id)= r_trials(cond_id) + allResponseTime2(i);
        c_n_trials(cond_id) = c_n_trials(cond_id) + 1;
    else
    end
end


% Compute the response for each coordinate
aveResponse = r_trials ./ c_n_trials;


% Reshape aveResponse to a grid
x = unique(TargetLocations(:, 1));
y = unique(TargetLocations(:, 2));
[X, Y] = meshgrid(x, y);
Z = nan(length(y), length(x));
for i = 1:n_targets;
    ix = find(x == TargetLocations(i, 1));
    iy = find(y == TargetLocations(i, 2));
    Z(iy, ix) = aveResponse(i);
end

% Flip the orientation of the y-axis
Z = flipud(Z);

%Define Grid of Accuracy 
aveGrid = Z


%non-smoothed heatmap of the CorrentResponse and Noresponse times

% Create a colormap with an extra color for NaN values
figure;
map = [gray; jet(64)];
colormap(map);



% Create a heatmap of the RT values

h = imagesc(x, y, Z, [0 500]); % specify the x, y, Z and the color range
colormap(jet);

% Add this line to make each cell square
axis equal;

cbar = colorbar('Ticks', [0 100 200 300 400 500]);
set(get(cbar, 'Label'), 'String', 'ms');
xlabel('Visual Degrees');
ylabel('Visual Degrees');
title('Response Time of CorrectResponse and NoResponse Trials');
% Flip the y-axis tick labels
yticklabels(flipud(get(gca, 'yticklabels')));


% Set NaN values to be transparent
set(h, 'alphadata', ~isnan(Z));

% Remove NaN values from the colorbar
set(cbar, 'ylim', [0 500]);

% ... (previous plotting code)

% Save the figure as SVG in the specified folder with the name '1'
saveas(gcf, '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Figure_Images/StimulusDetection/3.svg', 'svg');

%% RT basic stats
nanmean(aveGrid(:))
nanstd(aveGrid(:))


%% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  ANOVA S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard deviation of mean RT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% not necessary and somehow the same as the accuracy SDs in the folowing section

% Initialize variables
n_targets = size(TargetLocations, 1); % number of targets
ssd = zeros(n_targets, 1); % initialize the sum of squared differences to 0 for each target
n_trials = zeros(n_targets, 1); % initialize the number of trials to 0 for each target
allResponseTime2 = allResponseTime(2,:); %Extract only 2nd row


% Calculation the response time for every coordinate
for i = 1:length(allCondID) % loop through all trials
    cond_id = allCondID(i); % get the condition ID for the current trial
    if strcmp(allBehaviorOutcomesSTR(i), 'CorrectResponse') || strcmp(allBehaviorOutcomesSTR(i), 'NoResponse') % check if the trial was a correct response or no response
        n_trials(cond_id) = n_trials(cond_id) + 1; % increment the number of trials for the current target
        diff = allResponseTime2(i) - aveResponse(cond_id); % calculate the difference between the current response time and the average response time for the current target
        ssd(cond_id) = ssd(cond_id) + diff^2; % add the squared difference to the sum of squared differences for the current target
    end
end

% Compute the standard deviation for each coordinate
sdResponse = sqrt(ssd ./ (n_trials - 1)); % divide the sum of squared differences by (n_trials - 1) and take the square root to get the standard deviation for each target



% Reshape aveResponse to a grid
x = unique(TargetLocations(:, 1));
y = unique(TargetLocations(:, 2));
[X, Y] = meshgrid(x, y);
Z = nan(length(y), length(x));
for i = 1:n_targets;
    ix = find(x == TargetLocations(i, 1));
    iy = find(y == TargetLocations(i, 2));
    Z(iy, ix) = sdResponse(i);
end

% Flip the orientation of the y-axis
Z = flipud(Z);

%Define Grid of Accuracy 
RTsdGrid = Z



%non-smoothed heatmap of the CorrentResponse and Noresponse DSs

% Create a colormap with an extra color for NaN values
figure;
map = [gray; jet(64)];
colormap(map);

% Create a heatmap of the accuracy values

h = imagesc(x, y, Z, [0 500]); % specify the x, y, Z and the color range
colormap(jet);
cbar = colorbar('Ticks', [0 100 200 300 400 500]);
set(get(cbar, 'Label'), 'String', 'ms');
xlabel('Visual Degrees');
ylabel('Visual Degrees');
title('SD of CorrectResponse and NoResponse Trials');
% Flip the y-axis tick labels
yticklabels(flipud(get(gca, 'yticklabels')));


% Set NaN values to be transparent
set(h, 'alphadata', ~isnan(Z));

% Remove NaN values from the colorbar
set(cbar, 'ylim', [0 500]);


grid on

saveas(gcf, 'SD of CorrectResponse and NoResponse Trials smooth', 'jpeg');


%% Standard deviation accuracy values

% not necessary and probably not done correctly

% Initialize variables
n_targets = size(TargetLocations, 1); % number of targets
ssd = zeros(n_targets, 1); % initialize the sum of squared differences to 0 for each target
n_trials = zeros(n_targets, 1); % initialize the number of trials to 0 for each target
allResponseTime2 = allResponseTime(2,:); %Extract only 2nd row


% Calculation the response time for every coordinate
for i = 1:length(allCondID) % loop through all trials
    cond_id = allCondID(i); % get the condition ID for the current trial
    if strcmp(allBehaviorOutcomesSTR(i), 'CorrectResponse') || strcmp(allBehaviorOutcomesSTR(i), 'NoResponse') % check if the trial was a correct response or no response
        n_trials(cond_id) = n_trials(cond_id) + 1; % increment the number of trials for the current target
        diff = allResponseTime2(i) - aveResponse(cond_id); % calculate the difference between the current response time and the average response time for the current target
        ssd(cond_id) = ssd(cond_id) + diff^2; % add the squared difference to the sum of squared differences for the current target
    end
end

% Compute the standard deviation for each coordinate
sdResponse = sqrt(ssd ./ (n_trials - 1)); % divide the sum of squared differences by (n_trials - 1) and take the square root to get the standard deviation for each target



% Reshape aveResponse to a grid
x = unique(TargetLocations(:, 1));
y = unique(TargetLocations(:, 2));
[X, Y] = meshgrid(x, y);
Z = nan(length(y), length(x));
for i = 1:n_targets;
    ix = find(x == TargetLocations(i, 1));
    iy = find(y == TargetLocations(i, 2));
    Z(iy, ix) = sdResponse(i);
end

% Flip the orientation of the y-axis
Z = flipud(Z);

%Define Grid of Accuracy 
ACCsdGrid = Z



%non-smoothed heatmap of the CorrentResponse and Noresponse DSs

% Create a colormap with an extra color for NaN values
figure;
map = [gray; jet(64)];
colormap(map);

% Create a heatmap of the accuracy values

h = imagesc(x, y, Z, [0 500]); % specify the x, y, Z and the color range
colormap(jet);
cbar = colorbar('Ticks', [0 100 200 300 400 500]);
set(get(cbar, 'Label'), 'String', 'ms');
xlabel('Visual Degrees');
ylabel('Visual Degrees');
title('SD of CorrectResponse and NoResponse Trials');
% Flip the y-axis tick labels
yticklabels(flipud(get(gca, 'yticklabels')));


% Set NaN values to be transparent
set(h, 'alphadata', ~isnan(Z));

% Remove NaN values from the colorbar
set(cbar, 'ylim', [0 500]);


grid on

saveas(gcf, 'SD of CorrectResponse and NoResponse Trials smooth', 'jpeg');



%% Define groups

%Define observed group


% accGrid = percent of trials accurate per coordinate
% aveGrid = mean RT 
% cnGrid = number of trials
% ACCsdGrid = standard deviation of accuracy
% RTsdGrid = standard deviation of mean RT


% get the dimensions of the matrix
[nrows, ncols] = size(accGrid)

% calculate the index of the middle row/column
midrow = floor(nrows/2) + 1
midcol = floor(ncols/2) + 1


% Defining Control and Experimental group data


%defining ACCURACY groups
expAcc = accGrid(1:midrow-1,midcol+1:end)
conAcc = [accGrid(1:midrow-1,1:midcol-1);accGrid(midrow+1:end,midcol+1:end);accGrid(midrow+1:end,1:midcol-1)]
%Make 1d
expAcc = expAcc(:)
conAcc = conAcc(:)

%defining RT groups
expAve = aveGrid(1:midrow-1,midcol+1:end)
conAve = [aveGrid(1:midrow-1,1:midcol-1);aveGrid(midrow+1:end,midcol+1:end);aveGrid(midrow+1:end,1:midcol-1)]
%Make 1d
expAve = expAve(:)
conAve = conAve(:)

%defining count groups
expCn = cnGrid(1:midrow-1,midcol+1:end)
conCn = [cnGrid(1:midrow-1,1:midcol-1);cnGrid(midrow+1:end,midcol+1:end);cnGrid(midrow+1:end,1:midcol-1)]
%Make 1d
expCn = expCn(:)
conCn = conCn(:)

%defining accuracy SD groups
expSd = ACCsdGrid(1:midrow-1,midcol+1:end)
conSd = [ACCsdGrid(1:midrow-1,1:midcol-1);ACCsdGrid(midrow+1:end,midcol+1:end);ACCsdGrid(midrow+1:end,1:midcol-1)]
%Make 1d
expSd = expSd(:)
conSd = conSd(:)

%defining RT SD groups
expSd = RTsdGrid(1:midrow-1,midcol+1:end)
conSd = [RTsdGrid(1:midrow-1,1:midcol-1);RTsdGrid(midrow+1:end,midcol+1:end);RTsdGrid(midrow+1:end,1:midcol-1)]
%Make 1d
expSd = expSd(:)
conSd = conSd(:)


%% basic stats for exp and control group
%accuracy
mean(expAcc)
mean(conAcc)
std(expAcc)
std(conAcc)

mean(expAcc)-mean(conAcc)



% RT
mean(expAve)
mean(conAve)
std(expAve)
std(conAve)

mean(expAve)-mean(conAve)




%% Weighted means RTs

% Calculate the weighted means and standard deviations for the experimental group
n_exp = sum(expCn);
weights_exp = expCn / n_exp;
exp_mean = sum(weights_exp .* expAve);
exp_stdev = sqrt(sum(weights_exp .* (expAve - exp_mean).^2) ./ sum(weights_exp .* (expCn - 1) ./ expCn));

% Calculate the weighted means and standard deviations for the control group
n_con = sum(conCn);
weights_con = conCn / n_con;
con_mean = sum(weights_con .* conAve);
con_stdev = sqrt(sum(weights_con .* (conAve - con_mean).^2) ./ sum(weights_con .* (conCn - 1) ./ conCn));



%% Weighted Anova of RTs

% Concatenate the experimental and control data into a single variable
data = [expAve; conAve];

% Create a grouping variable indicating which observations belong to the experimental and control groups
group = [repmat({'Experimental'}, size(expAve)); repmat({'Control'}, size(conAve))];

% Perform a one-way ANOVA with the 'group' variable as the grouping variable
[p,tbl,stats] = anova1(data, group, 'varnames', {'Group', 'Measure'});

% Print the ANOVA table
disp(tbl);

% If the p-value is less than the chosen significance level (e.g., 0.05), reject the null hypothesis
if p < 0.05
    disp('There is a significant difference in RT between the experimental and control groups');
else
    disp('There is not a significant difference in RT between the experimental and control groups');
end


%% RT assumptions testing


% Assuming expAve and conAve are column vectors
rtData = [expAve; conAve];  % Combine experimental and control group data
groups = [ones(size(expAve)); 2 * ones(size(conAve))];  % Create a group variable (1 for experimental, 2 for control)

% Levene's Test for Homogeneity of Variances
% MATLAB does not have a built-in function for Levene's test, so you would typically use a custom function or a statistical toolbox that includes it.
% However, you can use the 'vartestn' function for a similar purpose, which is Bartlett's test by default but can perform Levene's test by setting 'TestType' to 'LeveneAbsolute'.

% Perform Levene's test using vartestn without extracting all outputs
[pLevene, tblLevene] = vartestn(rtData, groups, 'TestType', 'LeveneAbsolute', 'Display', 'off');

% Display the p-value from Levene's test
disp(['Levene''s Test p-value for RT data: ', num2str(pLevene)]);

%% Normality
% Kolmogorov-Smirnov Test for Normality
% You would typically perform this test on the residuals of a model, but you can also test your raw data.
% Here's how you might do it on the combined data:
[~, pKS, ksstat] = kstest(rtData - mean(rtData));

% Display the p-value from the KS test
disp(['Kolmogorov-Smirnov Test p-value for RT data: ', num2str(pKS)]);
disp(['Kolmogorov-Smirnov Test test statistic: ', num2str(ksstat)]);

%% weighted means Accuracy

% Calculate the weighted means and standard deviations for the experimental group's accuracy
n_exp = sum(expCn);
weights_exp = expCn / n_exp;
exp_mean = sum(weights_exp .* expAcc);  % Replaced expAve with expAcc for accuracy
exp_stdev = sqrt(sum(weights_exp .* (expAcc - exp_mean).^2) ./ sum(weights_exp .* (expCn - 1) ./ expCn));

% Calculate the weighted means and standard deviations for the control group's accuracy
n_con = sum(conCn);
weights_con = conCn / n_con;
con_mean = sum(weights_con .* conAcc);  % Replaced conAve with conAcc for accuracy
con_stdev = sqrt(sum(weights_con .* (conAcc - con_mean).^2) ./ sum(weights_con .* (conCn - 1) ./ conCn));


%% Weighted Anova Accuracy

data = [expAcc; conAcc];

% Create a grouping variable indicating which observations belong to the experimental and control groups
group = [repmat({'Experimental'}, size(expAcc)); repmat({'Control'}, size(conAcc))];

% Perform a one-way ANOVA with the 'group' variable as the grouping variable
[p,tbl,stats] = anova1(data, group, 'varnames', {'Group', 'Measure'});

% Print the ANOVA table
disp(tbl);

% If the p-value is less than the chosen significance level (e.g., 0.05), reject the null hypothesis
if p < 0.05
    disp('There is a significant difference in accuracy between the experimental and control groups');
else
    disp('There is not a significant difference in accuracy between the experimental and control groups');
end



%% Non-Parametric tests
% 
% 
%% RT Simple wilcoxon

% Perform the Mann-Whitney U test (Wilcoxon rank-sum test) for two groups
[p, h, stats] = ranksum(expAve, conAve);

% Display the test results
disp(['Mann-Whitney U test p-value for RT data: ', num2str(p)]);
if h == 1
    disp('There is a significant difference in RT between the experimental and control groups');
else
    disp('There is not a significant difference in RT between the experimental and control groups');
end

% Print the table
disp(stats);


%% Accuracy simple wilcoxon

% Perform the Mann-Whitney U test (Wilcoxon rank-sum test) for two groups
[p_mw, h_mw, stats_mw] = ranksum(expAcc, conAcc);

% Display the test results
disp(['Mann-Whitney U test p-value for RT data: ', num2str(p_mw)]);
if h_mw == 1
    disp('There is a significant difference in Accuracy between the experimental and control groups');
else
    disp('There is not a significant difference in Accuracy between the experimental and control groups');
end

% Print the ANOVA table
disp(stats_mw);


%% linear model

% Prepare data for weighted regression
% Concatenate your experimental and control data and their corresponding counts
rtData = [expAve; conAve];  % Reaction times
weights = [expCn; conCn];   % Corresponding weights for each observation

% Create a binary variable for group membership
% Assuming '1' for experimental and '0' for control
groupVariable = [ones(size(expAve)); zeros(size(conAve))];

% Perform weighted linear regression
% Note: The LinearModel.fit function allows for specification of weights
lm = fitlm(groupVariable, rtData, 'Weights', weights);

% Display regression results
disp(lm);


%% assumptions

% Perform weighted linear regression
lm = fitlm(groupVariable, rtData, 'Weights', weights);

% Check for linearity and homoscedasticity visually by plotting residuals vs. fitted values
figure;
plot(lm.Fitted, lm.Residuals.Raw, 'bo');
xlabel('Fitted values');
ylabel('Residuals');
title('Residuals vs Fitted Plot');
grid on;

% Add a reference line at 0
hold on;
refline(0);

% Check for normality of residuals with a Q-Q plot
figure;
qqplot(lm.Residuals.Raw);
title('Q-Q plot of the Residuals');

% If you have multiple predictors and want to check for multicollinearity, 
% you can look at the variance inflation factors (VIF)
if size(lm.DesignMatrix, 2) > 2  % Check if you have multiple predictors
    disp('Variance Inflation Factors (VIF):');
    disp(vif(lm));
end

% You can also use a formal test for homoscedasticity like the Breusch-Pagan test
% However, MATLAB does not have a built-in function for this test, 
% so you might need to implement it or use a custom function from the File Exchange.

% For normality, you can perform a statistical test like the Shapiro-Wilk test, 
% but MATLAB does not have a built-in function for this test either.
% As an alternative, you can use the Lilliefors test as a proxy, which is a 
% modification of the Kolmogorov-Smirnov test for normality.
[h, pValue] = lillietest(lm.Residuals.Raw);

% Display the result
disp(['Lilliefors test for normality, H = ', num2str(h), ', p-value = ', num2str(pValue)]);






%% RESULTS SECTION
%
% Statistical Analysis of Accuracy and Reaction Time:
%
% Two separate one-way analyses of variance (ANOVAs) were conducted to assess the differences
% in accuracy and reaction time (RT) between the upper right quadrant and the other quadrants in 
% the response space. For these analyses, responses from the upper right quadrant were assigned 
% to the experimental group, and responses from the other three quadrants were combined and assigned 
% to the control group.
%
% The first ANOVA focused on accuracy. The values along the x and y axes were excluded from the analysis 
% to avoid bias. The null hypothesis for this analysis posited no significant difference in accuracy 
% between the upper right quadrant and the other quadrants. The results provided significant evidence 
% to reject this null hypothesis, indicating a significant difference in accuracy between the two groups 
% (p < 0.05).
%
% The second ANOVA examined reaction time (RT). Like the accuracy analysis, RT values along the x and y axes 
% were excluded. The null hypothesis stated there would be no significant difference in RT between the upper 
% right quadrant and the other quadrants. Prior to the analysis, the weighted means of the reaction times were 
% computed for both the experimental and control groups. The ANOVA on RT yielded significant results, leading to 
% the rejection of the null hypothesis (p < 0.05).
%
% Collectively, these analyses suggest significant differences in both accuracy and reaction time between the 
% upper right quadrant and the other quadrants of the response space.
%}

