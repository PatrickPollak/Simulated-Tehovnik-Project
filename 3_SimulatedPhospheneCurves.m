clear  all
close all
clc

% Directory of data files 
SaveDir = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Simulated_Phosphenes/Figures';
FileDir = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Simulated_Phosphenes/Stimulus_Files';
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
concatinated_allFlagID = [];

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
    concatinated_allFlagID = [concatinated_allFlagID, data.allFlagID];
end

% Specify the filenames and the names of the arrays you want to concatenate
filenames = {fileList.name};
array_names = {'allBehaviorOutcomes', 'allBehaviorOutcomesSTR', 'allCondID', 'allFixationTime', 'allResponseTime', 'TargetLocations', 'allFlagID'};



% Fixing variable names

% Removing concatinated from variable names=
allBehaviorOutcomes = concatenated_allBehaviorOutcomes;
allBehaviorOutcomesSTR = concatenated_allBehaviorOutcomesSTR;
allCondID = concatenated_allCondID;
allFixationTime = concatenated_allFixationTime;
allResponseTime = concatenated_allResponseTime;
allFlagID = concatinated_allFlagID;

% Adding TargetLocations
TargetLocations = data.TargetLocations;

% Clearing unnecessary variables
%clearvars -except allBehaviorOutcomes allBehaviorOutcomesSTR allCondID allFixationTime allResponseTime SaveDir TargetLocations;



% TargetLocations:          A two-column matrix where each row corresponds to a unique target location (represented by x and y coordinates). 
% allCondID:                For each trial, an index indicating the row in 'TargetLocations' that corresponds to the target location used.
% allBehaviorOutcomes:      For each trial, a binary value indicating if the trial was successful (1) or not (0).
% allBehaviorOutcomesSTR:   For each trial, a string describing the specific outcome.
% allFixationTime:          For each trial, the duration of fixation prior to the stimulus.
% allResponseTime:          For each trial, two times indicating when the eye left the fixation window (1st row) and when it entered the target window (2nd row).
% allFlagID:                For each trial, a flag indicating the type; 0 for catch trials, 1 for no phosphene saccade trials, and 2 for phosphene saccade trials.




%% %%%%%%%%%%%%%%%%%% %%
%%% Basic Statistics %%%
%%%% %%%%%%%%%%%%%% %%%%

%Percent Correct trials
sum(allBehaviorOutcomes)/numel(allBehaviorOutcomes)




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram of trial outcomes

% Find the unique Trial outcomes
unique_outcomes = unique(allBehaviorOutcomesSTR);

% Count the occurrences of each unique string
counts = histcounts(categorical(allBehaviorOutcomesSTR), categorical(unique_outcomes));

% Compute total trials excluding 'PreFixBreak' and 'PreStimBreak'
totalTrials = sum(counts(~ismember(unique_outcomes, {'PreFixBreak', 'PreStimBreak'})));

% Create a bar plot to visualize the counts
barHandle = bar(counts);
xticklabels(unique_outcomes);
xlabel('Trial Outcomes');
ylabel('Counts');

% Add counts and percentage on top of each bar
for i = 1:length(counts)
    if ~ismember(unique_outcomes{i}, {'PreFixBreak', 'PreStimBreak'})
        percentage = (counts(i)/totalTrials)*100;
        textStr = sprintf('%d (%.2f%%)', counts(i), percentage);
        text(i, counts(i) + max(counts)*0.05, textStr, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

title('Trial Outcomes'); % Adding the title here

%% Percentages
% Assuming allFlagID, allBehaviorOutcomes, and allBehaviorOutcomesSTR are column vectors of equal length

% Filter out 'PreFixBreak' and 'PreStimBreak' from allBehaviorOutcomesSTR
valid_indices = ~(strcmp(allBehaviorOutcomesSTR, 'PreFixBreak') | strcmp(allBehaviorOutcomesSTR, 'PreStimBreak'));

filteredFlagID = allFlagID(valid_indices);
filteredBehaviorOutcomes = allBehaviorOutcomes(valid_indices);

% Calculation for allFlagID value of 0
trials_0 = find(filteredFlagID == 0);
success_trials_0 = find(filteredBehaviorOutcomes(trials_0) == 1);
percentage_0 = (length(success_trials_0) / length(trials_0)) * 100;

% Calculation for allFlagID value of 1
trials_1 = find(filteredFlagID == 1);
success_trials_1 = find(filteredBehaviorOutcomes(trials_1) == 1);
percentage_1 = (length(success_trials_1) / length(trials_1)) * 100;

% Calculation for allFlagID value of 2
trials_2 = find(filteredFlagID == 2);
success_trials_2 = find(filteredBehaviorOutcomes(trials_2) == 1);
percentage_2 = (length(success_trials_2) / length(trials_2)) * 100;

% Display the results
fprintf('Percentage of successful trials for allFlagID 0 (excluding specified outcomes): %.2f%%\n', percentage_0);
fprintf('Percentage of successful trials for allFlagID 1 (excluding specified outcomes): %.2f%%\n', percentage_1);
fprintf('Percentage of successful trials for allFlagID 2 (excluding specified outcomes): %.2f%%\n', percentage_2);


%% histogram of Fixation Times
histogram(allFixationTime,'BinWidth', 10, 'EdgeColor', 'black')

xlabel('Time (ms)');
ylabel('Number of Trials');
title('Fixation Time Distribution');


%% Histogram of all Response times
% Assuming 'allResponseTime' is a vector containing your data

% Filter out the values of 0 from the data
filteredResponseTime = allResponseTime(allResponseTime ~= 0);

% Create the histogram with the filtered data
histogram(filteredResponseTime, 'EdgeColor', 'black');

% Add labels and title
xlabel('Time (ms)');
ylabel('Number of Trials');
title('Response Time Distribution');


%% Histogram of coordinates for each trial
histogram(allCondID,'BinWidth', 1, 'EdgeColor', 'black')





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Heat Scatters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Finding averages

% Unique target location indices
uniqueTargets = unique(allCondID);

% Initialize variables
CorrectNoPhosRT = zeros(size(uniqueTargets));
CorrectPhosRT = zeros(size(uniqueTargets));
numCorrectNoPhos = zeros(size(uniqueTargets));
numCorrectPhos = zeros(size(uniqueTargets));

% Loop over unique target locations
for i = 1:length(uniqueTargets)
    % Get indices for correct no phosphene saccade trials for the current location
    idxNoPhos = allCondID == uniqueTargets(i) & allFlagID == 1 & allBehaviorOutcomes == 1;
    
    % Get indices for correct phosphene saccade trials for the current location
    idxPhos = allCondID == uniqueTargets(i) & allFlagID == 2 & allBehaviorOutcomes == 1;
    
    % Calculate averages and store them
    CorrectNoPhosRT(i) = mean(allResponseTime(1, idxNoPhos));
    CorrectPhosRT(i) = mean(allResponseTime(1, idxPhos));
    
    % Count the number of correct trials of each type and store them
    numCorrectNoPhos(i) = sum(idxNoPhos);
    numCorrectPhos(i) = sum(idxPhos);
end


% Calculate ratio of average response times
CorrectRatio = CorrectPhosRT ./ CorrectNoPhosRT;

%% Scattter of coordinates


OrginTargetLocations = vertcat(TargetLocations, [0 0]);



scatter(OrginTargetLocations(:,1),OrginTargetLocations(:,2))



% Add this line to make each cell square
axis equal;
%% scatterplot heat maps
close all


% Create a new figure
figure;

% Define marker size
markerSize = 100; % adjust as needed

% Create a scatter plot for No Phosphene Saccade trials
subplot(3,1,1);
scatter(TargetLocations(:,1), TargetLocations(:,2), markerSize, CorrectNoPhosRT, 'filled');
title('Average Response Time for Correct No Phosphene Saccade Trials');
colorbar; % Show color scale
xlabel('X-Coordinate (dva)');
ylabel('Y-Coordinate (dva)');
% Add trial counts as text annotations
text(TargetLocations(:,1), TargetLocations(:,2), string(numCorrectNoPhos), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% Create a scatter plot for Phosphene Saccade trials
subplot(3,1,2);
scatter(TargetLocations(:,1), TargetLocations(:,2), markerSize, CorrectPhosRT, 'filled');
title('Average Response Time for Correct Phosphene Saccade Trials');
colorbar; % Show color scale
xlabel('X-Coordinate (dva)');
ylabel('Y-Coordinate (dva)');
% Add trial counts as text annotations
text(TargetLocations(:,1), TargetLocations(:,2), string(numCorrectPhos), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% Create a scatter plot for ratio of Phosphene Saccade to No Phosphene Saccade trials
subplot(3,1,3);
scatter(TargetLocations(:,1), TargetLocations(:,2), markerSize, CorrectRatio, 'filled');
title('Ratio of Average Response Times for Correct Phosphene to No Phosphene Saccade Trials');
colorbar; % Show color scale
xlabel('X-Coordinate (dva)');
ylabel('Y-Coordinate (dva)');

saveas(gcf, 'HeatScatterRTs_ratio', 'jpeg');

%%
TrialCounts = [numCorrectPhos,numCorrectNoPhos]

mean(TrialCounts)
std(TrialCounts)

%% Curve Eccentrisity

MeridianIndex = [4 5 6 7 8 9 10]
EccentricityIndex = [1 2 3 7 11 12 13]

numpointsRadial =7 
numpointsArc = 7
lengthRadial = 4
lengthArc = 4

% The range starts from negative half of lengthRadial to positive half of lengthRadial
startValue = -lengthRadial / 2;
endValue = lengthRadial / 2;
% Create the distanceMeridian array with numpointsRadial number of points
distanceRadial = linspace(startValue, endValue, numpointsRadial);

figure;

% Select only the elements from 'CorrectPhosRT' specified in 'MeridianIndex'
meridian_CorrectPhosRT = CorrectPhosRT(MeridianIndex);

% 'k-' specifies a black line, 'MarkerFaceColor','k' specifies black filled markers,
% and 'LineWidth',2 specifies a line 2 points wide.
plot(distanceRadial, meridian_CorrectPhosRT, 'k-o', 'MarkerFaceColor','k', 'LineWidth',2);

xlabel('Meridian');
ylabel('Reaction Time');


%% Curve Meridian

figure;


% The range starts from negative half of lengthRadial to positive half of lengthRadial
startValue = -lengthArc / 2;
endValue = lengthArc / 2;
% Create the distanceMeridian array with numpointsRadial number of points
distanceArc = linspace(startValue, endValue, numpointsArc);



% Select only the elements from 'CorrectPhosRT' specified in 'EccentricityIndex'
radial_CorrectPhosRT = CorrectPhosRT(EccentricityIndex);

% 'k-' specifies a black line, 'MarkerFaceColor','k' specifies black filled markers,
% and 'LineWidth',2 specifies a line 2 points wide.
plot(distanceArc, radial_CorrectPhosRT, 'k-o', 'MarkerFaceColor','k', 'LineWidth',2);

xlabel('Eccentricity');
ylabel('Reaction Time');


%% Stats

%% extracting the values of each trial percondition and location

% Unique target location indices
uniqueTargets = unique(allCondID);

% Initialize cell arrays to store trial reaction times for each location
TrialNoPhosRTs = cell(size(uniqueTargets));
TrialPhosRTs = cell(size(uniqueTargets));

% Loop over unique target locations
for i = 1:length(uniqueTargets)
    % Get indices for correct no phosphene saccade trials for the current location
    idxNoPhos = allCondID == uniqueTargets(i) & allFlagID == 1 & allBehaviorOutcomes == 1;
    
    % Get indices for correct phosphene saccade trials for the current location
    idxPhos = allCondID == uniqueTargets(i) & allFlagID == 2 & allBehaviorOutcomes == 1;
    
    % Extract trial reaction times for both conditions and store them
    TrialNoPhosRTs{i} = allResponseTime(1, idxNoPhos);
    TrialPhosRTs{i} = allResponseTime(1, idxPhos);
end






% Assuming TrialPhosRTs and TrialNoPhosRTs are already created as cell arrays.

% Get the trial values for location 7
locationIndex = 7;
phosRTsForLocation7 = TrialPhosRTs{locationIndex};
noPhosRTsForLocation7 = TrialNoPhosRTs{locationIndex};



%% Tessting normality and varience

%normality
% Perform Jarque-Bera test for normality for data at location 7 with phosphene
locationIndex = 7;
phosRTsForLocation7 = TrialPhosRTs{locationIndex};
[hPhos, pValuePhos, ~, ~] = jbtest(phosRTsForLocation7);

% Perform Jarque-Bera test for normality for data at location 7 without phosphene
noPhosRTsForLocation7 = TrialNoPhosRTs{locationIndex};
[hNoPhos, pValueNoPhos, ~, ~] = jbtest(noPhosRTsForLocation7);

% Display the p-values for both groups
disp(['P-value (Phosphene): ', num2str(pValuePhos)]);
disp(['P-value (No Phosphene): ', num2str(pValueNoPhos)]);



%varience
% Perform F-test for equal variances
[~, pValueEqualVar, ~] = vartest2(phosRTsForLocation7, noPhosRTsForLocation7);

% Display the p-value
disp(['P-value (Equal Variance): ', num2str(pValueEqualVar)]);


%% Perform the two-sample t-teston orgin
[~, pValue, ~, stats] = ttest2(phosRTsForLocation7, noPhosRTsForLocation7);

% Display the p-value
disp(['P-value: ', num2str(pValue)]);

% Check if the means are significantly different at the 95% confidence level
alpha = 0.05;
if pValue < alpha
    disp('The means are significantly different.');
else
    disp('The means are not significantly different.');
end



%% 
%% all locaton assumption and ttests

% Unique target location indices
uniqueTargets = unique(allCondID);

% Initialize cell arrays to store results for each location
assumptionResults = cell(length(uniqueTargets), 1);
pValuesWilcoxon = zeros(length(uniqueTargets), 1);

% Loop over unique target locations
for i = 1:length(uniqueTargets)
    % Get the trial reaction times for the current location for both conditions
    phosRTsForLocation = TrialPhosRTs{i};
    noPhosRTsForLocation = TrialNoPhosRTs{i};
    
    % Perform Jarque-Bera test for normality for data at the current location with phosphene
    [~, pValuePhos, ~, ~] = jbtest(phosRTsForLocation);
    
    % Perform Jarque-Bera test for normality for data at the current location without phosphene
    [~, pValueNoPhos, ~, ~] = jbtest(noPhosRTsForLocation);
    
    % Perform F-test for equal variances
    [~, pValueEqualVar, ~] = vartest2(phosRTsForLocation, noPhosRTsForLocation);
    
    % Perform the Wilcoxon rank-sum test (Mann-Whitney U test)
    [pValueWilcoxon, ~, ~] = ranksum(phosRTsForLocation, noPhosRTsForLocation);
    
    % Save the p-values for assumption testing and Wilcoxon rank-sum test to a separate variable
    assumptionResults{i} = struct('Location', uniqueTargets(i), ...
                                  'PValueNormalityPhos', pValuePhos, ...
                                  'PValueNormalityNoPhos', pValueNoPhos, ...
                                  'PValueEqualVariance', pValueEqualVar);
    pValuesWilcoxon(i) = pValueWilcoxon;
end

% Save the p-values for the Wilcoxon rank-sum test to a separate variable
disp('Wilcoxon Rank-Sum Test P-values:');
disp(pValuesWilcoxon);

% Display the p-values for assumption testing at each location
for i = 1:length(uniqueTargets)
    disp(['P-values for Assumption Testing at Location ', num2str(uniqueTargets(i)), ':']);
    disp(assumptionResults{i});
    disp(' ');
end

% Alternatively, you can save the assumptionResults and pValuesWilcoxon variables to a file for later use
% save('assumption_results.mat', 'assumptionResults', 'pValuesWilcoxon');


%%

% Assuming phosRTsForLocation7 and noPhosRTsForLocation7 contain the data for the two groups.

% Perform the Wilcoxon rank-sum test (Mann-Whitney U test)
[pValueWilcoxon, ~, stats] = ranksum(phosRTsForLocation7, noPhosRTsForLocation7);

% Display the p-value
disp(['P-value (Wilcoxon rank-sum test): ', num2str(pValueWilcoxon)]);
stats

%% extracting locations

% Unique target location indices
uniqueTargets = unique(allCondID);

% Initialize variables
CorrectNoPhosRT = zeros(size(uniqueTargets));
CorrectPhosRT = zeros(size(uniqueTargets));
numCorrectNoPhos = zeros(size(uniqueTargets));
numCorrectPhos = zeros(size(uniqueTargets));

% Loop over unique target locations
for i = 1:length(uniqueTargets)
    % Get indices for correct no phosphene saccade trials for the current location
    idxNoPhos = allCondID == uniqueTargets(i) & allFlagID == 1 & allBehaviorOutcomes == 1;
    
    % Get indices for correct phosphene saccade trials for the current location
    idxPhos = allCondID == uniqueTargets(i) & allFlagID == 2 & allBehaviorOutcomes == 1;
    
    % Calculate averages and store them
    CorrectNoPhosRT(i) = mean(allResponseTime(1, idxNoPhos));
    CorrectPhosRT(i) = mean(allResponseTime(1, idxPhos));
    
    % Count the number of correct trials of each type and store them
    numCorrectNoPhos(i) = sum(idxNoPhos);
    numCorrectPhos(i) = sum(idxPhos);
end

% Calculate ratio of average response times
CorrectRatio = CorrectPhosRT ./ CorrectNoPhosRT;

% Extract the p-values from the Wilcoxon rank-sum test for each location
pValuesWilcoxon = zeros(size(uniqueTargets));
for i = 1:length(uniqueTargets)
    % Get the trial reaction times for the current location for both conditions
    phosRTsForLocation = TrialPhosRTs{i};
    noPhosRTsForLocation = TrialNoPhosRTs{i};
    
    % Perform the Wilcoxon rank-sum test (Mann-Whitney U test)
    [pValueWilcoxon, ~, ~] = ranksum(phosRTsForLocation, noPhosRTsForLocation);
    
    % Save the p-value for the current location
    pValuesWilcoxon(i) = pValueWilcoxon;
end

%% all curves together

MeridianIndex = [4 5 6 7 8 9 10];
EccentricityIndex = [1 2 3 7 11 12 13];

numpointsRadial = 7;
numpointsArc = 7;
lengthRadial = 4;
lengthArc = 4;

% The range starts from negative half of lengthRadial to positive half of lengthRadial
startValue = -lengthRadial / 2;
endValue = lengthRadial / 2;
% Create the distanceMeridian array with numpointsRadial number of points
distanceRadial = linspace(startValue, endValue, numpointsRadial);

% Select only the elements from 'CorrectPhosRT' specified in 'MeridianIndex'
meridian_CorrectPhosRT = CorrectPhosRT(MeridianIndex);

% Select only the elements from 'CorrectNoPhosRT' specified in 'MeridianIndex'
meridian_CorrectNoPhosRT = CorrectNoPhosRT(MeridianIndex);

% Create a new figure for the combined plot
figure;

% Plot the first set of data (Meridian)
subplot(2, 1, 1); % Two rows, one column, first subplot
plot(distanceRadial, meridian_CorrectPhosRT, 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2);
hold on;
plot(distanceRadial, meridian_CorrectNoPhosRT, 'k--o', 'MarkerFaceColor', 'k', 'LineWidth', 2);

% Add '*' above the points with significant p-value in the Meridian plot for CorrectPhosRT
alpha = 0.05;
for i = 1:length(MeridianIndex)
    if pValuesWilcoxon(MeridianIndex(i)) < alpha
        plot(distanceRadial(i), meridian_CorrectPhosRT(i), 'k*', 'MarkerSize', 12); % Increased MarkerSize to 12
        text(distanceRadial(i), meridian_CorrectPhosRT(i), '*', 'FontWeight', 'bold', 'FontSize', 12, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

hold off;
xlabel('Meridian (deg.)');
ylabel('Reaction Time (ms)');

% Add a legend to the top subplot
legend('Phosphene', 'No Phosphene', 'Significant', 'Location', 'best'); % Added 'Significant' entry to the legend

% Plot the second set of data (Eccentricity)
subplot(2, 1, 2); % Two rows, one column, second subplot

% The range starts from negative half of lengthArc to positive half of lengthArc
startValue = -lengthArc / 2;
endValue = lengthArc / 2;
% Create the distanceArc array with numpointsArc number of points
distanceArc = linspace(startValue, endValue, numpointsArc);

% Select only the elements from 'CorrectPhosRT' specified in 'EccentricityIndex'
radial_CorrectPhosRT = CorrectPhosRT(EccentricityIndex);

% Select only the elements from 'CorrectNoPhosRT' specified in 'EccentricityIndex'
radial_CorrectNoPhosRT = CorrectNoPhosRT(EccentricityIndex);

plot(distanceArc, radial_CorrectPhosRT, 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2);
hold on;
plot(distanceArc, radial_CorrectNoPhosRT, 'k--o', 'MarkerFaceColor', 'k', 'LineWidth', 2);

% Add '*' above the points with significant p-value in the Eccentricity plot for CorrectPhosRT
for i = 1:length(EccentricityIndex)
    if pValuesWilcoxon(EccentricityIndex(i)) < alpha
        plot(distanceArc(i), radial_CorrectPhosRT(i), 'k*', 'MarkerSize', 12); % Increased MarkerSize to 12
        text(distanceArc(i), radial_CorrectPhosRT(i), '*', 'FontWeight', 'bold', 'FontSize', 12, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

hold off;
xlabel('Eccentricity (deg.)');
ylabel('Reaction Time (ms)');


%% Diffrence curve



% Select only the elements from 'CorrectPhosRT' specified in 'MeridianIndex'
meridian_CorrectPhosRT = CorrectPhosRT(MeridianIndex);

% Select only the elements from 'CorrectNoPhosRT' specified in 'MeridianIndex'
meridian_CorrectNoPhosRT = CorrectNoPhosRT(MeridianIndex);

% Finding the diffrence in latency
meridian_DiffRT = meridian_CorrectPhosRT-meridian_CorrectNoPhosRT

% Create a new figure for the combined plot
figure;

% Plot the first set of data (Meridian)
subplot(2, 1, 1); % Two rows, one column, first subplot
plot(distanceRadial, meridian_DiffRT, 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2);
hold on;


hold off;
xlabel('Meridian (deg.)');
ylabel('Latency Difference (ms)');

% Add a legend to the top subplot
legend('Latency Difference', 'No Phosphene', 'Significant', 'Location', 'best'); % Added 'Significant' entry to the legend

% Plot the second set of data (Eccentricity)
subplot(2, 1, 2); % Two rows, one column, second subplot


% Select only the elements from 'CorrectPhosRT' specified in 'EccentricityIndex'
radial_CorrectPhosRT = CorrectPhosRT(EccentricityIndex);

% Select only the elements from 'CorrectNoPhosRT' specified in 'EccentricityIndex'
radial_CorrectNoPhosRT = CorrectNoPhosRT(EccentricityIndex);


% Finding the diffrence in latency
radial_DiffRT = radial_CorrectPhosRT-radial_CorrectNoPhosRT

plot(distanceArc, radial_DiffRT, 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2);
hold on;



hold off;
xlabel('Eccentricity (deg.)');
ylabel('Latency Difference (ms)');












%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    ALLLL GRAPHS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


offset = 5; % You can adjust this value based on your graph's scale



% Initialize variables
SE_NoPhos = zeros(size(uniqueTargets));
SE_Phos = zeros(size(uniqueTargets));

% Loop over unique target locations
for i = 1:length(uniqueTargets)
    % Indices similar as before
    idxNoPhos = allCondID == uniqueTargets(i) & allFlagID == 1 & allBehaviorOutcomes == 1;
    idxPhos = allCondID == uniqueTargets(i) & allFlagID == 2 & allBehaviorOutcomes == 1;
    
    % Calculate Standard Deviations
    SD_NoPhos = std(allResponseTime(1, idxNoPhos));
    SD_Phos = std(allResponseTime(1, idxPhos));
    
    % Calculate Standard Errors (which is SEM)
    SE_NoPhos(i) = SD_NoPhos / sqrt(numCorrectNoPhos(i));
    SE_Phos(i) = SD_Phos / sqrt(numCorrectPhos(i));
end

% Extracting Standard Errors for Meridian and Eccentricity
meridian_SE_NoPhos = SE_NoPhos(MeridianIndex);
meridian_SE_Phos = SE_Phos(MeridianIndex);
eccentricity_SE_NoPhos = SE_NoPhos(EccentricityIndex);
eccentricity_SE_Phos = SE_Phos(EccentricityIndex);

% Plotting with Error Bars (SEM)
figure;

% Meridian Plot with Error Bars
subplot(2, 2, 1);
errorbar(distanceRadial, meridian_CorrectPhosRT, meridian_SE_Phos, 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2);
hold on;
errorbar(distanceRadial, meridian_CorrectNoPhosRT, meridian_SE_NoPhos, 'k--o', 'MarkerFaceColor', 'k', 'LineWidth', 2);
hold off;
xlabel('Meridian (deg.)');
ylabel('Reaction Time (ms)');
title('All Curves Together - Meridian');

% Eccentricity Plot with Error Bars
subplot(2, 2, 3);
errorbar(distanceArc, radial_CorrectPhosRT, eccentricity_SE_Phos, 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2);
hold on;
errorbar(distanceArc, radial_CorrectNoPhosRT, eccentricity_SE_NoPhos, 'k--o', 'MarkerFaceColor', 'k', 'LineWidth', 2);
hold off;
xlabel('Eccentricity (deg.)');
ylabel('Reaction Time (ms)');
title('All Curves Together - Eccentricity');

% ... [other subplots] ...


% THIRD SUBPLOT - top right (difference curve - Meridian)
subplot(2, 2, 2);
plot(distanceRadial, meridian_DiffRT, 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2);
hold on;


% Add significant points for Meridian if necessary
alpha = 0.05;
for i = 1:length(MeridianIndex)
    if pValuesWilcoxon(MeridianIndex(i)) < alpha
        plot(distanceRadial(i), meridian_DiffRT(i)+offset, 'k*', 'MarkerSize', 12);
    end
end


hold off;
xlabel('Meridian (deg.)');
ylabel('Latency Difference (ms)');
title('Difference Curve - Meridian');

% FOURTH SUBPLOT - bottom right (difference curve - Eccentricity)
subplot(2, 2, 4);
plot(distanceArc, radial_DiffRT, 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2);
hold on;

% Add significant points for Eccentricity if necessary
for i = 1:length(EccentricityIndex)
    if pValuesWilcoxon(EccentricityIndex(i)) < alpha
        plot(distanceArc(i), radial_DiffRT(i)+offset, 'k*', 'MarkerSize', 12);
    end
end

hold off;
xlabel('Eccentricity (deg.)');
ylabel('Latency Difference (ms)');
title('Difference Curve - Eccentricity');
















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Heat Maps  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Heat Map of Trial Counts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_targets = length(TargetLocations);
n_trials = length(allResponseTime);
c_n_trials = zeros(n_targets, 1);  %How many Correct of no Responce trials at each target


%Calculation the responce time for every coordinate
for i = 1:n_trials;
    cond_id = allCondID(i);
    if strcmp(allBehaviorOutcomesSTR(i), 'CorrectResponse') == 1 || strcmp(allBehaviorOutcomesSTR(i), 'NoResponse') == 1 || strcmp(allBehaviorOutcomesSTR(i), 'FlyingOverTime') == 1|| strcmp(allBehaviorOutcomesSTR(i), 'WrongChoice') == 1|| strcmp(allBehaviorOutcomesSTR(i), 'EyeCrossRespWindow') == 1;
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


%non-smoothed heatmap of the CorrentResponse and Noresponse counts

% Create a colormap with an extra color for NaN values
figure;
map = [gray; jet(64)];
colormap(map);

% Create a heatmap of the trial count values
m = max(Z(:))
h = imagesc(x, y, Z, [0 m]); % specify the x, y, Z and the color range
colormap(jet);
cbar = colorbar('Ticks', [0 m*.2 m*.4 m*.6 m*.8 m]);
set(get(cbar, 'Label'), 'String', 'Trial Count');
xlabel('Visual Degrees');
ylabel('Visual Degrees');
title('Trial counts CorrectResponse and NoResponse Trials');
% Flip the y-axis tick labels
yticklabels(flipud(get(gca, 'yticklabels')));


% Set NaN values to be transparent
set(h, 'alphadata', ~isnan(Z));

% Remove NaN values from the colorbar
set(cbar, 'ylim', [0 6]);

saveas(gcf, 'Trial counts of CorrectResponse and NoResponse Trials smooth', 'jpeg');
%% Black and white


% Create a black and white colormap
colormap(gray);

% Create a heatmap of the trial count values
m = max(Z(:));
h = imagesc(x, y, Z, [0 m]); % specify the x, y, Z, and the color range
colorbar('Ticks', [0 m*.2 m*.4 m*.6 m*.8 m]);
set(get(colorbar, 'Label'), 'String', 'Trial Count');
xlabel('Visual Degrees');
ylabel('Visual Degrees');
title('Trial counts CorrectResponse and NoResponse Trials');
yticklabels(flipud(get(gca, 'yticklabels')));

% Set NaN values to be transparent
set(h, 'alphadata', ~isnan(Z));

saveas(gcf, 'Trial counts of CorrectResponse and NoResponse Trials blackwhite', 'jpeg');


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


grid on

saveas(gcf, 'Accuracy CorrectResponse vs NoResponse smooth', 'jpeg');





%% Heat Map of mean response time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_targets = length(TargetLocations);
n_trials = length(allResponseTime);
r_trials = zeros(n_targets, 1);
c_n_trials = zeros(n_targets, 1);
allResponseTime2 = allResponseTime(2,:);

%Calculation the responce time for every coordinate
for i = 1:n_trials;
    cond_id = allCondID(i);
    if strcmp(allBehaviorOutcomesSTR(i), 'CorrectResponse') == 1 || strcmp(allBehaviorOutcomesSTR(i), 'NoResponse') == 1|| strcmp(allBehaviorOutcomesSTR(i), 'FlyingOverTime') == 1;
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

% Create a heatmap of the accuracy values

h = imagesc(x, y, Z, [0 500]); % specify the x, y, Z and the color range
colormap(jet);
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

saveas(gcf, 'Mean Response Time of CorrectResponse and NoResponse Trials smooth', 'jpeg');


%% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  ANOVA S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard deviation of mean per coordiante %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
sdGrid = Z



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


%% Standard deviation per coordinate of accuracy values


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
sdGrid = Z



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
% aveGrid = mean per coordinate
% cnGrid = number of trials
% sdGrid = standard deviation

% get the dimensions of the matrix
[nrows, ncols] = size(accGrid)

% calculate the index of the middle row/column
midrow = floor(nrows/2) + 1
midcol = floor(ncols/2) + 1


%%

% %defining  quadrants
% UR = 1:midrow-1,midcol+1:end
% UL = 1:midrow-1,1:midcol-1
% LR = midrow+1:end,midcol+1:end
% LL = midrow+1:end,1:midcol-1




%% Defining Control and Experimental group data


%defining ACCURACY groups
expAcc = accGrid(1:midrow-1,midcol+1:end)
conAcc = [accGrid(1:midrow-1,1:midcol-1);accGrid(midrow+1:end,midcol+1:end);accGrid(midrow+1:end,1:midcol-1)]
%Make 1d
expAcc = expAcc(:)
conAcc = conAcc(:)

%defining mean groups
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

%defining SD groups
expSd = sdGrid(1:midrow-1,midcol+1:end)
conSd = [sdGrid(1:midrow-1,1:midcol-1);sdGrid(midrow+1:end,midcol+1:end);sdGrid(midrow+1:end,1:midcol-1)]
%Make 1d
expSd = expSd(:)
conSd = conSd(:)


%Data variable names for anova

% expAcc
% conAcc
% 
% expAve
% conAve
% 
% expCn
% conCn
% 
% expSd
% conSd

%% Weighted means of mean times

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

%% Anova of mean times conducted


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
    disp('There is a significant difference between the experimental and control groups');
else
    disp('There is not a significant difference between the experimental and control groups');
end


%%


% Create a 2x2 contingency table of the number of correct and incorrect responses for each group
contingency_table = [sum(expAcc), sum(expCn) - sum(expAcc); sum(conAcc), sum(conCn) - sum(conAcc)];

% Perform a chi-squared test
[h,p,stats] = fishertest(contingency_table);

% If the p-value is less than the chosen significance level (e.g., 0.05), reject the null hypothesis
if p < 0.05
    disp('There is a significant difference in accuracy between the experimental and control groups');
else
    disp('There is not a significant difference in accuracy between the experimental and control groups');
end


%%
% Calculate the standard errors for the experimental and control groups
expSE = expSd ./ sqrt(expCn);
conSE = conSd ./ sqrt(conCn);

% Concatenate the standard errors into a single variable
SE = [expSE; conSE];

% Perform a one-way ANOVA with the 'group' variable as the grouping variable and 'SE' as the weights
[p,tbl,stats] = anova1(data, group, 'varnames', {'Group', 'Measure'}, 'weights', SE);

% Print the ANOVA table
disp(tbl);

% If the p-value is less than the chosen significance level (e.g., 0.05), reject the null hypothesis
if p < 0.05
    disp('There is a significant difference between the experimental and control groups');
else
    disp('There is not a significant difference between the experimental and control groups');
end
