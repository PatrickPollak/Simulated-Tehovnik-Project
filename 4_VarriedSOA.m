clear  all
close all
clc

% TargetLocations:          A two-column matrix where each row corresponds to a unique target location (represented by x and y coordinates). 
% allCondID:                For each trial, an index indicating the row in 'TargetLocations' that corresponds to the target location used.
% allBehaviorOutcomes:      For each trial, a binary value indicating if the trial was successful (1) or not (0).
% allBehaviorOutcomesSTR:   For each trial, a string describing the specific outcome.
% allFixationTime:          For each trial, the duration of fixation prior to the stimulus.
% allResponseTime:          For each trial, two times indicating when in ms after the stimulus the eye left the fixation window (1st row) and when it entered the target window (2nd row).
% allFlagID:                For each trial, a flag indicating the type; 0 for catch trials, 1 for no phosphene saccade trials, and 2 for phosphene saccade trials.
% allSOATime                For each trial, the stimulus onset asyncrony time indicating how many ms before the saccade target the phosphene was shown


% Directory of data files 
SaveDir = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Varried_SOA/Figures';
FileDir = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Varried_SOA/Stimulus_Files';
cd(SaveDir);

% Get a list of all files in the directory
fileList = dir(fullfile(FileDir, '*.mat'));
numFiles = numel(fileList);

% Initialize variables to store the concatenated data
allBehaviorOutcomes = [];
allBehaviorOutcomesSTR = [];
allCondID = [];
allFixationTime = [];
allResponseTime = [];
TargetLocations = [];
allFlagID = [];
allSOATime = [];

% Iterate over the files
for i = 1:numFiles
    % Get the file name and full path
    fileName = fileList(i).name;
    filePath = fullfile(FileDir, fileName);
    
    % Load the data from the file
    data = load(filePath);
    
    % Concatenate the arrays
    allBehaviorOutcomes = [allBehaviorOutcomes, data.allBehaviorOutcomes];
    allBehaviorOutcomesSTR = [allBehaviorOutcomesSTR, data.allBehaviorOutcomesSTR];
    allCondID = [allCondID, data.allCondID];
    allFixationTime = [allFixationTime, data.allFixationTime];
    allResponseTime = [allResponseTime, data.allResponseTime];
    TargetLocations = [TargetLocations, data.TargetLocations];
    allFlagID = [allFlagID, data.allFlagID];
    allSOATime = [allSOATime, data.allSOATime];

end



% Specify the arrays you want to extract but not concatenate
ConsistentVariables = {'TargetLocations', 'TrialConditions', 'lengthArc', 'lengthRadial', 'numRfs', 'numpointsArc', 'numpointsRadial', 'numpointsStraight', 'NumRepeats', 'RFs', 'Times'};

% Load the consistent variables from the first file only
firstFilePath = fullfile(FileDir, fileList(1).name);
data = load(firstFilePath, ConsistentVariables{:});

% Add these variables to the workspace
for j = 1:length(ConsistentVariables)
    assignin('base', ConsistentVariables{j}, data.(ConsistentVariables{j}));
end




% Clearing unnecessary variables
%clearvars -except allBehaviorOutcomes allBehaviorOutcomesSTR allCondID allFixationTime allResponseTime SaveDir TargetLocations;



%% Building allSOATimes
% only for blocks originally w/o allSOATimes

% if ~exist('allSOATimes', 'var')
% 
% % Create a logical index of strings to keep
% toKeep = ~strcmp(allBehaviorOutcomesSTR, 'PreFixBreak') & ~strcmp(allBehaviorOutcomesSTR, 'PreStimBreak');
% 
% % Use this logical index to create the new version without unwanted strings
% newBehaviorOutcomesSTR = allBehaviorOutcomesSTR(toKeep);
% 
% % Taking the soa times
% newSOATime = TrialConditions(5,:);
% dips('AAAAAAAAAH')
% 
% if NumRemainingTrials > 1
% % Only for incomplete blocks
% newSOATime = newSOATime(1:end-NumRemainingTrials);
% else
% end
% 
% % Start with an array of NaN values for allSOATimes
% allSOATimes = NaN(size(allBehaviorOutcomesSTR));
% 
% % Use the toKeep index to insert the newSOATimes into the appropriate locations
% allSOATime(toKeep) = newSOATime;
% end
% 


%% %%%%%%%%%%%%%%%%%% %%
%%% Basic Statistics %%%
%%%% %%%%%%%%%%%%%% %%%%

%Percent Correct trials
sum(allBehaviorOutcomes)/numel(allBehaviorOutcomes)




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram of trial outcomes for varied allFlagID values





% KEY FOR ORDERED_OUTCOMES:
% 'CorrectResponse': Correct response to Saccade Target (ST)
% 'CorrectCatchResponse': Fixation maintained fixation throughout catch trial
% 'PreFixBreak' : Broke fixation in the first 300ms
% 'PreStimBreak' : Broke fixation in the next 200ms
% 'NoStimPhosBreak' : Fixation broke in No phos trial during time taht there would have been (SP)
% 'SimPhosBreak': Fixation broke during display of Simulated Phosphene (SP)
% 'BrokeBlankSOA': Fixation broke during time between SP disapearing and ST appearing or the 500ms of fixation in Catch Trials
% 'CatchBreak': Fixation broke during last 500ms of Catch Trial
% 'RespondTooEarly': Fixation broke within 70ms of ST appearing
% 'NoResponse': Fixation not broken with 500ms of ST appearing
% 'FlyingOvertime': Gaze did not enter ST window within 30ms
% 'WrongChoice': Similar to 'FlyingOvertime' eccept the incorrect response detected by DasCheck function
% 'EyeCrossRespWindow': Eye entered ST window but did not return within acceptable time defined by DasCheck function


% Define relevant flagIDs and corresponding titles
flagIDs = [0, 1, 2];
titles = {'Catch Trial Outcomes', 'No Phosphene Trial Outcomes', 'Phosphene Trial Outcomes'};

% Predefined order for filtered_outcomes
ordered_outcomes = {'CorrectResponse', 'CorrectCatchResponse', 'SimPhosBreak', 'BrokeBlankSOA', ...
                    'CatchBreak', 'RespondTooEarly', 'NoResponse', 'FlyingOvertime', 'WrongChoice', 'EyeCrossRespWindow'};

% Create the figure
figure;

TheYlim = 790

for idx = 1:length(flagIDs)
    % Extracting relevant indices based on current allFlagID value
    relevant_indices = allFlagID == flagIDs(idx) & ...
                       ~(strcmp(allBehaviorOutcomesSTR, 'PreFixBreak') | ...
                         strcmp(allBehaviorOutcomesSTR, 'PreStimBreak') | ...
                         strcmp(allBehaviorOutcomesSTR, 'NoPhosBreak'));
    
    % Count occurrences of each unique string for the current allFlagID
    counts = histcounts(categorical(allBehaviorOutcomesSTR(relevant_indices)), categorical(ordered_outcomes));

    % Compute total trials
    totalTrials = sum(counts);
    
    % Create a subplot for the current allFlagID value
    subplot(3, 1, idx);
    barHandle = bar(counts);
    xticklabels(ordered_outcomes);
    xlabel('Trial Outcomes');
    ylabel('Counts');
    title(titles{idx});
    ylim([0 TheYlim]); % Set y-axis limit
    
    % Rotate xtick labels for better visualization (given there are many categories)
    xtickangle(45);
    
    % Add counts and percentage on top of each bar
    for i = 1:length(counts)
        percentage = (counts(i)/totalTrials)*100;
        textStr = sprintf('%d (%.2f%%)', counts(i), percentage);
        text(i, counts(i) + max(counts)*0.05, textStr, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end


%% Histograms with just correct and not correct

% Define relevant flagIDs and corresponding titles
flagIDs = [0, 1, 2];
titles = {'Catch Trial Outcomes', 'No Phosphene Trial Outcomes', 'Phosphene Trial Outcomes'};

% Create the figure
figure;

% Predefined outcomes for aggregation
correct_outcomes = {'CorrectResponse', 'CorrectCatchResponse'};
incorrect_outcomes = {'SimPhosBreak', 'BrokeBlankSOA', 'CatchBreak', 'RespondTooEarly', 'NoResponse', 'FlyingOvertime', 'WrongChoice', 'EyeCrossRespWindow'};

TheYlim = 790;  % Adjust as needed

for idx = 1:length(flagIDs)
    % Extracting relevant indices based on current allFlagID value
    relevant_indices = allFlagID == flagIDs(idx) & ...
                       ~(strcmp(allBehaviorOutcomesSTR, 'PreFixBreak') | ...
                         strcmp(allBehaviorOutcomesSTR, 'PreStimBreak') | ...
                         strcmp(allBehaviorOutcomesSTR, 'NoPhosBreak'));
    
    % Extracting the data corresponding to the relevant indices
    relevant_data = allBehaviorOutcomesSTR(relevant_indices);
    
    % Counting the occurrences of correct and incorrect outcomes
    correct_count = sum(ismember(relevant_data, correct_outcomes));
    incorrect_count = sum(ismember(relevant_data, incorrect_outcomes));
    aggregated_counts = [correct_count, incorrect_count];
    
    % Compute total trials
    totalTrials = sum(aggregated_counts);
    
    % Create a subplot for the current allFlagID value
    subplot(3, 1, idx);
    barHandle = bar(aggregated_counts);
    xticklabels({'Correct Response', 'Incorrect Response'});
    % xlabel('Trial Outcomes');
    ylabel('Counts');
    title(titles{idx});
    ylim([0 TheYlim]); % Set y-axis limit

    % Add counts and percentage on top of each bar
    for i = 1:length(aggregated_counts)
        percentage = (aggregated_counts(i)/totalTrials)*100;
        textStr = sprintf('%d (%.2f%%)', aggregated_counts(i), percentage);
        text(i, aggregated_counts(i) + max(aggregated_counts)*0.05, textStr, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end



%% Percentages based on allBehaviorOutcomes
% Define relevant flagIDs and corresponding titles
flagIDs = [0, 1, 2];
titles = {'Catch Trial Outcomes', 'No Phosphene Trial Outcomes', 'Phosphene Trial Outcomes'};

% Predefined order for outcomes
ordered_outcomes = {'CorrectResponse', 'CorrectCatchResponse', 'SimPhosBreak', 'BrokeBlankSOA', ...
                    'CatchBreak', 'RespondTooEarly', 'NoResponse', 'FlyingOvertime', 'WrongChoice', 'EyeCrossRespWindow'};

% Initialize results storage
percentage_results = zeros(1, length(flagIDs));

% Iterate over each allFlagID type
for idx = 1:length(flagIDs)
    % Extracting relevant indices based on current allFlagID value
    relevant_indices = allFlagID == flagIDs(idx) & ismember(allBehaviorOutcomesSTR, ordered_outcomes);
    
    % Extract the subset of allBehaviorOutcomes based on relevant indices
    subset_outcomes = allBehaviorOutcomes(relevant_indices);
    
    % Calculate the percentage of these trials that equal 1
    percentage_results(idx) = mean(subset_outcomes == 1) * 100;
end

% Display results
for idx = 1:length(flagIDs)
    disp([titles{idx} ': ' sprintf('%.2f%%', percentage_results(idx))]);
end



%% Percentages based on allBehaviorOutcomesSTR

% Define relevant flagIDs and corresponding titles
flagIDs = [0, 1, 2];
titles = {'Catch Trial Outcomes', 'No Phosphene Trial Outcomes', 'Phosphene Trial Outcomes'};

% Predefined order for outcomes
ordered_outcomes = {'CorrectResponse', 'CorrectCatchResponse', 'SimPhosBreak', 'BrokeBlankSOA', ...
                    'CatchBreak', 'RespondTooEarly', 'NoResponse', 'FlyingOvertime', 'WrongChoice', 'EyeCrossRespWindow'};

% The outcomes we're interested in calculating percentages for
desired_outcomes = {'CorrectResponse', 'CorrectCatchResponse'};

% Initialize results storage
percentage_results = zeros(1, length(flagIDs));

% Iterate over each allFlagID type
for idx = 1:length(flagIDs)
    % Extracting relevant indices based on current allFlagID value
    relevant_indices = allFlagID == flagIDs(idx) & ismember(allBehaviorOutcomesSTR, ordered_outcomes);
    
    % Extract the subset of allBehaviorOutcomesSTR based on relevant indices
    subset_outcomes_str = allBehaviorOutcomesSTR(relevant_indices);
    
    % Calculate the percentage of these trials that are either 'CorrectResponse' or 'CorrectCatchResponse'
    percentage_results(idx) = mean(ismember(subset_outcomes_str, desired_outcomes)) * 100;
end

% Display results
for idx = 1:length(flagIDs)
    disp([titles{idx} ': ' sprintf('%.2f%%', percentage_results(idx))]);
end


%% Which values of allBehavioralOutcomes give a allBehavioralOutcomes equal to 1


% Find indices where allBehaviorOutcomes is 1
indices = find(allBehaviorOutcomes == 1);

% Extract corresponding values from allBehaviorOutcomesSTR
corresponding_values = allBehaviorOutcomesSTR(indices);

% Display unique values and their counts
[unique_values, ~, id] = unique(corresponding_values);
counts = accumarray(id, 1);
for i = 1:length(unique_values)
    disp([unique_values{i} ': ' num2str(counts(i))]);
end



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
%%%%%  Restructuring Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Restructureinf the values Location > SOA

% `TrialRTs` Structure Description:
% 1. NoPhos/Phos: Top-level split representing trials without and with phosphene saccades.
% 2. Location: Sub-categorized by target locations, Fields like 'Loc_1', 'Loc_2', etc
% 3. SOA: Further categorized by Stimulus Onset Asynchrony for each location. Fields like 'SOA_50', 'SOA_100', etc.
% 4. Reaction Times: Vector within each SOA category containing reaction times for trials.
% Access example: `TrialRTs.NoPhos{i}.(sprintf('SOA_%d', SOA))` retrieves reaction times for NoPhos trials, location `i`, and specific SOA value.


% Unique target location indices
uniqueTargets = unique(allCondID);

% Unique SOA times excluding NaN values
uniqueSOATimes = unique(allSOATime(~isnan(allSOATime)));

% Initialize the top-level structure
TrialRTs = struct('NoPhos', [], 'Phos', []);

% Iterate over each unique target location
for i = 1:length(uniqueTargets)
    locationLabel = sprintf('Loc_%d', uniqueTargets(i));
    
    % Initialize cell arrays for current location for each SOA time
    NoPhos_SOAs = cell(size(uniqueSOATimes));
    Phos_SOAs = cell(size(uniqueSOATimes));
    
    for j = 1:length(uniqueSOATimes)
        % Get indices for correct no phosphene saccade trials for the current location and SOA time
        idxNoPhos = allCondID == uniqueTargets(i) & allFlagID == 1 & allBehaviorOutcomes == 1 & allSOATime == uniqueSOATimes(j);
    
        % Get indices for correct phosphene saccade trials for the current location and SOA time
        idxPhos = allCondID == uniqueTargets(i) & allFlagID == 2 & allBehaviorOutcomes == 1 & allSOATime == uniqueSOATimes(j);
    
        % Extract trial reaction times for both conditions
        NoPhos_SOAs{j} = allResponseTime(1, idxNoPhos);
        Phos_SOAs{j} = allResponseTime(1, idxPhos);
    end
    
    % Store the extracted data in the TrialRTs structure
    TrialRTs.NoPhos.(locationLabel) = cell2struct(NoPhos_SOAs, strcat('SOA_', arrayfun(@num2str, uniqueSOATimes, 'UniformOutput', false)), 2);
    TrialRTs.Phos.(locationLabel) = cell2struct(Phos_SOAs, strcat('SOA_', arrayfun(@num2str, uniqueSOATimes, 'UniformOutput', false)), 2);
end



% Restructure the values SOA > LOCATION
% `TrialRTs_Location` Structure Description:
% 1. NoPhos/Phos: Top-level split representing trials without and with phosphene saccades.
% 2. SOA: Sub-categorized by Stimulus Onset Asynchrony. Fields like 'SOA_50', 'SOA_100', etc.
% 3. Location: Further categorized by target locations for each SOA. Fields like 'Loc_1', 'Loc_2', etc.
% 4. Reaction Times: Vector within each location category containing reaction times for trials.
% Access example: `TrialRTs_Location.NoPhos.(sprintf('SOA_%d', SOA)).Loc_1` retrieves reaction times for NoPhos trials, specific SOA value, and location 1.

% Unique target location indices
uniqueTargets = unique(allCondID);

% Unique SOA times excluding NaN values
uniqueSOATimes = unique(allSOATime(~isnan(allSOATime)));

% Initialize the top-level structure
TrialRTs_Location = struct('NoPhos', [], 'Phos', []);

% Iterate over each unique SOA time
for j = 1:length(uniqueSOATimes)
    SOALabel = sprintf('SOA_%d', uniqueSOATimes(j));
    
    % Initialize cell arrays for current SOA for each target location
    NoPhos_Locs = cell(size(uniqueTargets));
    Phos_Locs = cell(size(uniqueTargets));

    for i = 1:length(uniqueTargets)
        locationLabel = sprintf('Loc_%d', uniqueTargets(i));
        
        % Get indices for correct no phosphene saccade trials for the current SOA time and location
        idxNoPhos = allCondID == uniqueTargets(i) & allFlagID == 1 & allBehaviorOutcomes == 1 & allSOATime == uniqueSOATimes(j);
    
        % Get indices for correct phosphene saccade trials for the current SOA time and location
        idxPhos = allCondID == uniqueTargets(i) & allFlagID == 2 & allBehaviorOutcomes == 1 & allSOATime == uniqueSOATimes(j);
    
        % Extract trial reaction times for both conditions
        NoPhos_Locs{i} = allResponseTime(1, idxNoPhos);
        Phos_Locs{i} = allResponseTime(1, idxPhos);
    end
    
    % Store the extracted data in the TrialRTs_Location structure
    TrialRTs_Location.NoPhos.(SOALabel) = cell2struct(NoPhos_Locs, strcat('Loc_', arrayfun(@num2str, uniqueTargets, 'UniformOutput', false)), 2);
    TrialRTs_Location.Phos.(SOALabel) = cell2struct(Phos_Locs, strcat('Loc_', arrayfun(@num2str, uniqueTargets, 'UniformOutput', false)), 2);
end

%% findng the counts



%% Heat Scatterplot of RTs and 

% Define the SOA you want to visualize for heat scatterplot
SOA = 1000;

% Assuming `TrialRTs` has fields for each unique target
locationFields = fieldnames(TrialRTs.NoPhos);

% Loop over unique target locations
for i = 1:length(locationFields)
    
    location = locationFields{i}; % Extract current location name
    
    % Extract reaction times for the specific SOA and target location for No Phosphene trials
    rt_NoPhos = TrialRTs.NoPhos.(location).(sprintf('SOA_%d', SOA));
    
    % Extract reaction times for the specific SOA and target location for Phosphene trials
    rt_Phos = TrialRTs.Phos.(location).(sprintf('SOA_%d', SOA));

    % Calculate averages and store them
    CorrectRT_NoPhos(i) = mean(rt_NoPhos);
    CorrectRT_Phos(i) = mean(rt_Phos);

    % Count the number of correct trials of each type and store them
    numCorrect_NoPhos(i) = length(rt_NoPhos);
    numCorrect_Phos(i) = length(rt_Phos);
end


% Calculate ratio of average response times
CorrectDiffrence = CorrectRT_Phos - CorrectRT_NoPhos;


% Heat scatterplot based on the designated SOA value 
close all

% Create a new figure
figure;

% Define marker size
markerSize = 100; % adjust as needed

% Create a scatter plot for No Phosphene Saccade trials
subplot(3,1,1);
scatter(TargetLocations(:,1), TargetLocations(:,2), markerSize, CorrectRT_NoPhos, 'filled');
title(['Average Response Time for Correct No Phosphene Saccade Trials - SOA: ' num2str(SOA) 'ms']);
colorbar; % Show color scale
xlabel('X-Coordinate (dva)');
ylabel('Y-Coordinate (dva)');
% Add trial counts as text annotations
text(TargetLocations(:,1), TargetLocations(:,2), string(numCorrect_NoPhos), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% Create a scatter plot for Phosphene Saccade trials
subplot(3,1,2);
scatter(TargetLocations(:,1), TargetLocations(:,2), markerSize, CorrectRT_Phos, 'filled');
title(['Average Response Time for Correct Phosphene Saccade Trials - SOA: ' num2str(SOA) 'ms']);
colorbar; % Show color scale
xlabel('X-Coordinate (dva)');
ylabel('Y-Coordinate (dva)');
% Add trial counts as text annotations
text(TargetLocations(:,1), TargetLocations(:,2), string(numCorrect_Phos), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% Create a scatter plot for Diffrence of Phosphene Saccade to No Phosphene Saccade trials
subplot(3,1,3);
scatter(TargetLocations(:,1), TargetLocations(:,2), markerSize, CorrectDiffrence, 'filled');
title(['Diffrence of Average Response Times for Correct Phosphene to No Phosphene Saccade Trials - SOA: ' num2str(SOA) 'ms']);
colorbar; % Show color scale
xlabel('X-Coordinate (dva)');
ylabel('Y-Coordinate (dva)');

% Save the figure
saveas(gcf, ['HeatScatterRTs_Diffrence_SOA_' num2str(SOA) 'ms'], 'jpeg');


%% trial average
% Collect all trial counts into a single array
all_trial_counts = [];

% Iterate over each condition in the structure
conditions = fieldnames(TrialRTs);
for i = 1:length(conditions)
    condition = conditions{i};
    locations = fieldnames(TrialRTs.(condition));
    
    % Iterate over locations
    for j = 1:length(locations)
        location = locations{j};
        soas = fieldnames(TrialRTs.(condition).(location));
        
        % Iterate over SOAs
        for k = 1:length(soas)
            soa = soas{k};
            rt_data = TrialRTs.(condition).(location).(soa);
            
            % Append the number of trials for the current SOA to the array
            all_trial_counts = [all_trial_counts, length(rt_data)];
        end
    end
end

% Calculate the overall average and standard deviation for all trial counts
overall_average = mean(all_trial_counts);
overall_std_deviation = std(all_trial_counts);

% Display the overall average and standard deviation
disp(['Overall average number of trials: ', num2str(overall_average)]);
disp(['Overall standard deviation of trial counts: ', num2str(overall_std_deviation)]);



%% Define the Location, Each SOA RTs plotted

close all

% Define the location you want to visualize
Loc_ = 7; % Change this to specify the location.

% Convert the location to string format for field referencing
locationField = sprintf('Loc_%d', Loc_);

% Get SOA values and calculate the middle location index
uniqueSOAs = flip(unique(allSOATime(~isnan(allSOATime))));

% Extract data for the specified location
data_Phos = TrialRTs.Phos.(locationField);
data_NoPhos = TrialRTs.NoPhos.(locationField);

% Initialization
meanRTs_Phos = zeros(size(uniqueSOAs));
semRTs_Phos = zeros(size(uniqueSOAs));
meanRTs_NoPhos = zeros(size(uniqueSOAs));
semRTs_NoPhos = zeros(size(uniqueSOAs));

for i = 1:length(uniqueSOAs)
    currentSOA = uniqueSOAs(i);
    
    % For Phos trials
    rt_Phos = data_Phos.(sprintf('SOA_%d', currentSOA));
    meanRTs_Phos(i) = mean(rt_Phos);
    semRTs_Phos(i) = std(rt_Phos) / sqrt(length(rt_Phos));
    
    % For NoPhos trials
    rt_NoPhos = data_NoPhos.(sprintf('SOA_%d', currentSOA));
    meanRTs_NoPhos(i) = mean(rt_NoPhos);
    semRTs_NoPhos(i) = std(rt_NoPhos) / sqrt(length(rt_NoPhos));
end

% Plotting

figure;

% Define positions for each SOA on x-axis
xPos = 1:length(uniqueSOAs);

% Combined Phos and NoPhos Plot
subplot(2, 1, 1);
errorbar(xPos, meanRTs_Phos, semRTs_Phos, '-o', 'Color', 'r', 'LineWidth', 2); % Phos in red
hold on;
errorbar(xPos, meanRTs_NoPhos, semRTs_NoPhos, '-o', 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--'); % NoPhos in black dashed line
set(gca, 'XTick', xPos, 'XTickLabel', uniqueSOAs);
xlabel('SOA (ms)');
ylabel('Reaction Time (ms)');
title(sprintf('Mean Reaction Time for Trials at Location %d per SOA', Loc_));
legend('Phosphene', 'No Phosphene');

hold off;

% Differences Plot
subplot(2, 1, 2);
diffRTs = meanRTs_Phos - meanRTs_NoPhos;
plot(xPos, diffRTs, '-o', 'Color', 'k', 'LineWidth', 2);
set(gca, 'XTick', xPos, 'XTickLabel', uniqueSOAs);
xlabel('SOA (ms)');
ylabel('Reaction Time (ms)');
title('Difference between Phosphene Trials and No Phosphene Trials');


% Rank-sum test and asterisk placement
for k = 1:length(uniqueSOAs)
    p = ranksum(TrialRTs.Phos.(locationField).(sprintf('SOA_%d', uniqueSOAs(k))), ...
            TrialRTs.NoPhos.(locationField).(sprintf('SOA_%d', uniqueSOAs(k))));

    if p < 0.05
        text(xPos(k), diffRTs(k) + max(diffRTs) * 0.05, '*', 'FontSize', 30, 'HorizontalAlignment', 'center', 'Color', 'k');
    end
end

% Generate the filename dynamically
filename = sprintf('SOA_Loc_%d.jpg', Loc_);

% Save the current figure as a JPEG with the dynamic filename
saveas(gcf, filename, 'jpeg');



%% figure description at the intersection point

% Figure Y. Temporal Dynamics of Reaction Times at the Simulated Phosphene Location.

% Top Panel: Mean reaction times (RTs) at the location where the simulated phosphene is presented, coinciding with the intersection of the meridian and eccentricity coordinates in the visual field, across various stimulus onset asynchronies (SOAs). Two conditions are contrasted: trials with phosphene-induced saccades (represented by a solid red line with circular markers) and non-phosphene trials (depicted as a dashed black line with circular markers). Error bars indicate the standard error of the mean. The x-axis labels showcase the unique SOAs used in the study.

% Bottom Panel: This section displays the difference in RTs between the phosphene and non-phosphene trials across the span of SOAs. Significant rank-sum test results (p < 0.05) between the two trial types for specific SOAs are denoted by an asterisk, pointing to a statistically significant difference in distributions between conditions where a phosphene was and wasn't present.

%% Reaction Time per Location

%% Define SOA
SOA = 100; % Define your SOA value here

MeridianIndex = [4 5 6 7 8 9 10];
EccentricityIndex = [1 2 3 7 11 12 13];

% Set up the x-axis coordinates for Meridian
distanceRadial = linspace(-lengthRadial / 2, lengthRadial / 2, length(MeridianIndex));

% Set up the x-axis coordinates for Eccentricity
distanceArc = linspace(-lengthArc / 2, lengthArc / 2, length(EccentricityIndex));

% Data extraction for meridian and eccentricity
SOAField = sprintf('SOA_%d', SOA);

% Initialization
meridian_Phos = zeros(size(MeridianIndex));
sem_meridian_Phos = zeros(size(MeridianIndex));
meridian_NoPhos = zeros(size(MeridianIndex));
sem_meridian_NoPhos = zeros(size(MeridianIndex));

eccentricity_Phos = zeros(size(EccentricityIndex));
sem_eccentricity_Phos = zeros(size(EccentricityIndex));
eccentricity_NoPhos = zeros(size(EccentricityIndex));
sem_eccentricity_NoPhos = zeros(size(EccentricityIndex));

% Extraction using direct indexing
for i = 1:length(MeridianIndex)
    currentLoc = MeridianIndex(i);
    locField = sprintf('Loc_%d', currentLoc);

    % For Phos trials - Meridian
    rt_Phos = TrialRTs_Location.Phos.(SOAField).(locField);
    meridian_Phos(i) = mean(rt_Phos);
    sem_meridian_Phos(i) = std(rt_Phos) / sqrt(length(rt_Phos));

    % For NoPhos trials - Meridian
    rt_NoPhos = TrialRTs_Location.NoPhos.(SOAField).(locField);
    meridian_NoPhos(i) = mean(rt_NoPhos);
    sem_meridian_NoPhos(i) = std(rt_NoPhos) / sqrt(length(rt_NoPhos));
end

for i = 1:length(EccentricityIndex)
    currentLoc = EccentricityIndex(i);
    locField = sprintf('Loc_%d', currentLoc);

    % For Phos trials - Eccentricity
    rt_Phos = TrialRTs_Location.Phos.(SOAField).(locField);
    eccentricity_Phos(i) = mean(rt_Phos);
    sem_eccentricity_Phos(i) = std(rt_Phos) / sqrt(length(rt_Phos));

    % For NoPhos trials - Eccentricity
    rt_NoPhos = TrialRTs_Location.NoPhos.(SOAField).(locField);
    eccentricity_NoPhos(i) = mean(rt_NoPhos);
    sem_eccentricity_NoPhos(i) = std(rt_NoPhos) / sqrt(length(rt_NoPhos));
end

% Plotting Data
figure;

% Meridian subplot
subplot(2, 1, 1);

% Define line handles
lineHandlesMeridian = gobjects(1, 2);

% Plot NoPhos data first, then Phos to ensure it's on top
lineHandlesMeridian(2) = errorbar(distanceArc, meridian_NoPhos, sem_meridian_NoPhos, 'k--o', 'LineWidth', 2); % No Phosphene Trials
hold on;
lineHandlesMeridian(1) = errorbar(distanceArc, meridian_Phos, sem_meridian_Phos, 'r-o', 'LineWidth', 2); % Phosphene Trials
hold off;

title(sprintf('Meridian SOA %dms', SOA));
xlabel('Meridian (deg.)');
ylabel('Reaction Time (ms)');

legend(lineHandlesMeridian, {'Phosphene Trials', 'No Phosphene Trials'}, 'Location', 'eastoutside');

% Eccentricity subplot
subplot(2, 1, 2);

% Define line handles
lineHandlesEccentricity = gobjects(1, 2);

% Plot NoPhos data first, then Phos to ensure it's on top
lineHandlesEccentricity(2) = errorbar(distanceRadial, eccentricity_NoPhos, sem_eccentricity_NoPhos, 'k--o', 'LineWidth', 2); % No Phosphene Trials
hold on;
lineHandlesEccentricity(1) = errorbar(distanceRadial, eccentricity_Phos, sem_eccentricity_Phos, 'r-o', 'LineWidth', 2); % Phosphene Trials
hold off;

title(sprintf('Eccentricity SOA %dms', SOA));
xlabel('Eccentricity (deg.)');
ylabel('Reaction Time (ms)');

legend(lineHandlesEccentricity, {'Phosphene Trials', 'No Phosphene Trials'}, 'Location', 'eastoutside');

% Adjust y-axis limits
meridianYLim = get(subplot(2,1,1), 'YLim');
eccentricityYLim = get(subplot(2,1,2), 'YLim');

commonYLim = [min(meridianYLim(1), eccentricityYLim(1)), max(meridianYLim(2), eccentricityYLim(2))];
set(subplot(2,1,1), 'YLim', commonYLim);
set(subplot(2,1,2), 'YLim', commonYLim);

% Saving the figure with dynamic filename
filename = sprintf('Latency_SOA_%d', SOA);
saveas(gcf, filename, 'png');

%% Plot description

%Figure X. Reaction Time Profiles Across Visual Field

%This figure displays the reaction times (RTs) for two conditions: trials with phosphene-induced saccades ("Phosphene Trials", denoted by solid red lines with circular markers) and trials without phosphene-induced saccades ("No Phosphene Trials", denoted by dashed black lines with circular markers).

%Top Panel (Meridian SOA Xms): Represents the reaction times along the meridian of the visual field. The x-axis (Meridian) spans from negative to positive degrees, capturing the relative orientation. The y-axis denotes the measured reaction times in milliseconds.
%Bottom Panel (Eccentricity SOA Xms): Demonstrates the reaction times based on eccentricity in the visual field. Here, the x-axis (Eccentricity) ranges between negative and positive degrees, representing the deviation from the focal point. Reaction times, in milliseconds, are on the y-axis.
%In both panels, error bars indicate the standard error of the mean for the reaction times at each point in the visual field. A consistent y-axis scale between the panels allows for direct comparison of reaction times across the different sections of the visual field.




%% meridian eccentricity with diffrence underneath

% Plotting Data
figure;

% Meridian subplot
subplot(2, 2, 1);

% Define line handles
lineHandlesMeridian = gobjects(1, 2);

% Plot NoPhos data first, then Phos to ensure it's on top
lineHandlesMeridian(2) = errorbar(distanceArc, meridian_NoPhos, sem_meridian_NoPhos, 'k--', 'LineWidth', 2); % No Phosphene Trials
hold on;
lineHandlesMeridian(1) = errorbar(distanceArc, meridian_Phos, sem_meridian_Phos, 'r-', 'LineWidth', 2); % Phosphene Trials
hold off;

title(sprintf('Meridian SOA %dms', SOA));
xlabel('Meridian (deg.)');
ylabel('Reaction Time (ms)');

% Eccentricity subplot
subplot(2, 2, 2);

% Define line handles
lineHandlesEccentricity = gobjects(1, 2);

% Plot NoPhos data first, then Phos to ensure it's on top
lineHandlesEccentricity(2) = errorbar(distanceRadial, eccentricity_NoPhos, sem_eccentricity_NoPhos, 'k--', 'LineWidth', 2); % No Phosphene Trials
hold on;
lineHandlesEccentricity(1) = errorbar(distanceRadial, eccentricity_Phos, sem_eccentricity_Phos, 'r-', 'LineWidth', 2); % Phosphene Trials
hold off;

title(sprintf('Eccentricity SOA %dms', SOA));
xlabel('Eccentricity (deg.)');
ylabel('Reaction Time (ms)');








% Rank-sum test and difference curve plotting for Meridian
subplot(2, 2, 3);
difference_meridian = meridian_Phos - meridian_NoPhos;  % Calculate the difference
plot(distanceArc, difference_meridian, 'k-', 'LineWidth', 2);  % Plot the difference curve
hold on;
for i = 1:length(MeridianIndex)
    currentLoc = MeridianIndex(i);
    locField = sprintf('Loc_%d', currentLoc);

    rt_Phos = TrialRTs_Location.Phos.(SOAField).(locField);
    rt_NoPhos = TrialRTs_Location.NoPhos.(SOAField).(locField);
    p = ranksum(rt_Phos, rt_NoPhos);
    
    % If p-value < 0.05, place an asterisk above the corresponding data point
    if p < 0.05 
        text(distanceArc(i), difference_meridian(i) + (diffYMax - diffYMin)*0.05, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
    end
end
hold off;
title(sprintf('Difference Meridian SOA %dms', SOA));
xlabel('Meridian (deg.)');
ylabel('Difference in Reaction Time (ms)');

% Rank-sum test and difference curve plotting for Eccentricity
subplot(2, 2, 4);
difference_eccentricity = eccentricity_Phos - eccentricity_NoPhos;  % Calculate the difference
plot(distanceRadial, difference_eccentricity, 'k-', 'LineWidth', 2);  % Plot the difference curve
hold on;
for i = 1:length(EccentricityIndex)
    currentLoc = EccentricityIndex(i);
    locField = sprintf('Loc_%d', currentLoc);

    rt_Phos = TrialRTs_Location.Phos.(SOAField).(locField);
    rt_NoPhos = TrialRTs_Location.NoPhos.(SOAField).(locField);
    p = ranksum(rt_Phos, rt_NoPhos);
    
    % If p-value < 0.05, place an asterisk above the corresponding data point
    if p < 0.05
        text(distanceRadial(i), difference_eccentricity(i) + (diffYMax - diffYMin)*0.05, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
    end
end
hold off;
title(sprintf('Difference Eccentricity SOA %dms', SOA));
xlabel('Eccentricity (deg.)');
ylabel('Difference in Reaction Time (ms)');






% Getting the y-axis limits for top two subplots
subplot(2, 2, 1);
yLimMeridian = ylim;
subplot(2, 2, 2);
yLimEccentricity = ylim;

% Determine the minimum and maximum y values for the top two subplots
topYMin = min(yLimMeridian(1), yLimEccentricity(1));
topYMax = max(yLimMeridian(2), yLimEccentricity(2));

% Applying the synchronized y-axis limits to the top two subplots
subplot(2, 2, 1);
pbaspect([86 72 1]);
ylim([topYMin, topYMax]);
subplot(2, 2, 2);
pbaspect([86 72 1]);
ylim([topYMin, topYMax]);

% Getting the y-axis limits for bottom two subplots
subplot(2, 2, 3);
yLimDiffMeridian = ylim;
subplot(2, 2, 4);
yLimDiffEccentricity = ylim;

% Determine the minimum and maximum y values for the bottom two subplots
bottomYMin = min(yLimDiffMeridian(1), yLimDiffEccentricity(1));
bottomYMax = max(yLimDiffMeridian(2), yLimDiffEccentricity(2));

% Applying the synchronized y-axis limits to the bottom two subplots
subplot(2, 2, 3);
pbaspect([86 72 1]);
ylim([bottomYMin, -bottomYMin]);
subplot(2, 2, 4);
pbaspect([86 72 1]);
ylim([bottomYMin, -bottomYMin]);




%% COMBINED Difference

%% Extract unique SOAs and distances
%% Define distances and get all unique SOAs excluding NaNs
MeridianIndex = [4 5 6 7 8 9 10];
EccentricityIndex = [1 2 3 7 11 12 13];
distanceRadial = linspace(-lengthRadial / 2, lengthRadial / 2, length(MeridianIndex));
distanceArc = linspace(-lengthArc / 2, lengthArc / 2, length(EccentricityIndex));

allSOAs = unique(allSOATime(~isnan(allSOATime)));

legendEntries = cell(1, length(allSOAs));
for idx = 1:length(allSOAs)
    legendEntries{idx} = sprintf('SOA %d ms', allSOAs(idx));
end

% Define the initial and final HSV values
HSV_start = [240, 1, .6];
HSV_end = [40, 0.45, 01];

% Number of colors
nColors = numel(allSOAs);

% Create linearly spaced HSV values
H = linspace(HSV_start(1), HSV_end(1), nColors);
S = linspace(HSV_start(2), HSV_end(2), nColors);
V = linspace(HSV_start(3), HSV_end(3), nColors);

colors = zeros(nColors, 3);

% Convert each HSV value to RGB
for i = 1:nColors
    colors(i, :) = hsv2rgb([H(i)/360, S(i), V(i)]);
end


% Placeholder for line handles
lineHandlesMeridian = gobjects(1, length(allSOAs));
lineHandlesEccentricity = gobjects(1, length(allSOAs));

figure;

% Loop through each unique SOA in reverse
for j = length(allSOAs):-1:1
    SOA = allSOAs(j);
    SOAField = sprintf('SOA_%d', SOA);
    diff_meridian = zeros(size(MeridianIndex));
    diff_eccentricity = zeros(size(EccentricityIndex));
    
    se_meridian = zeros(size(MeridianIndex));
    se_eccentricity = zeros(size(EccentricityIndex));
    
    % Calculate differences for Meridian and Eccentricity
    % Calculate differences for Meridian and Eccentricity
    for i = 1:length(MeridianIndex)
        currentLoc = MeridianIndex(i);
        locField = sprintf('Loc_%d', currentLoc);
        rt_Phos = TrialRTs_Location.Phos.(SOAField).(locField);
        rt_NoPhos = TrialRTs_Location.NoPhos.(SOAField).(locField);
        diff_meridian(i) = mean(rt_Phos) - mean(rt_NoPhos);
        
        % Standard error for meridian
        se_meridian(i) = sqrt((std(rt_Phos)^2/length(rt_Phos)) + (std(rt_NoPhos)^2/length(rt_NoPhos)));
    end

    for i = 1:length(EccentricityIndex)
        currentLoc = EccentricityIndex(i);
        locField = sprintf('Loc_%d', currentLoc);
        rt_Phos = TrialRTs_Location.Phos.(SOAField).(locField);
        rt_NoPhos = TrialRTs_Location.NoPhos.(SOAField).(locField);
        diff_eccentricity(i) = mean(rt_Phos) - mean(rt_NoPhos);
        
        % Standard error for eccentricity
        se_eccentricity(i) = sqrt((std(rt_Phos)^2/length(rt_Phos)) + (std(rt_NoPhos)^2/length(rt_NoPhos)));
    end


    % Meridian subplot
    subplot(1, 2, 1);
    hold on;
    lineHandlesMeridian(j) = errorbar(distanceRadial, diff_meridian, se_meridian, '-', 'LineWidth', 2, 'Color', colors(j, :));
    pbaspect([86 150 1]); % Set the aspect ratio



    % Eccentricity subplot
    subplot(1, 2, 2);
    hold on;
    lineHandlesEccentricity(j) = errorbar(distanceArc, diff_eccentricity, se_eccentricity, '-', 'LineWidth', 2, 'Color', colors(j, :));
    pbaspect([86 150 1]); % Set the aspect ratio
end

% After plotting all the data, add labels, asterisks, and adjust limits

% For Meridian
subplot(1, 2, 1);
title('Meridian Difference');
xlabel('Meridian (deg.)');
ylabel('Reaction Time (ms)');
yline(0, 'k--'); % Add horizontal line at y = 0
legend(lineHandlesMeridian, legendEntries, 'Location', 'eastoutside');
pbaspect([108 188 1]); % Set the aspect ratio


% For Eccentricity
subplot(1, 2, 2);
title('Eccentricity Difference');
xlabel('Eccentricity (deg.)');
ylabel('Reaction Time (ms)');
yline(0, 'k--'); % Add horizontal line at y = 0
legend(lineHandlesEccentricity, legendEntries, 'Location', 'eastoutside');
pbaspect([108 188 1]); % Set the aspect ratio


% Getting the y-axis limits for both subplots
subplot(1, 2, 1);
yLimMeridian = ylim;
subplot(1, 2, 2);
yLimEccentricity = ylim;

% Determining the global minimum and maximum y values
globalYMin = min(yLimMeridian(1), yLimEccentricity(1));
globalYMax = max(yLimMeridian(2), yLimEccentricity(2));

% Determining the global minimum and maximum y values
globalYMin = -50;
globalYMax = max(yLimMeridian(2), yLimEccentricity(2));


% Applying the global y-axis limits to both subplots
subplot(1, 2, 1);
ylim([globalYMin, globalYMax]);
subplot(1, 2, 2);
ylim([globalYMin, globalYMax]);



%% Figure description for this plot

% Figure 1. Difference in reaction times as a function of location and Stimulus Onset Asynchrony (SOA).

% A. Reaction time differences along the meridian. Differences are calculated as the mean reaction time with phosphene saccades minus the mean reaction time without phosphene saccades for each SOA. Positive values indicate faster reactions with phosphene saccades. The x-axis represents spatial location along the meridian in degrees, with the fovea at 0.

% B. Reaction time differences across eccentricities. As with (A), differences are based on phosphene saccades' influence on reaction time. The x-axis represents spatial eccentricity in degrees from the fovea.

% In both (A) and (B), error bars represent the standard error of the difference between independent means for the two conditions. Colors of the lines represent different SOAs, with the spectrum from dark blue to light blue indicating increasing SOA values. The dotted horizontal line represents zero difference, where reaction times are the same for both conditions.



%% Diffrence plots per location

% Define Meridian and Eccentricity Index
MeridianIndex = [4 5 6 7 8 9 10];
EccentricityIndex = [1 2 3 7 11 12 13];

% Calculate distances
distanceRadial = linspace(-lengthRadial / 2, lengthRadial / 2, length(EccentricityIndex));
distanceArc = linspace(-lengthArc / 2, lengthArc / 2, length(MeridianIndex));

% Sort the distances for plotting in descending order
[sortedDistRadial, sortIdxRadial] = sort(distanceRadial, 'descend');
[sortedDistArc, sortIdxArc] = sort(distanceArc, 'descend');

% Define the initial and final HSV values for the color scheme
HSV_start = [240, 1, .6];
HSV_end = [40, 0.45, 1];
nColors = length(EccentricityIndex);  % The number of colors needed

% Generate the colors
H = linspace(HSV_start(1), HSV_end(1), nColors);
S = linspace(HSV_start(2), HSV_end(2), nColors);
V = linspace(HSV_start(3), HSV_end(3), nColors);
colors = zeros(nColors, 3);
for i = 1:nColors
    colors(i, :) = hsv2rgb([H(i)/360, S(i), V(i)]);
end

% Unique SOAs for x-axis
uniqueSOAs = flip(unique(allSOATime(~isnan(allSOATime))));
xPos = 1:length(uniqueSOAs);

% Start plotting
figure;

% Plotting Arc differences (Meridian)
subplot(1, 2, 1);

for idx = 1:length(MeridianIndex)
    i = sortIdxArc(idx);  % Use the sorted indices to plot in descending order
    locationField = sprintf('Loc_%d', MeridianIndex(i));

    diffRTs = zeros(size(uniqueSOAs));
    semDiffRTs = zeros(size(uniqueSOAs));
    for j = 1:length(uniqueSOAs)
        currentSOA = uniqueSOAs(j);
        rt_Phos = TrialRTs.Phos.(locationField).(sprintf('SOA_%d', currentSOA));
        rt_NoPhos = TrialRTs.NoPhos.(locationField).(sprintf('SOA_%d', currentSOA));
        diffRTs(j) = mean(rt_Phos) - mean(rt_NoPhos);
        semDiffRTs(j) = sqrt((var(rt_Phos) + var(rt_NoPhos)) / 2); % Assuming equal number of Phos and NoPhos trials
    end

    errorbar(xPos, diffRTs, semDiffRTs, '-', 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', sprintf('%.1f°', sortedDistArc(idx)));
    hold on;
end

title('Time Differences at all Meridian Locations');
xlabel('SOA (ms)');
ylabel('Reaction Time (ms)');
legend('show');
hold off;
pbaspect([108 188 1]); % Set the aspect ratio
set(gca, 'XTick', xPos, 'XTickLabel', uniqueSOAs); % Set x-tick labels to SOA values

% Plotting Radial differences (Eccentricity)
subplot(1, 2, 2);

for idx = 1:length(EccentricityIndex)
    i = sortIdxRadial(idx);  % Use the sorted indices to plot in descending order
    locationField = sprintf('Loc_%d', EccentricityIndex(i));

    diffRTs = zeros(size(uniqueSOAs));
    semDiffRTs = zeros(size(uniqueSOAs));
    for j = 1:length(uniqueSOAs)
        currentSOA = uniqueSOAs(j);
        rt_Phos = TrialRTs.Phos.(locationField).(sprintf('SOA_%d', currentSOA));
        rt_NoPhos = TrialRTs.NoPhos.(locationField).(sprintf('SOA_%d', currentSOA));
        diffRTs(j) = mean(rt_Phos) - mean(rt_NoPhos);
        semDiffRTs(j) = sqrt((var(rt_Phos) + var(rt_NoPhos)) / 2); % Assuming equal number of Phos and NoPhos trials
    end

    errorbar(xPos, diffRTs, semDiffRTs, '-', 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', sprintf('%.1f°', sortedDistRadial(idx)));
    hold on;
end

title('Time Differences at all Eccentricity Locations');
xlabel('SOA (ms)');
ylabel('Reaction Time (ms)');
legend('show');
hold off;
pbaspect([108 188 1]); % Set the aspect ratio
set(gca, 'XTick', xPos, 'XTickLabel', uniqueSOAs); % Set x-tick labels to SOA values


% Getting the y-axis limits for both subplots
subplot(1, 2, 1);
yLimMeridian = ylim;
subplot(1, 2, 2);
yLimEccentricity = ylim;

% Determining the global minimum and maximum y values
globalYMin = min(yLimMeridian(1), yLimEccentricity(1));
globalYMax = max(yLimMeridian(2), yLimEccentricity(2));


% Determining the global minimum and maximum y values
globalYMin = -70
globalYMax = max(yLimMeridian(2), yLimEccentricity(2));


% Applying the global y-axis limits to both subplots
subplot(1, 2, 1);
ylim([globalYMin, globalYMax]);
subplot(1, 2, 2);
ylim([globalYMin, globalYMax]);





%% Center vs Others


% Locations categorization
CenterLoc = [7];
OtherLoc = [1,2,3,4,5,6,8,9,10,11,12,13];

% Unique SOA times excluding NaN values
uniqueSOATimes = unique(allSOATime(~isnan(allSOATime)));

% Initialize structures to hold combined data
CombinedRTs = struct('NoPhos', struct('Center', [], 'Other', []), 'Phos', struct('Center', [], 'Other', []));

% Iterate over each unique SOA time
for j = 1:length(uniqueSOATimes)
    SOALabel = sprintf('SOA_%d', uniqueSOATimes(j));
    
    % Combine for CenterLoc
    rt_NoPhos_Center = [];
    rt_Phos_Center = [];
    for loc = CenterLoc
        locationLabel = sprintf('Loc_%d', loc);
        
        rt_NoPhos_Center = [rt_NoPhos_Center, TrialRTs.NoPhos.(locationLabel).(SOALabel)];
        rt_Phos_Center = [rt_Phos_Center, TrialRTs.Phos.(locationLabel).(SOALabel)];
    end
    CombinedRTs.NoPhos.Center.(SOALabel) = rt_NoPhos_Center;
    CombinedRTs.Phos.Center.(SOALabel) = rt_Phos_Center;
    
    % Combine for OtherLoc
    rt_NoPhos_Other = [];
    rt_Phos_Other = [];
    for loc = OtherLoc
        locationLabel = sprintf('Loc_%d', loc);
        
        rt_NoPhos_Other = [rt_NoPhos_Other, TrialRTs.NoPhos.(locationLabel).(SOALabel)];
        rt_Phos_Other = [rt_Phos_Other, TrialRTs.Phos.(locationLabel).(SOALabel)];
    end
    CombinedRTs.NoPhos.Other.(SOALabel) = rt_NoPhos_Other;
    CombinedRTs.Phos.Other.(SOALabel) = rt_Phos_Other;
end







%% plots

figure;



locations = {'Center', 'Other'};
colors = {'r', 'b'};  % You can adjust these colors as needed

for locIdx = 1:length(locations)
    loc = locations{locIdx};
    
    diffRTs = zeros(size(uniqueSOAs));
    semDiffRTs = zeros(size(uniqueSOAs));
    for j = 1:length(uniqueSOAs)
        SOALabel = sprintf('SOA_%d', uniqueSOAs(j));
        rt_Phos = CombinedRTs.Phos.(loc).(SOALabel);
        rt_NoPhos = CombinedRTs.NoPhos.(loc).(SOALabel);
        diffRTs(j) = mean(rt_Phos) - mean(rt_NoPhos);
        semDiffRTs(j) = sqrt((var(rt_Phos)/length(rt_Phos) + var(rt_NoPhos)/length(rt_NoPhos)) / 2);
    end
    
    errorbar(xPos, diffRTs, semDiffRTs, '-', 'Color', colors{locIdx}, 'LineWidth', 2, 'DisplayName', loc);
    hold on;
end

title('Time Differences for Center vs Other (Eccentricity)');
xlabel('SOA (ms)');
ylabel('Reaction Time Difference (ms)');
legend('show');
pbaspect([108 188 1]); % Set the aspect ratio
hold off;




set(gca, 'XTick', 1:length(uniqueSOAs), 'XTickLabel', uniqueSOAs);



%%


% Specify the SOA value of interest
SOA = 500;
SOALabel = sprintf('SOA_%d', SOA);

% Calculate the average reaction time for NoPhos trials
avgRT_NoPhos_Center = mean(CombinedRTs.NoPhos.Center.(SOALabel));
avgRT_NoPhos_Other = mean(CombinedRTs.NoPhos.Other.(SOALabel));

% Subtract the NoPhos average from each Phos trial
adjustedRT_Phos_Center = CombinedRTs.Phos.Center.(SOALabel) - avgRT_NoPhos_Center;
adjustedRT_Phos_Other = CombinedRTs.Phos.Other.(SOALabel) - avgRT_NoPhos_Other;

% Perform the ranksum test
[p_value,~,~] = ranksum(adjustedRT_Phos_Center, adjustedRT_Phos_Other);

% Display the result
fprintf('At SOA %d ms, the p-value for the comparison between Center and Other is %f\n', SOA, p_value);
disp(stats)
%% Is there a Phos NoPhos dirrfence for other at 500 SOA
ranksum(CombinedRTs.Phos.Other.SOA_500, CombinedRTs.NoPhos.Other.SOA_500)

disp(stats)
%% Stats

%% rank sum to see if center diffrenece diffrent at an SOA

% Desired SOA value
targetSOA = 250;

% Retrieve the reaction times for the desired SOA
rt_NoPhos_Center_avg = mean(CombinedRTs.NoPhos.Center.(sprintf('SOA_%d', targetSOA)));
rt_NoPhos_Other_avg = mean(CombinedRTs.NoPhos.Other.(sprintf('SOA_%d', targetSOA)));

rt_Phos_Center_adjusted = CombinedRTs.Phos.Center.(sprintf('SOA_%d', targetSOA)) - rt_NoPhos_Center_avg;
rt_Phos_Other_adjusted = CombinedRTs.Phos.Other.(sprintf('SOA_%d', targetSOA)) - rt_NoPhos_Other_avg;

% Run rank sum test to compare adjusted Phos values between Center and Other
[p, ~] = ranksum(rt_Phos_Center_adjusted, rt_Phos_Other_adjusted);

% Display p-value
disp(['The p-value comparing adjusted Phos values for Center and Other at SOA = ' num2str(targetSOA) ' is: ' num2str(p)]);


%% Runs Ranksum test at SOA and Location

% Specify the desired values directly in the script
Loc_ =  7;
SOA  = 500;

% Convert the location and SOA to string format for field referencing
locationField = sprintf('Loc_%d', Loc_);
soaField = sprintf('SOA_%d', SOA);

% Extract data for the specified location and SOA
data_Phos = TrialRTs.Phos.(locationField).(soaField);
data_NoPhos = TrialRTs.NoPhos.(locationField).(soaField);

% Rank-sum test
[p,h,stats] = ranksum(data_Phos, data_NoPhos);

% Check and display significance
disp(['p=',num2str(p)])
disp(stats)

%% Get z statistic
% Perform the rank-sum test
[p,h,stats] = ranksum(data_Phos, data_NoPhos);

% The W statistic is already provided in stats
W = stats.ranksum;
disp(['W statistic: ', num2str(W)]);

% Calculate the sample sizes
n1 = numel(data_Phos);
n2 = numel(data_NoPhos);

% Calculate the expected mean and standard deviation of W under H0
meanW = n1*(n1+n2+1)/2;
stdW = sqrt(n1*n2*(n1+n2+1)/12);

% Calculate the Z statistic using the normal approximation
Z = (W - meanW) / stdW;
disp(['Z statistic: ', num2str(Z)]);

% Display the p-value
disp(['p-value: ', num2str(p)]);

% If stats also includes a zval, display it for comparison
if isfield(stats, 'zval')
    disp(['Z statistic from ranksum: ', num2str(stats.zval)]);
end


%% calcualte U

% Perform the rank-sum test
[p,h,stats] = ranksum(data_Phos, data_NoPhos);

% The W statistic is already provided in stats
W = stats.ranksum;

% Calculate the sample sizes for both groups
n1 = numel(data_Phos); % Size of the first sample (smaller group)
n2 = numel(data_NoPhos); % Size of the second sample (larger group)

% Calculate U1 and U2
U1 = W - (n1*(n1+1))/2;
U2 = n1*n2 - U1;

% Report the smaller U value
U = min(U1, U2);

% Display the results
disp(['U statistic: ', num2str(U)]);
disp(['p-value: ', num2str(p)]);


%% Diffrence regression

    % restructure data in to TrialRTs_Loc_Diff

% Initialize the new structure
TrialRTs_Loc_Diff = struct();

uniqueTargets = unique(allCondID)

LocationCount = numel(uniqueTargets)

% Loop through each SOA
for soa = uniqueSOATimes % assuming uniqueSOATimes is an array of all SOA values
    soaField = sprintf('SOA_%d', soa);
    TrialRTs_Loc_Diff.(soaField) = struct();
    
    % Loop through each location
    for loc = 1:LocationCount 
        locField = sprintf('Loc_%d', loc);
        
        % Calculate the mean NoPhos value for this SOA and location
        meanNoPhos = mean(TrialRTs_Location.NoPhos.(soaField).(locField));
        
        % Subtract this mean from each Phos value
        adjustedPhos = TrialRTs_Location.Phos.(soaField).(locField) - meanNoPhos;
        
        % Store in the new structure
        TrialRTs_Loc_Diff.(soaField).(locField) = adjustedPhos;
    end
end



%% regression absolute


% Define SOA and Axis
SOA = 500;  % Can be changed
axis = 'meridian';  % Can be 'meridian' or 'eccentricity'


% Define Meridian and Eccentricity Index locations
MeridianIndex = [4 5 6 7 8 9 10];
EccentricityIndex = [1 2 3 7 11 12 13];

% Calculate distances
distanceArc = linspace(-lengthArc / 2, lengthArc / 2, length(MeridianIndex));
distanceRadial = linspace(-lengthRadial / 2, lengthRadial / 2, length(EccentricityIndex));

% Select the appropriate index and distance based on axis
if strcmp(axis, 'meridian')
    Index = MeridianIndex;
    distances = abs(distanceArc);
else  % 'eccentricity'
    strcmp(axis, 'eccentricity')
    Index = EccentricityIndex;
    distances = abs(distanceRadial);
end

% Initialize arrays for regression data
RT_diff = [];
dist = [];

% Extract data for regression
SOALabel = sprintf('SOA_%d', SOA);
for i = 1:length(Index)
    locField = sprintf('Loc_%d', Index(i));
    RT_diff = [RT_diff; TrialRTs_Loc_Diff.(SOALabel).(locField)'];
    dist = [dist; repmat(distances(i), numel(TrialRTs_Loc_Diff.(SOALabel).(locField)), 1)];
end

% Run linear regression
lm = fitlm(dist, RT_diff);

% Display the results
disp(lm);


figure;
scatter(dist, RT_diff)


%% regression not absolute


% Define SOA and Axis
SOA = 500;  % Can be changed
axis = 'meridian';  % Can be 'meridian' or 'eccentricity'


% Define Meridian and Eccentricity Index locations
MeridianIndex = [4 5 6 7 8 9 10];
EccentricityIndex = [1 2 3 7 11 12 13];

% Calculate distances
distanceArc = linspace(-lengthArc / 2, lengthArc / 2, length(MeridianIndex));
distanceRadial = linspace(-lengthRadial / 2, lengthRadial / 2, length(EccentricityIndex));

% Select the appropriate index and distance based on axis
if strcmp(axis, 'meridian')
    Index = MeridianIndex;
    distances = distanceArc;
else  % 'eccentricity'
    strcmp(axis, 'eccentricity')
    Index = EccentricityIndex;
    distances = distanceRadial;
end

% Initialize arrays for regression data
RT_diff = [];
dist = [];

% Extract data for regression
SOALabel = sprintf('SOA_%d', SOA);
for i = 1:length(Index)
    locField = sprintf('Loc_%d', Index(i));
    RT_diff = [RT_diff; TrialRTs_Loc_Diff.(SOALabel).(locField)'];
    dist = [dist; repmat(distances(i), numel(TrialRTs_Loc_Diff.(SOALabel).(locField)), 1)];
end

% Run linear regression
lm = fitlm(dist, RT_diff);

% Display the results
disp(lm);


figure;
scatter(dist, RT_diff)


%% Center virsus other per direction

centRTdirection = RT_diff(dist == 0);
otherRTdirection = RT_diff(dist ~= 0);


% Get z statistic
% Perform the rank-sum test
[p,h,stats] = ranksum(centRTdirection,otherRTdirection);

% The W statistic is already provided in stats
W = stats.ranksum;
disp(['W statistic: ', num2str(W)]);

% Calculate the sample sizes
n1 = numel(data_Phos);
n2 = numel(data_NoPhos);

% Calculate the expected mean and standard deviation of W under H0
meanW = n1*(n1+n2+1)/2;
stdW = sqrt(n1*n2*(n1+n2+1)/12);

% Calculate the Z statistic using the normal approximation
Z = (W - meanW) / stdW;
disp(['Z statistic: ', num2str(Z)]);

% Display the p-value
disp(['p-value: ', num2str(p)]);

% If stats also includes a zval, display it for comparison
if isfield(stats, 'zval')
    disp(['Z statistic from ranksum: ', num2str(stats.zval)]);
end
%% Linearity Assumption


x= dist;
y=RT_diff;

% Assuming 'x' is your independent variable and 'y' is your dependent variable
lm = fitlm(x, y);  % Linear model

% Plotting residuals against predicted values
figure;
plot(lm.Fitted, lm.Residuals.Raw, '.');
xlabel('Fitted values');
ylabel('Residuals');
title('Residuals vs Fitted');

%% Independence Assumption

% Extract residuals and predictor variables from the LinearModel object
residuals = lm.Residuals.Raw;
predictors = table2array(lm.Variables(:, 2:end));  % Convert table to array

% Perform the Durbin-Watson test
[pValue, dwStat] = dwtest(residuals, predictors);
disp(['Durbin-Watson stat: ', num2str(dwStat)]);
disp(['p-value: ', num2str(pValue)]);


%% Homnoscidastisity Assumption

% Plot for homoscedasticity
figure;
plot(lm.Fitted, lm.Residuals.Raw, 'bo');
xlabel('Fitted values');
ylabel('Residuals');
title('Homoscedasticity Check');

%% Normality of residuals assumpion

% Q-Q plot
figure;
qqplot(lm.Residuals.Raw);
title('Q-Q Plot for Residuals');


% Assuming 'lm' is your LinearModel object
residuals = lm.Residuals.Raw;

% Shapiro-Wilk test for normality
% Note: swtest may need to be downloaded from MATLAB File Exchange
[pValue, W] = swtest(residuals);

% Display the results
disp(['Shapiro-Wilk W statistic: ', num2str(W)]);
disp(['p-value for normality: ', num2str(pValue)]);


%% Multicollinearity Check (For Multiple Regression) Assumption

% Calculating VIF - for multiple regression
% Assuming 'X' is your design matrix
[~,~,~,~,stats] = regress(y, X);
vif = 1./(1 - stats(1:end-1));
disp('VIF values:');
disp(vif);
