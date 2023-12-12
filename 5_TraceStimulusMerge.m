

%==========================================================================
% SUMMARY:
%==========================================================================
% This code focuses on analyzing eye tracking data during an experimental task.
% 
% EXPERIMENT OVERVIEW:
% The task involves a monkey fixating on a visual cue and responding to various stimuli.
% There are three types of trials:
% 1. Catch Trial
% 2. No Phosphene Saccade Trial
% 3. Phosphene Saccade Trial
% 
% The trials are characterized by the presence of a photodiode signal, 
% specifically looking for events where the code changes a particular screen region's color 
% to white, signaling the monkey's response.
% 
% ERROR CODES:
% Each trial can result in various outcomes or 'error codes'. These error codes are associated 
% with the number of photodiode signals detected:

% zeroSignalCodes = {'PreFixBreak'};
% singleSignalCodes = {'PreStimBreak', 'NoPhosBreak', 'SimPhosBreak', 'BrokeBlankSOA'};
% doubleSignalCodes = {'CorrectResponse', 'CorrectCatchResponse', 'CatchBreak', 'RespondTooEarly', 'NoResponse', 'FlyingOvertime', 'WrongChoice', 'EyeCrossRespWindow'};

% STIMULUS VARIABLES

% TargetLocations:          A two-column matrix where each row corresponds to a unique target location (represented by x and y coordinates). 
% allCondID:                For each trial, an index indicating the row in 'TargetLocations' that corresponds to the target location used.
% allBehaviorOutcomes:      For each trial, a binary value indicating if the trial was successful (1) or not (0).
% allBehaviorOutcomesSTR:   For each trial, a string describing the specific outcome.
% allFixationTime:          For each trial, the duration of fixation prior to the stimulus.
% allResponseTime:          For each trial, two times indicating when in ms after the stimulus the eye left the fixation window (1st row) and when it entered the target window (2nd row).
% allFlagID:                For each trial, a flag indicating the type; 0 for catch trials, 1 for no phosphene saccade trials, and 2 for phosphene saccade trials.
% allSOATime     
% 
% SIGNAL PROCESSING:
% The photodiode data is processed to extract the onset timestamps of each trial. 
% This is achieved by thresholding the raw photodiode data, identifying onset and offset events,
% and segmenting the data into individual trials.
%
% Once the signal timestamps are extracted, the goal is to associate them with 
% the error codes from the trial. Two arrays are generated:
% TrialTimeStamps(1,:): Timestamps of the first photodiode signal for each trial
% TrialTimeStamps(2, nonNanIndices): Timestamps of the second photodiode signal for double-signal trials
% 
% OBJECTIVES:
% 1. SEGMENT EYE TRACES: Extract individual trial eye traces based on trial outcomes.
% 2. TRACE SCALING: Derive scaling factors from correct saccade trials and apply to all traces.
%==========================================================================


function [transformedTraces, stimulus] = processEyetrace(rootDir, traceFile, stimulusFile)


%% run this for use as not a function
% rootDir = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis'
% cd(rootDir);
% traceFile = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Varried_SOA/Eyetrace/Martin20230807_Block-1'
% stimulusFile = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data/Varried_SOA/Stimulus_Files/SimulatedPhosphene_20230807_164420.mat'
%%

% Define paths and SDK
TDTMatlabSDK = fullfile(rootDir, 'TDTSDK');
addpath(genpath(TDTMatlabSDK));


% Load stimulus file
stimulus = load(stimulusFile, 'allBehaviorOutcomes', 'allBehaviorOutcomesSTR', 'allCondID', 'allFixationTime', 'allResponseTime', 'allSOATime', 'allFlagID', 'TargetLocations', 'Times');

% Load trace file
trace = TDTbin2mat(traceFile);

channels = trace.streams.Eye_.data;
diode = channels(3, :);
samplerate = 1000; % Hz

%% Create a subplot for each channel
% 
% numRows = size(channels, 1); % Get the number of rows in the data
% 
% figure;
% for i = 1:numRows
%     subplot(numRows, 1, i);
%     plot(data.streams.Eye_.data(i, 1:50000));
%     title(['Channel ' num2str(i)]);
% end
%%
% Time stamps of diode events
% This code finds the onset timestamps of each trial and cuts eye data into trials
% using photodiode data from an experiment recording. It also removes the baseline 
% eye position from each trial. 

% KEY:
% - photodiodeData: raw data from the photodiode
% - threshold: threshold for identifying onset of each trial
% - idx: binary vector indicating onset (1) and offset (0) of photodiode signal below threshold
% - idx2: binary vector indicating where the signal drops below threshold
% - idx3: indices where the signal drops below threshold
% - idx4: differences between consecutive indices in idx3
% - idx5: binary vector indicating where the difference in idx4 is greater than 600
% - idx6: indices where the difference in idx4 is greater than 600
% - idx7: binary vector indicating where trials start (1) and end (0)


threshold = -3.3; % for the photo diode

idx = diode < threshold; % threshold crossing
% 
% %plot
% figure
% plot(diode((1:50000)))
% hold on
% plot(1:50000,idx((1:50000))*threshold,'o')
% 
% 
idx2 = diff(idx) == 1; % all threshold crossing time point
idx2 = [0,idx2];  

% %plot
% plot(1:50000,idx2((1:50000))*threshold,'x')

idx3 = find(idx2); 
idx3 = [0,idx3]; % all rising phase
idx4 = diff(idx3);

%Distance clusters are appart to be concidered independant trials (s)
secondsDistance = .05;
signalDistance = samplerate*secondsDistance;


idx5 = idx4 > signalDistance;
idx5 = logical([0,idx5]);

idx6 = idx3(idx5);
idx7 = zeros(size(idx2));
idx7(idx6) = threshold;

% %plot
% plot(1:50000,idx7((1:50000)),'d')

SignalOnsetTimestamps = find(idx7==threshold);

% figure
% hold on
% plot(diode)
% plot(SignalOnsetTimestamps,ones(size(SignalOnsetTimestamps))*threshold,'o','MarkerSize', 8,'LineWidth', 2)
% 
% 
% 


% Define the error codes and the number of signals associated with each

zeroSignalCodes = {'PreFixBreak'};
singleSignalCodes = {'PreStimBreak', 'NoPhosBreak', 'SimPhosBreak', 'BrokeBlankSOA'};
doubleSignalCodes = {'CorrectResponse', 'CorrectCatchResponse', 'CatchBreak', 'RespondTooEarly', 'NoResponse', 'FlyingOvertime', 'WrongChoice', 'EyeCrossRespWindow'};



% Filter out values in stimulus.allBehaviorOutcomesSTR that are in zeroSignalCodes
singleSignalBehaviorOutcomesSTR = stimulus.allBehaviorOutcomesSTR(~ismember(stimulus.allBehaviorOutcomesSTR, zeroSignalCodes));

% Get trials that have error codes associated with double signals
doubleSignalOnlyCodes = [doubleSignalCodes, singleSignalCodes]; % Include singleSignalCodes to remove them
doubleSignalBehaviorOutcomesSTR = stimulus.allBehaviorOutcomesSTR(~ismember(stimulus.allBehaviorOutcomesSTR, doubleSignalOnlyCodes));


% Get time stamps for PreStim And Stim

% Initialize TrialTimeStamps with empty columns
TrialTimeStamps = [];

% Keep track of our position in the SignalOnsetTimestamps
signalIndex = 1;

for i = 1:length(singleSignalBehaviorOutcomesSTR)
    % If the current behavior outcome is in singleSignalCodes
    if ismember(singleSignalBehaviorOutcomesSTR{i}, singleSignalCodes)
        % Add a column with the current timestamp on top and NaN at the bottom
        TrialTimeStamps = [TrialTimeStamps, [SignalOnsetTimestamps(signalIndex); NaN]];
        
        % Move to the next timestamp
        signalIndex = signalIndex + 1;
    % If the current behavior outcome is in doubleSignalCodes
    elseif ismember(singleSignalBehaviorOutcomesSTR{i}, doubleSignalCodes)
        % Add a column with the two consecutive timestamps
        TrialTimeStamps = [TrialTimeStamps, SignalOnsetTimestamps(signalIndex:signalIndex+1)'];
        
        % Move two positions ahead in the SignalOnsetTimestamps
        signalIndex = signalIndex + 2;
    end
end




% %TIMESTAMPS PLOTTINT
% figure
% hold on;
% 
% % Plot the diode data
% plot(diode);  % Assuming 'k' for black or choose any color you like
% 
% % Plot PreStimTimeStamps (blue circles)
% plot(TrialTimeStamps(1, :), ones(size(TrialTimeStamps(1, :))) * threshold, 'bo', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'PreStim');
% 
% % Plot StimTimeStamps (red circles)
% nonNanIndices = ~isnan(TrialTimeStamps(2, :));
% plot(TrialTimeStamps(2, nonNanIndices), ones(size(TrialTimeStamps(2, nonNanIndices))) * threshold, 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Stim');
% 
% xlabel('Time');
% ylabel('Diode Signal');
% title('PreStim and Stim TimeStamps over Diode Data');
% legend;
% hold off;


% cut eye data into Trials

% extract signals 
eyedata.x = channels(1, :);
eyedata.y = channels(2, :);

% Given data and assumptions
sampleRate = 1000;  % Assuming a sample rate of 1kHz; adjust as necessary
trialDuration = 2.5 * sampleRate;  % 2.5 seconds after PreStim onset

% Number of trials
nTrials = length(TrialTimeStamps(1,:));

% Preallocate memory for segmented data
eyedata.trial.x = NaN(nTrials, trialDuration);
eyedata.trial.y = NaN(nTrials, trialDuration);

% Loop through each trial and extract corresponding eye trace
for trialIdx = 1:nTrials
    % Convert the timestamp to an index
    onsetIdx = TrialTimeStamps(1, trialIdx);
    
    % Extract a segment from the onset of PreStim to 2 seconds after
    endIndex = onsetIdx + trialDuration-1;  % "-1" ensures that the duration is exactly 2 seconds
    
    eyedata.trial.x(trialIdx, :) = eyedata.x(onsetIdx:endIndex);
    eyedata.trial.y(trialIdx, :) = eyedata.y(onsetIdx:endIndex);
   
end

% defining section to zero data to 
baseDurationSamples = round(stimulus.Times.PreStimT/1000 * sampleRate);

% data zero'd 
for i = 1:nTrials
    eyedata.trial.x(i,:) = eyedata.trial.x(i,:) - mean(eyedata.trial.x(i,1:baseDurationSamples));
    eyedata.trial.y(i,:) = eyedata.trial.y(i,:) - mean(eyedata.trial.y(i,1:baseDurationSamples));
end


% extract only correct trials

% 1. Remove 'PreFixBreak' trials
validIdx = ~strcmp(stimulus.allBehaviorOutcomesSTR, 'PreFixBreak');

validBehaviorOutcomesSTR = stimulus.allBehaviorOutcomesSTR(validIdx);

% 2. Find indices of 'CorrectResponse' trials
correctIdx = strcmp(validBehaviorOutcomesSTR, 'CorrectResponse');

% 3. Extract relevant segments using the correct indices
eyedata.correct.x = eyedata.trial.x(correctIdx, :);
eyedata.correct.y = eyedata.trial.y(correctIdx, :);




validCondID = stimulus.allCondID(validIdx);
correctCondID = validCondID(correctIdx);



% The next thing is to cut down the correct trials to include just 100ms prior to stimulus and 100ms after eye entered the target window

% subtract first diode from stimulus onset
diodeSecondDelay = TrialTimeStamps(2,:) - TrialTimeStamps(1,:);

% find only correct delays
correctSecondDelay = diodeSecondDelay(correctIdx);

% Time before stimulus and affter saccade to cut trial
trialHeadTail = [.1,.1]; % seconds


% Syncing Response times
allSaccadeTime = stimulus.allResponseTime(2,:);
validSaccadeTime = allSaccadeTime(validIdx);
correctSaccadeTime = validSaccadeTime(correctIdx);

% Initialize vaiable with varying trial lengths
eyedata.correct.cut.x = nan(size(eyedata.correct.x));
eyedata.correct.cut.y = nan(size(eyedata.correct.y));

numCorrectTrials = sum(correctIdx);

for i = 1:numCorrectTrials
    startTrace = round(correctSecondDelay(i) - trialHeadTail(1) * samplerate);
    endTrace = round(correctSecondDelay(i) + correctSaccadeTime(i) + trialHeadTail(2) * samplerate);
    
    traceLength = endTrace - startTrace + 1;
    
    eyedata.correct.cut.x(i, 1:traceLength) = eyedata.correct.x(i, startTrace:endTrace);
    eyedata.correct.cut.y(i, 1:traceLength) = eyedata.correct.y(i, startTrace:endTrace);
end



% Plot first 10 individual traces

% figure;
% hold on;
% 
% % Plot the first 10 correct trials as an example
% for trialNum = 1:10
%     plot(eyedata.correct.cut.x(trialNum, :), eyedata.correct.cut.y(trialNum, :))
% end
% 
% title('Eye traces for first 10 correct trials')
% xlabel('X-coordinate')
% ylabel('Y-coordinate')
% legend(arrayfun(@(x) ['Trial ' num2str(x)], 1:10, 'UniformOutput', false))
% grid on;
% axis equal;
% hold off;

% FIND AVERAGE TRACES

numTargets = size(stimulus.TargetLocations, 1);

% Create storage for average traces
avgTraceX = nan(numTargets, size(eyedata.correct.cut.x, 2));
avgTraceY = nan(numTargets, size(eyedata.correct.cut.y, 2));


% Iterate through each condition
for condID = 1:numTargets
    % Find trials for the current condition
    trialsForCurrentCond = find(correctCondID == condID);
    
    % Compute the average trace for the current condition
    avgTraceX(condID, :) = nanmean(eyedata.correct.cut.x(trialsForCurrentCond, :), 1);
    avgTraceY(condID, :) = nanmean(eyedata.correct.cut.y(trialsForCurrentCond, :), 1);
    
end




% Plot average traces

% Define scaling factors
XscaleCon = 10;
YscaleCon = 10;

% Prepare the figure for plotting
close all;
figure;
hold on;


% Generate a colormap for the unique conditions
colors = lines(numTargets);

% Store handles for the traces (for the legend)
traceHandles = zeros(numTargets, 1);

% Iterate through each condition
for condID = 1:numTargets
    % Find trials for the current condition
    trialsForCurrentCond = find(correctCondID == condID);
    
    % Plot the average trace with scaling applied and store its handle
    traceHandles(condID) = plot(avgTraceX(condID, :) * XscaleCon, avgTraceY(condID, :) * YscaleCon, 'Color', colors(condID, :), 'LineWidth', 1.5);
    
    % Plot the target location for this condition using the same color
    plot(stimulus.TargetLocations(condID, 1), stimulus.TargetLocations(condID, 2), 'o', 'Color', colors(condID, :));
end


close all

% Enhancements for the plot
title('Average Eye Traces for Each Condition');
xlabel('Scaled X-coordinate');
ylabel('Scaled Y-coordinate');
xlim([-1 , 8 ]);  % Adjust as per your needs
ylim([-8, 1 ]);  % Adjust as per your needs
grid on;

% Add a legend to differentiate the traces, only using the trace handles
%legend(traceHandles, arrayfun(@(id) ['Cond ' num2str(id)], 1:numTargets, 'UniformOutput', false));

% EXTRACTING ENDPOINTS

% Get the number of targets
numTargets = size(stimulus.TargetLocations, 1);
samplerate = 1000;  % Assuming you've defined this earlier in your script

% Calculate the start and end indices based on the head and tail time
startSamples = round(trialHeadTail(1) * samplerate);
endSamples = round(trialHeadTail(1) * samplerate);

% Initialize the StartEnd structure
StartEnd.startCoor.x = nan(numTargets, startSamples);
StartEnd.startCoor.y = nan(numTargets, startSamples);
StartEnd.endCoor.x = nan(numTargets, endSamples);
StartEnd.endCoor.y = nan(numTargets, endSamples);

StartEnd.Average.startCoor.x = nan(numTargets, 1);
StartEnd.Average.startCoor.y = nan(numTargets, 1);
StartEnd.Average.endCoor.x = nan(numTargets, 1);
StartEnd.Average.endCoor.y = nan(numTargets, 1);

% Loop through each target
for i = 1:numTargets
    % Extract the beginning and end coordinates for each average trace
    StartEnd.startCoor.x(i, :) = avgTraceX(i, 1:startSamples);
    StartEnd.startCoor.y(i, :) = avgTraceY(i, 1:startSamples);
    
    % For the end coordinates, we'll take the last non-NaN segment
    validEndSegment = find(~isnan(avgTraceX(i, :)), endSamples, 'last');
    StartEnd.endCoor.x(i, :) = avgTraceX(i, validEndSegment-(endSamples-1):validEndSegment);
    StartEnd.endCoor.y(i, :) = avgTraceY(i, validEndSegment-(endSamples-1):validEndSegment);
    
    % Compute the average values for each section
    StartEnd.Average.startCoor.x(i) = nanmean(StartEnd.startCoor.x(i, :));
    StartEnd.Average.startCoor.y(i) = nanmean(StartEnd.startCoor.y(i, :));
    StartEnd.Average.endCoor.x(i) = nanmean(StartEnd.endCoor.x(i, :));
    StartEnd.Average.endCoor.y(i) = nanmean(StartEnd.endCoor.y(i, :));
end


% Plotting average points

figure;
hold on;

% Plot average start points
scatter(StartEnd.Average.startCoor.x, StartEnd.Average.startCoor.y, 'b', 'filled', 'DisplayName', 'Start Avg');
scatter(StartEnd.Average.endCoor.x, StartEnd.Average.endCoor.y, 'r', 'filled', 'DisplayName', 'End Avg');

% Adding indices to the end points
arrayfun(@(idx) text(StartEnd.Average.endCoor.x(idx), StartEnd.Average.endCoor.y(idx), ...
         sprintf('%d', idx), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom'), ...
         1:numTargets);

xlabel('X-coordinate');
ylabel('Y-coordinate');
title('Average Start and End Points for All Target Locations');
legend('show');
grid on;

close all

% Trans form points and Maintain indicies so that can connect the the stimulus location with calculated endpoint even after outliers removed

% 1. Subtract original start points from the endpoints
StartEnd.Average.endCoor.x = StartEnd.Average.endCoor.x - StartEnd.Average.startCoor.x;
StartEnd.Average.endCoor.y = StartEnd.Average.endCoor.y - StartEnd.Average.startCoor.y;

% 2. Zero out the start points
StartEnd.Average.startCoor.x = zeros(size(StartEnd.Average.startCoor.x));
StartEnd.Average.startCoor.y = zeros(size(StartEnd.Average.startCoor.y));

% 3. Prepare data for transformation
correctX = [stimulus.TargetLocations(:,1); zeros(length(StartEnd.Average.startCoor.x), 1)];
correctY = [stimulus.TargetLocations(:,2); zeros(length(StartEnd.Average.startCoor.y), 1)];

movingX = [StartEnd.Average.endCoor.x; StartEnd.Average.startCoor.x];
movingY = [StartEnd.Average.endCoor.y; StartEnd.Average.startCoor.y];

movingPoints = [movingX(:), movingY(:)];
fixedPoints = [correctX(:), correctY(:)];

% 4. Set RANSAC parameters and perform transformation
maxDistance = 10;
maxIterations = 1000;

[tform, inlierPoints] = estimateGeometricTransform(movingPoints, fixedPoints, 'similarity', ...
    'MaxDistance', maxDistance, 'MaxNumTrials', maxIterations);

transformedPoints = transformPointsForward(tform, movingPoints);
transformedX = transformedPoints(:, 1);
transformedY = transformedPoints(:, 2);



% Plot Lines Inbetween

% Plot the original and transformed points
figure;
h1 = plot(correctX, correctY, 'ro', 'MarkerSize', 8);
hold on;
h2 = plot(transformedX, transformedY, 'gx', 'MarkerSize', 8);

% Plot black lines between coordinate pairs
for i = 1:numel(correctX)
    plot([correctX(i), transformedX(i)], [correctY(i), transformedY(i)], 'b-');
end

% Set legend labels
legend([h1, h2], 'Stimulus Coordinates', 'Transformed Eyetrace EndPoints');

title('RANSAC Transformation, Outliers Removed ');

% Set the plot limits
xLimits = [-20, 20];
yLimits = [-20, 20];

close all


% % Transform eyedata.trial.x and eyedata.trial.y & avgTraceX and avgTraceY into transformedTraces

% Apply the transformation to all trials in eyedata.trial.x and eyedata.trial.y
[numTrials, numPointsPerTrial] = size(eyedata.trial.x);
transformedTraces.trials.x = NaN(numTrials, numPointsPerTrial);
transformedTraces.trials.y = NaN(numTrials, numPointsPerTrial);

for trialIdx = 1:numTrials
    trialPoints = [eyedata.trial.x(trialIdx, :)', eyedata.trial.y(trialIdx, :)'];
    transformedPoints = transformPointsForward(tform, trialPoints);
    transformedTraces.trials.x(trialIdx, :) = transformedPoints(:, 1)';
    transformedTraces.trials.y(trialIdx, :) = transformedPoints(:, 2)';
end

% Apply the transformation to avgTraceX and avgTraceY
[numTargets, numPointsPerTarget] = size(avgTraceX);
transformedTraces.average.x = NaN(numTargets, numPointsPerTarget);
transformedTraces.average.y = NaN(numTargets, numPointsPerTarget);

for targetIdx = 1:numTargets
    avgPoints = [avgTraceX(targetIdx, :)', avgTraceY(targetIdx, :)'];
    transformedAvgPoints = transformPointsForward(tform, avgPoints);
    transformedTraces.average.x(targetIdx, :) = transformedAvgPoints(:, 1)';
    transformedTraces.average.y(targetIdx, :) = transformedAvgPoints(:, 2)';
end

%Save transformed points
% Determine the split point (halfway)
splitIndex = numel(transformedX) / 2;

% Store the first half of the values into transformedTraces.endpoint
transformedTraces.endpoint.x = transformedX(1:splitIndex);
transformedTraces.endpoint.y = transformedY(1:splitIndex);

% Store the second half of the values into transformedTraces.startpoint
transformedTraces.startpoint.x = transformedX(splitIndex+1:end);
transformedTraces.startpoint.y = transformedY(splitIndex+1:end);



% Plot transformedTraces.average
% Prepare the figure for plotting
close all;
figure;
hold on;

% Generate a colormap for the unique conditions
colors = lines(numTargets);

% Store handles for the traces (for the legend)
traceHandles = zeros(numTargets, 1);

% Iterate through each condition
for condID = 1:numTargets
    % Plot the transformed average trace and store its handle
    traceHandles(condID) = plot(transformedTraces.average.x(condID, :), transformedTraces.average.y(condID, :), 'Color', colors(condID, :), 'LineWidth', 1.5);
    
    % Plot the target location for this condition using the same color
    plot(stimulus.TargetLocations(condID, 1), stimulus.TargetLocations(condID, 2), 'o', 'Color', colors(condID, :));
end

% Enhancements for the plot
title('Transformed Average Eye Traces for Each Condition');
xlabel('X-coordinate');
ylabel('Y-coordinate');
xlim([-1 , 8 ]);  % Adjust as per your needs
ylim([-8, 1 ]);  % Adjust as per your needs
grid on;

% Add a legend to differentiate the traces, only using the trace handles
legend(traceHandles, arrayfun(@(id) ['Cond ' num2str(id)], 1:numTargets, 'UniformOutput', false));

close all


% Plot Lines Inbetween and average traces

% Plot the original and transformed points
figure;
h1 = plot(correctX, correctY, 'ro', 'MarkerSize', 8);
hold on;
h2 = plot(transformedX, transformedY, 'gx', 'MarkerSize', 8);

% Plot black lines between coordinate pairs
for i = 1:numel(correctX)
    plot([correctX(i), transformedX(i)], [correctY(i), transformedY(i)], 'b-');
end

% Add the transformed average traces
colors = lines(numTargets); % Generate a colormap for the unique conditions

% Store handles for the traces (for the legend)
traceHandles = zeros(numTargets, 1);

% Iterate through each condition
for condID = 1:numTargets
    % Plot the transformed average trace and store its handle
    traceHandles(condID) = plot(transformedTraces.average.x(condID, :), transformedTraces.average.y(condID, :), 'Color', colors(condID, :), 'LineWidth', 1.5);
end

% Set legend labels
%legend([h1, h2, traceHandles], ['Stimulus Coordinates', 'Transformed Eyetrace EndPoints', arrayfun(@(id) ['Transformed Trace ' num2str(id)], 1:numTargets, 'UniformOutput', false)]);

title('RANSAC Transformation, Outliers Removed');

% Set the plot limits
xLimits = [-20, 20];
yLimits = [-20, 20];
xlim(xLimits);
ylim(yLimits);


% Assuming validIdx is either loaded with the stimulus file or defined somewhere above in this function
% Filter stimulus data for only valid trials
stimulus.allBehaviorOutcomes = stimulus.allBehaviorOutcomes(validIdx);
stimulus.allBehaviorOutcomesSTR = stimulus.allBehaviorOutcomesSTR(validIdx);
stimulus.allCondID = stimulus.allCondID(validIdx);
stimulus.allFixationTime = stimulus.allFixationTime(validIdx);
stimulus.allResponseTime = stimulus.allResponseTime(validIdx);
stimulus.allSOATime = stimulus.allSOATime(validIdx);
stimulus.allFlagID = stimulus.allFlagID(validIdx);

% adding more to treasformedTraces

transformedTraces.diodeSecondDelay = diodeSecondDelay;

close all

end