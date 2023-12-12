clear all
close all
clc

rootDir = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis'

% Load variables from the saved file
load('/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/allMartinData_20230824.mat');


%% trace results


samplerate = 1000; % 1000Hz

% Defining points to see when fixation broke and where


fixend = 50


% Find trials that have behavior outcomes 'SimPhosBreak' or 'BrokeBlankSOA'
falseCatchIndices = find(ismember(allBehaviorOutcomesSTR, {'SimPhosBreak', 'BrokeBlankSOA'}));

uniqueSOAs = unique(allSOATime(falseCatchIndices));
nSOAs = length(uniqueSOAs);

% Store the extracted traces per SOA
traceSegmentsX = cell(1, nSOAs);
traceSegmentsY = cell(1, nSOAs);

% Extract trace segments for the identified trials
for i = 1:length(falseCatchIndices)
    idx = falseCatchIndices(i);
    
    startTimeIdx = ((Times.PreStimT/1000-fixend/1000) * samplerate);
    endTimeIdx = (Times.PreStimT/1000 + allSOATime(idx)/1000 + Times.CatchFixT/1000) * samplerate;
    
    currentSOA = allSOATime(idx);
    soaIdx = find(uniqueSOAs == currentSOA);
    
    traceSegmentsX{soaIdx} = [traceSegmentsX{soaIdx}; allEyetrace.trials.x(idx, startTimeIdx:endTimeIdx)];
    traceSegmentsY{soaIdx} = [traceSegmentsY{soaIdx}; allEyetrace.trials.y(idx, startTimeIdx:endTimeIdx)];
end


%
% Zeroing the data based on the mean of the first fixend samples
for soaIdx = 1:nSOAs
    nTracesAtCurrentSOA = size(traceSegmentsX{soaIdx}, 1);
    
    for j = 1:nTracesAtCurrentSOA
        traceSegmentsX{soaIdx}(j,:) = traceSegmentsX{soaIdx}(j,:) - mean(traceSegmentsX{soaIdx}(j,1:fixend));
        traceSegmentsY{soaIdx}(j,:) = traceSegmentsY{soaIdx}(j,:) - mean(traceSegmentsY{soaIdx}(j,1:fixend));
    end
end

% remove fixend 
for soaIdx = 1:nSOAs
    traceSegmentsX{soaIdx} = traceSegmentsX{soaIdx}(:, fixend/1000*samplerate:end);
    traceSegmentsY{soaIdx} = traceSegmentsY{soaIdx}(:, fixend/1000*samplerate:end);
end



%%


% Compute the average trace for each SOA
averageTraceX = cellfun(@mean, traceSegmentsX, 'UniformOutput', false);
averageTraceY = cellfun(@mean, traceSegmentsY, 'UniformOutput', false);

% Plot the average traces for each SOA
figure;
hold on;
colors = jet(nSOAs); % Generate distinct colors for each SOA

for i = 1:nSOAs
    plot(averageTraceX{i}, averageTraceY{i}, 'Color', colors(i, :), 'DisplayName', ['SOA = ', num2str(uniqueSOAs(i)), 'ms']);
end

legend('Location', 'best');
title('Average Traces for Different SOAs');
xlabel('X');
ylabel('Y');
hold off;

%% plot all traces per SOA


% Set the axes limits
xlimits = [-3, 12];
ylimits = [-12, 3];

for i = 1:nSOAs
    figure;
    hold on;
    
    % Plot X vs Y for each trace at this SOA
    for j = 1:size(traceSegmentsX{i}, 1)
        plot(traceSegmentsX{i}(j, :), traceSegmentsY{i}(j, :));
    end
    
    % Add red dots at specified locations
    plot(0, 0, 'ro');
    plot(5, -5, 'ro');
    
    title(['SOA = ', num2str(uniqueSOAs(i)), 'ms']);
    xlabel('X Position');
    ylabel('Y Position');
    
    xlim(xlimits);
    ylim(ylimits);
    
    hold off;
end




%% distances

distancesFromOrigin = cell(1, nSOAs);

for soaIdx = 1:nSOAs
    nTracesAtCurrentSOA = size(traceSegmentsX{soaIdx}, 1);
    distancesFromOrigin{soaIdx} = zeros(nTracesAtCurrentSOA, size(traceSegmentsX{soaIdx}, 2));
    
    for j = 1:nTracesAtCurrentSOA
        distancesFromOrigin{soaIdx}(j,:) = sqrt(traceSegmentsX{soaIdx}(j,:).^2 + traceSegmentsY{soaIdx}(j,:).^2);
    end
end



%% when leaves fixation

fixationWindowRadius = 2;
firstExceedIndex = cell(1, nSOAs);

for soaIdx = 1:nSOAs
    nTracesAtCurrentSOA = size(distancesFromOrigin{soaIdx}, 1);
    firstExceedIndex{soaIdx} = zeros(nTracesAtCurrentSOA, 1);
    
    for j = 1:nTracesAtCurrentSOA
        exceedIndices = find(distancesFromOrigin{soaIdx}(j,:) > fixationWindowRadius, 1, 'first'); % finds the first index where distance > 2
        if ~isempty(exceedIndices)
            firstExceedIndex{soaIdx}(j) = exceedIndices;
        else
            firstExceedIndex{soaIdx}(j) = NaN; % If no such point exists, store NaN
        end
    end
end

figure;

for soaIdx = 1:nSOAs
    subplot(nSOAs, 1, soaIdx);
    histogram(firstExceedIndex{soaIdx}, 'BinWidth', 1); % You can adjust 'BinWidth' as necessary
    
    title(['SOA: ', num2str(uniqueSOAs(soaIdx))]);
    xlabel('First Exceed Index');
    ylabel('Frequency');
    xlim([0 max(cellfun(@max, firstExceedIndex))+5]); % Ensure consistent x limits across subplots
end


for soaIdx = 1:nSOAs
    firstExceedIndex{soaIdx} = firstExceedIndex{soaIdx}';
end


%% when trials enter target window

firstEnterIndex = cell(1, nSOAs);
targetCoord = [5, -5];
targetWindowRadius = 5;

for soaIdx = 1:nSOAs
    tracesX = traceSegmentsX{soaIdx};
    tracesY = traceSegmentsY{soaIdx};
    nTraces = size(tracesX, 1);
    firstEnterIndex{soaIdx} = zeros(1, nTraces);
    
    for traceIdx = 1:nTraces
        distanceFromTarget = sqrt((tracesX(traceIdx, :) - targetCoord(1)).^2 + (tracesY(traceIdx, :) - targetCoord(2)).^2);
        exceedIdx = find(distanceFromTarget < targetWindowRadius, 1);  % find the first time point
        
        if ~isempty(exceedIdx)
            firstEnterIndex{soaIdx}(traceIdx) = exceedIdx;
        else
            firstEnterIndex{soaIdx}(traceIdx) = NaN;  % If never enters within 5 units, set to NaN
        end
    end
end


%% Trials that didn't leave to early and didn't enter too late Lenient

mustLeave = Times.RespTooEarlyTW/1000*samplerate
mustEnter = (Times.RespTooEarlyTW + Times.RespTW)/1000*samplerate


% Initialize a cell array to store binary indices for each SOA
falseAlarmTrialsPerSOA = cell(1, nSOAs);

for soaIdx = 1:nSOAs
    exceedCondition = firstExceedIndex{soaIdx} > mustLeave;
    enterCondition = firstEnterIndex{soaIdx} < mustEnter;
    
    % Store the combined condition (logical AND) for this SOA
    falseAlarmTrialsPerSOA{soaIdx} = exceedCondition & enterCondition;
end



%% And must meet EyeFlyingTime restriction

mustLeave = Times.RespTooEarlyTW/1000 * samplerate;
mustEnter = (Times.RespTooEarlyTW + Times.RespTW) / 1000 * samplerate;
mustDiff = Times.EyeFlyingTime/1000 * samplerate; % Maximum allowable difference between enter and exceed times

% Initialize a cell array to store binary indices for each SOA
falseAlarmTrialsPerSOA = cell(1, nSOAs);

for soaIdx = 1:nSOAs
    exceedCondition = firstExceedIndex{soaIdx} > mustLeave;
    enterCondition = firstEnterIndex{soaIdx} < mustEnter;
    
    % Check if the difference between the enter and exceed times is less than or equal to mustDiff
    diffCondition = (firstEnterIndex{soaIdx} - firstExceedIndex{soaIdx}) <= mustDiff;
    
    % Store the combined conditions (logical AND across all) for this SOA
    falseAlarmTrialsPerSOA{soaIdx} = exceedCondition & enterCondition & diffCondition;
end


%% plot the false alarm trials


figure;

extraTime = 10/1000*samplerate; % Time in samples after enterCondition to extend the trace

% Set the axes limits
xlimits = [-3, 12];
ylimits = [-12, 3];

% Number of horizontal subplots in each row
nSubplotsPerRow = 3;

for soaIdx = 1:nSOAs
    % Calculate the correct subplot index for horizontal arrangement
    subplot(ceil(nSOAs/nSubplotsPerRow), nSubplotsPerRow, soaIdx);
    
    % Force create plotting area by using scatter
    scatter(0, 0, 50, 'r', 'o');  % This is the red filled dot at the origin
    hold on;
    scatter(5, -5, 50, 'b', 'o', 'filled'); % This is the blue filled dot at (5,-5)


% Red circle at (0,0) with a radius of 2
rectangle('Position', [-2,-2,4,4], 'Curvature', [1,1], 'EdgeColor', 'r', 'LineWidth', 1.5);

% Blue circle at (5,-5) with a radius of 5
rectangle('Position', [0,-10,10,10], 'Curvature', [1,1], 'EdgeColor', 'b', 'LineWidth', 1.5);


    
    % Retrieve the indices of the false alarm trials for the current SOA
    currentFalseAlarmIndices = find(falseAlarmTrialsPerSOA{soaIdx});
    
    % Loop through each false alarm trial and plot it up to the enterCondition time point plus extraTime
    if ~isempty(currentFalseAlarmIndices)
        for j = 1:length(currentFalseAlarmIndices)
            trialIdx = currentFalseAlarmIndices(j);
            
            % Determine the time point (index) to stop plotting based on enterCondition plus extraTime
            stopIdx = min(length(traceSegmentsX{soaIdx}(trialIdx, :)), firstEnterIndex{soaIdx}(trialIdx) + extraTime);
            
            plot(traceSegmentsX{soaIdx}(trialIdx, 1:stopIdx), traceSegmentsY{soaIdx}(trialIdx, 1:stopIdx));
            hold on;
        end
    end
    
    % Always set the limits and plot the points
    set(gca, 'XLim', xlimits);
    set(gca, 'YLim', ylimits);
    
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title(['Traces for SOA ' num2str(soaIdx)]);
end




%% Calculate false alarm counts

nSOAs = length(falseAlarmTrialsPerSOA);
falseAlarms = zeros(1, nSOAs);

for soaIdx = 1:nSOAs
    falseAlarms(soaIdx) = sum(falseAlarmTrialsPerSOA{soaIdx});
end


%% calculate signal theory trial result percentages

uniqueSOAs = unique(allSOATime(~isnan(allSOATime)));
nSOAs = length(uniqueSOAs);

% Preallocate arrays for storing percentages
CatchCorrectPercents = zeros(1, nSOAs);
CatchFalsePercents = zeros(1, nSOAs);
NoPhosHitPercents = zeros(1, nSOAs);
NoPhosMissPercents = zeros(1, nSOAs);
PhosHitPercents = zeros(1, nSOAs);
PhosMissPercents = zeros(1, nSOAs);

% Pre-compute masks
catchMasks = allFlagID == 0;
noPhosMasks = allFlagID == 1;
phosMasks = allFlagID == 2;

for i = 1:nSOAs
    currentSOAMask = allSOATime == uniqueSOAs(i);
    
    % Catch Trials
    currentCatchMask = catchMasks & currentSOAMask;
    CatchCorrect = sum(strcmp(allBehaviorOutcomesSTR(currentCatchMask), 'CorrectCatchResponse'));
    CatchFalse = falseAlarms(i);
    
    CatchTotal = CatchCorrect + CatchFalse;
    CatchCorrectPercents(i) = CatchCorrect / CatchTotal * 100;
    CatchFalsePercents(i) = CatchFalse / CatchTotal * 100;

    % No Phosphene Trials 
    currentNoPhosMask = noPhosMasks & currentSOAMask;
    NoPhosHit = sum(strcmp(allBehaviorOutcomesSTR(currentNoPhosMask), 'CorrectResponse'));
    NoPhosMiss = sum(ismember(allBehaviorOutcomesSTR(currentNoPhosMask), {'NoResponse', 'FlyingOvertime', 'WrongChoice', 'EyeCrossRespWindow'}));
    NoPhosTotal = NoPhosHit + NoPhosMiss;
    NoPhosHitPercents(i) = NoPhosHit / NoPhosTotal * 100;
    NoPhosMissPercents(i) = NoPhosMiss / NoPhosTotal * 100;

    % Phosphene Trials 
    currentPhosMask = phosMasks & currentSOAMask;
    PhosHit = sum(strcmp(allBehaviorOutcomesSTR(currentPhosMask), 'CorrectResponse'));
    PhosMiss = sum(ismember(allBehaviorOutcomesSTR(currentPhosMask), {'NoResponse', 'FlyingOvertime', 'WrongChoice', 'EyeCrossRespWindow'}));
    PhosTotal = PhosHit + PhosMiss;
    PhosHitPercents(i) = PhosHit / PhosTotal * 100;
    PhosMissPercents(i) = PhosMiss / PhosTotal * 100;
end

%% Dot plots

figure;

xValues = 1:length(uniqueSOAs); % Create equally spaced x values

% Flip the order of the percentages to match the flipped SOAs
CatchCorrectPercentsFlipped = fliplr(CatchCorrectPercents);
CatchFalsePercentsFlipped = fliplr(CatchFalsePercents);
NoPhosHitPercentsFlipped = fliplr(NoPhosHitPercents);
NoPhosMissPercentsFlipped = fliplr(NoPhosMissPercents);
PhosHitPercentsFlipped = fliplr(PhosHitPercents);
PhosMissPercentsFlipped = fliplr(PhosMissPercents);

% Catch Trials
subplot(1,3,3);
plot(xValues, CatchCorrectPercentsFlipped, '-og', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(xValues, CatchFalsePercentsFlipped, '-or', 'LineWidth', 1.5, 'MarkerSize', 8);
set(gca, 'XTick', xValues, 'XTickLabel', flip(uniqueSOAs));  % Set the tick labels to the reversed SOAs
legend('Correct Rejection', 'False Alarm', 'Location', 'northeast');
title('Catch Trials');
xlabel('SOAs');
ylabel('Percentage');
grid on;
ylim([0, 100]);
pbaspect([86 72 1]); % Set the aspect ratio to 4:3

% No Phosphene Trials
subplot(1,3,2);
plot(xValues, NoPhosHitPercentsFlipped, '-og', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(xValues, NoPhosMissPercentsFlipped, '-or', 'LineWidth', 1.5, 'MarkerSize', 8);
set(gca, 'XTick', xValues, 'XTickLabel', flip(uniqueSOAs));  % Set the tick labels to the reversed SOAs
legend('Hit', 'Miss', 'Location', 'northeast');
title('No Phosphene Trials');
xlabel('SOAs');
ylabel('Percentage');
grid on;
ylim([0, 100]);
pbaspect([86 72 1]); % Set the aspect ratio to 4:3

% Phosphene Trials
subplot(1,3,1);
plot(xValues, PhosHitPercentsFlipped, '-og', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(xValues, PhosMissPercentsFlipped, '-or', 'LineWidth', 1.5, 'MarkerSize', 8);
set(gca, 'XTick', xValues, 'XTickLabel', flip(uniqueSOAs));  % Set the tick labels to the reversed SOAs
legend('Hit', 'Miss', 'Location', 'northeast');
title('Phosphene Trials');
xlabel('SOAs');
ylabel('Percentage');
grid on;
ylim([0, 100]);
pbaspect([86 72 1]); % Set the aspect ratio to 4:3



