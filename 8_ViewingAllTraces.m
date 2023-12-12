clc;close all;clear all

%%
% Load variables from the saved file
load('/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/allMartinData_20230824.mat');

%%
plot(allTransformedTraces(3).average.x(1,:) , allTransformedTraces(3).average.y(1,:) )


%%



% Number of sessions
numSessions = length(allTransformedTraces);


EndLength = 400

% Initialize matrices to store the averages for x and y
averageX = zeros(13, EndLength);
averageY = zeros(13, EndLength);

% Calculate the average for each target and each time point
for targetIdx = 1:13
    for timePoint = 1:EndLength
        xSum = 0;
        ySum = 0;
        for sessionIdx = 1:numSessions
            xSum = xSum + allTransformedTraces(sessionIdx).average.x(targetIdx, timePoint);
            ySum = ySum + allTransformedTraces(sessionIdx).average.y(targetIdx, timePoint);
        end
        averageX(targetIdx, timePoint) = xSum / numSessions;
        averageY(targetIdx, timePoint) = ySum / numSessions;
    end
end

%




baseDurationSamples = 2

% data zero'd 
for i = 1:13
    averageX(i,:) = averageX(i,:) - mean(averageX(i,1:baseDurationSamples));
    averageY(i,:) = averageY(i,:) - mean(averageY(i,1:baseDurationSamples));
end




scale = .57
% Assuming 'TargetLocations' is a 13x2 matrix with target coordinates

% Plot the average traces for all targets
figure;
hold on; % Allows multiple plots in the same figure
for i = 1:13 
    plot(averageX(i, :)*scale, averageY(i, :)*scale);
    % Optionally, you can add markers for start and end points of each trace
    % plot(averageX(i, 1), averageY(i, 1), 'go'); % Start point in green
    % plot(averageX(i, end), averageY(i, end), 'ro'); % End point in red
end

% Plot target locations
plot(TargetLocations(:, 1), TargetLocations(:, 2), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
% 'kx' plots the points in black color with 'x' marker

title('Average Traces per Targets');
xlabel('Visual Degrees');
ylabel('Y Coordinate');
legend('Traces', 'Target Locations', 'Location', 'best');
hold off; % Release the figure for other plots
