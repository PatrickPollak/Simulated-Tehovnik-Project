
% ---------------------
% VARIABLE DESCRIPTIONS
% ---------------------
%
% allBehaviorOutcomes:
%   A matrix containing outcomes of all the behavioral tasks.
%
% allBehaviorOutcomesSTR:
%   A matrix of strings describing the outcomes of all behavioral tasks.
%
% allCondID:
%   A matrix specifying the condition ID for each trial.
%
% allDiodeSecondDelay:
%   A vector indicating the delay of the diode for each trial.
%
% allEyetrace:
%   .trials.x: 
%       A matrix of X-coordinates for the eyetrace for each trial.
%   .trials.y:
%       A matrix of Y-coordinates for the eyetrace for each trial.
%
% allFixationTime:
%   A matrix indicating the time taken for fixation in each trial.
%
% allFlagID:
%   A matrix specifying the flag ID for each trial.
%
% allResponseTime:
%   A matrix indicating the response time for each trial.
%
% allSOATime:
%   A matrix indicating the stimulus onset asynchrony time for each trial.
%
% allStimulus:
%   A structure containing information related to stimulus for each trial.
%
% allTransformedTraces:
%   A structure containing transformed traces for each trial.
%
% Times:
%   A structure containing various time points of interest.
%
% TargetLocations:
%   A matrix specifying the X and Y target locations.



%%


clear all
close all;

% Define the root path
rootPath = '/Users/patrickpollak/Documents/Classes/Year2/Visual_Prosthesis/Martin_Data';

% Define paths to trace and stimulus files
traceFilesDir = fullfile(rootPath, 'Varried_SOA/Eyetrace/');
stimulusFilesDir = fullfile(rootPath, 'Varried_SOA/Stimulus_Files/*.mat');

% Get a list of the trace and stimulus files
traceFiles = dir(traceFilesDir);
stimulusFiles = dir(stimulusFilesDir);

% Filter out any non-folders from the traceFiles list
traceFiles = traceFiles([traceFiles.isdir]);

% Remove the "." and ".." entries
traceFiles = traceFiles(~ismember({traceFiles.name}, {'.', '..'}));

% Ensure the number of trace files and stimulus files are the same
if length(traceFiles) ~= length(stimulusFiles)
    error('Mismatch in number of trace and stimulus files');
end

% Iterate through each file and process them
allTransformedTraces = [];
allStimulus = [];

for i = 1:length(traceFiles)
    % Get full path to the current trace and stimulus file
    traceFilePath = fullfile(traceFiles(i).folder, traceFiles(i).name);
    stimulusFilePath = fullfile(stimulusFiles(i).folder, stimulusFiles(i).name);
    
    % Process the current pair of files
    [transformedTraces, stimulus] = TraceStimulusMerge(rootPath, traceFilePath, stimulusFilePath);
    
    % Concatenate or accumulate the results as necessary
    allTransformedTraces = [allTransformedTraces; transformedTraces];
    allStimulus = [allStimulus; stimulus];  % Adjust as necessary if `stimulus` is not a simple structure or array
end

%% merging into single matrices

% Initialize empty arrays for the variables to be concatenated
allBehaviorOutcomes = [];
allBehaviorOutcomesSTR = [];
allCondID = [];
allFixationTime = [];
allResponseTime = [];
allSOATime = [];
allFlagID = [];

% Iterate through the allStimulus structure array and concatenate the variables
for i = 1:length(allStimulus)
    allBehaviorOutcomes = [allBehaviorOutcomes, allStimulus(i).allBehaviorOutcomes];
    allBehaviorOutcomesSTR = [allBehaviorOutcomesSTR, allStimulus(i).allBehaviorOutcomesSTR];
    allCondID = [allCondID, allStimulus(i).allCondID];
    allFixationTime = [allFixationTime, allStimulus(i).allFixationTime];
    allResponseTime = [allResponseTime, allStimulus(i).allResponseTime];
    allSOATime = [allSOATime, allStimulus(i).allSOATime];
    allFlagID = [allFlagID, allStimulus(i).allFlagID];
end

% Extract a single instance of 'TargetLocations' and 'Times'
TargetLocations = allStimulus(1).TargetLocations;
Times = allStimulus(1).Times;

%% combine eye traces data

% Initialize empty arrays
allEyetrace.trials.x = [];
allEyetrace.trials.y = [];
allDiodeSecondDelay = [];

% Loop through allTransformedTraces and concatenate data
for i = 1:numel(allTransformedTraces)
    allEyetrace.trials.x = [allEyetrace.trials.x; allTransformedTraces(i).trials.x];
    allEyetrace.trials.y = [allEyetrace.trials.y; allTransformedTraces(i).trials.y];
    allDiodeSecondDelay = [allDiodeSecondDelay, allTransformedTraces(i).diodeSecondDelay];
end





%% parts to save 


% Get the current date
currentDate = datestr(now, 'yyyymmdd');

% File name to be saved
fileName = ['allMartinData_' currentDate '.mat'];

% Save all variables at once
save(fileName, 'allBehaviorOutcomes', 'allBehaviorOutcomesSTR', 'allCondID', ...
     'allDiodeSecondDelay', 'allEyetrace', 'allFixationTime', ...
     'allFlagID', 'allResponseTime', 'allSOATime', 'allStimulus', ...
     'allTransformedTraces', 'Times', 'TargetLocations');

