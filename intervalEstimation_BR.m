function data = intervalEstimation_BR(modality, interval_type,test_subject, save_flag)
%% Help
% modality='auditory' or 'visual';  defines what modality the task will
% run.
%% Clear the screen and figures
close all;
sca;

%% Initializing a different seed for the random number generator every session
rng('shuffle');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Screen('Preference', 'SkipSyncTests', 1); %ONLY IF RUNNING IN SHITTY PC!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Session details are defined here
if strcmp(modality,'Visual') || strcmp(modality,'visual')
    isvisual =1;
    modality_text = 'visual';
elseif strcmp(modality,'Auditory') || strcmp(modality,'auditory')
    isvisual =0;
    modality_text = 'auditory';
else
    error('Modality unknown, the program only has Visual and Auditory modalities.')
end

if strcmp(interval_type,'Long') || strcmp(interval_type,'long')
    islong =1;
    interval_type_text = 'long';
elseif strcmp(interval_type,'Short') || strcmp(interval_type,'short')
    islong =0;
    interval_type_text = 'short';
else
    error('Time Interval Distribution unknown, the program only has Short and Long distributions.')
end

n_trials = 10;


%% Time intervals definitions
pre_stim_dist = 1.000; % Let us use a gaussion distribution around the mean

short_time_def = [.500 .850]; %seconds
long_time_def = [.850 1.200]; %seconds

if islong
   curr_time_dist =  long_time_def;
elseif ~islong
   curr_time_dist =  short_time_def;
end

%% variable initialization
data.estimate = zeros(n_trials,1);
data.pre_stim = zeros(n_trials,1);
data.time = zeros(n_trials,1);



%% PsychToolbox initializations

PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');
% Draw to the external screen if avaliable
% screenNumber = max(screens);
screenNumber = 1;

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Setup the text type for the window
Screen('TextFont', window, 'Ariel');
Screen('TextSize', window, 36);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Here we set the size of the arms of our fixation cross
fixCrossDimPix = 40;

% Now we set the coordinates (these are all relative to zero we will let
% the drawing routine center the cross in the center of our monitor for us)
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];

% Set the line width for our fixation cross
lineWidthPix = 4;

KbWait;

% Draw the fixation cross in white, set it to the center of our screen and
% set good quality antialiasing
Screen('DrawLines', window, allCoords,...
    lineWidthPix, white, [xCenter yCenter], 2);

% Flip to the screen
Screen('Flip', window);

%% Trial starts

for trl = 1:n_trials
    curr_pre_stim = pre_stim_dist + .100*randn(1);
    curr_time = curr_time_dist(1) + (curr_time_dist(2)-curr_time_dist(1))*rand(1);
 
    %% Auditory Block
    if ~isvisual
        
    end
    
    %% Visual Block
    if isvisual
        
        %wait from trial start
        WaitSecs(curr_pre_stim);
        
        
        
        dotColor = [1 0 0];
        dotXpos = rand * screenXpixels;
        dotYpos = rand * screenYpixels;
        dotSizePix = 20;
        Screen('DrawDots', window, [dotXpos dotYpos], dotSizePix, dotColor, [], 2);
        Screen('Flip', window);
        
        %Wait between stimuli
        WaitSecs(curr_time);
        
        dotColor2 = [1 1 0];
        Screen('DrawDots', window, [dotXpos dotYpos], dotSizePix, dotColor2, [], 2);
        Screen('Flip', window);
        
        %wait for Keystroke
        [a,b,c]=KbWait;
        curr_estimate = curr_time_dist(1) + (curr_time_dist(2)-curr_time_dist(1))*rand(1); %fake data!!!
        
        dotColor2 = [0 1 0];
        Screen('DrawDots', window, [dotXpos dotYpos], dotSizePix, dotColor2, [], 2);
        Screen('Flip', window);
    end
    
%% allocating relevant data to Structure

data.pre_stim(trl) = curr_pre_stim;
data.time(trl) = curr_time;
data.estimate(trl) = curr_estimate;
data.time_dist = curr_time_dist;


%% Trial end
end
%% PsychtToolbox - closing instances and clearing buffers

Screen('CloseAll');

%% Saving Data
if save_flag
    
    if ~exist([save_path filesep modality_text], 'dir')
        mkdir([save_path filesep modality_text]);
    end
    
    save([save_path filesep modality_text filesep modality_text '_' interval_type '_' test_subject '.mat'],'data')
end

end