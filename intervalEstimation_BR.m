function data = intervalEstimation_BR(modality, interval_type,test_subject, n_trials, experiment, save_path, save_flag)

%% Help
% modality='auditory' or 'visual';  defines what modality the task will
% data = intervalEstimation_BR('auditory', 'long','bb', 3,'/Users/baylorbrangers/Desktop', 1)
% run.
% experiment = 1 or 2 - 1 means auditory first; 2 means visual first
% data = intervalEstimation_BR('visual','long','rod_test_1',100,1,'C:\Users\Rodrigo\Documents\INDP2015\Project\DATA',1)
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

%% Time intervals definitions
pre_stim_dist = 1.000; % Let us use a gaussion distribution around the mean

short_time_def = .45:.05 :.90; %seconds
long_time_def = .80:0.05: 1.250; %seconds
cumprob = 0:.1:1;
if islong
    curr_time_dist =  long_time_def;
elseif ~islong
    curr_time_dist =  short_time_def;
end

%% variable initialization
data.estimate = zeros(n_trials,1);
data.pre_stim = zeros(n_trials,1);
data.time = zeros(n_trials,1);
data.trial_time = zeros(n_trials,1);
data.correct=false(n_trials, 1);
data.abs_err = zeros(n_trials,1);
data.stim_presentation_time = zeros(n_trials,1);

score=0;

%% PsychToolbox initializations

PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');
% Draw to the external screen if avaliable
screenNumber = max(screens);
% screenNumber = 1;

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

%% Instructions
if isvisual
    line1 = 'Hello Subject';
    line2 = '\n\n\n Welcome to our Time Interval Estimation Task - Visual modality. ';
    line3 = '\n\n\n A white fixation cross will be displayed in the center of the screen.';
    line4 = '\n\n\n You will start each trial by pressing any key.';
    line5 = '\n\n\n After a random time interval, the Ready Stimulus (blue square) will be presented.';
    line6 = '\n\n\n After another random time interval the Set Stimulus (yellow square) will be shown.';
    line7 = '\n\n\n Your task is to replicate the time interval between the Ready and the Set Stimulus by';
    line8 = '\n\n pressing any key after what you think was the time that passed.';
    line9 = '\n\n\n You will be provided feedback on wether your estimation was close enough to the actual time interval.';
    line10 = '\n\n\n A white circle means you were correct, whilst a red circle means you were not good.';
    line11 = '\n\n The diameter of the circle is proportional to the absolute time difference of your estimate.';
    
    Endline = '\n\n\n\n Press any key to begin the session. Press again to start each trial.';
    Screen('TextSize', window, 20 );
    DrawFormattedText(window, [line1 line2 line3 line4 line5 line6 line7 line8 line9 line10 line11 Endline],...
        'center', screenYpixels * 0.1, white);
    Screen('Flip', window);
    
else  %If Auditory
    
    % Initialize Sounddriver
    InitializePsychSound(1);
    % Number of channels and Frequency of the sound
    nrchannels = 2;
    freq = 44000;
    % How many times to we wish to play the sound
    repetitions = 1;
    % Length of the beep
    beepLengthSecs = 0.1;
    % Length of the pause between beeps
    beepPauseTime = 1;
    % Start immediately (0 = immediately)
    startCue = 0;
    % Should we wait for the device to really start (1 = yes)
    waitForDeviceStart = 1;
    pahandle = PsychPortAudio('Open', [], 1, 1, freq, nrchannels);
    
    PsychPortAudio('Volume', pahandle, 0.5);
    %diff tones
    
    startBeep = MakeBeep(450, beepLengthSecs, freq);
    cueBeep = MakeBeep(600,beepLengthSecs, freq);
    correctBeep = MakeBeep(750,beepLengthSecs, freq);
    incorrectBeep = MakeBeep(300,beepLengthSecs, freq);
    
    % Instructions
    line1 = 'Hello Subject';
    line2 = '\n\n\n Welcome to our Time Interval Estimation Task - Auditory modality. ';
    line3 = '\n\n\n A white fixation cross will be displayed in the center of the screen.';
    line4 = '\n\n\n You will start each trial by pressing any key.';
    line5 = '\n\n\n After a random time interval, the Ready Stimulus (low frequency tone) will be played.';
    line6 = '\n\n\n After another random time interval the Set Stimulus (higher frequency tone) will be shown.';
    line7 = '\n\n\n Your task is to replicate the time interval between the Ready and the Set Stimulus by ';
    line8 = '\n\n pressing any key after what you think was the time that passed.';
    line9 = '\n\n\n You will be provided feedback on wether your estimation was close enough to the actual time interval.';
    line10 = '\n\n\n A high pitched sound means you were very close; ';
    line11 = '\n\n A low pitch means you were not good. (press any key to hear both)';
    
    Screen('TextSize', window, 20 );
    DrawFormattedText(window, [line1 line2 line3 line4 line5 line6 line7 line8 line9 line10 line11],...
        'center', screenYpixels * 0.10, white);
    Screen('Flip', window);
    
    KbWait;
    
    PsychPortAudio('FillBuffer', pahandle, [correctBeep; correctBeep]);
    PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
    
    WaitSecs(1);
    
    PsychPortAudio('FillBuffer', pahandle, [incorrectBeep; incorrectBeep]);
    PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
    
    WaitSecs(1);
    
    line12 = '\n\n\n You will also be given visual feedback.';
    line13 = '\n\n\n A white circle means you were correct, whilst a red circle means you were not good.';
    line14 = '\n\n The diameter of the circle is proportional to the absolute time difference of your estimate.';
    Endline = '\n\n\n\n Press any key to begin the session. Press again to start each trial.';
    Screen('TextSize', window, 25 );
    DrawFormattedText(window, [line12 line13 line14 Endline],...
        'center', screenYpixels * 0.25, white);
    Screen('Flip', window);
end
KbWait;

%% Trial starts
KbQueueCreate; % initialize the queue to get accurate timings

for trl = 1:n_trials
    stim_presentation_time = GetSecs;
    curr_pre_stim = pre_stim_dist + .100*randn(1);
    
    rand_val = rand;
    for i=1:10
        if rand_val>cumprob(i) &&  rand_val <=cumprob(i+1)
            current_value=i;
        end
    end   
    curr_time = curr_time_dist(current_value);

    error_threshold=(0.1*curr_time);  %10% of the stimulus interval
        
    Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
    Screen('Flip', window);
    
    KbQueueWait;
    
    %% Auditory Block
    if ~isvisual
        %wait from trial start
        WaitSecs(curr_pre_stim);
        KbQueueFlush;  

        PsychPortAudio('FillBuffer', pahandle, [startBeep; startBeep]);
        %%%%play ready sound
        PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
        %Wait between stimuli
        
        WaitSecs(curr_time);
        PsychPortAudio('FillBuffer', pahandle, [cueBeep; cueBeep]);
        %%%%%play set sound
        PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
        
        %Time reference
        stim_presentation_time=GetSecs;
        KbQueueStart;
        
        %wait for Keystroke
        curr_time_dif=KbQueueWait;
        
        curr_estimate = curr_time_dif - stim_presentation_time;
        curr_abs_err = abs(curr_estimate-curr_time);

        if curr_abs_err > error_threshold
            PsychPortAudio('FillBuffer', pahandle, [incorrectBeep; incorrectBeep]);
            PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
            dotColor = [1 0 0];  % WRONG > 10% error
            if curr_abs_err <=(0.2*curr_time)
                score = score + 0.5;
            end
        else
            PsychPortAudio('FillBuffer', pahandle, [correctBeep; correctBeep]);
            PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
            dotColor = [1 1 1];  % CORRECT     
            data.correct(trl)= true;
            score = score + 1;
        end
        
        scale_factor=500;
        baseRect_err = [xCenter-curr_abs_err*scale_factor yCenter-curr_abs_err*scale_factor xCenter+curr_abs_err*scale_factor yCenter+curr_abs_err*scale_factor];
        Screen('FillOval', window, dotColor, baseRect_err );
        Screen('Flip', window);        
        
        %%%%%play feedback sound
        WaitSecs(0.75);   
    
    %% Visual Block
    elseif isvisual
        %wait from trial start
        WaitSecs(curr_pre_stim);
        KbQueueFlush;  

        rectColor = [0 0 1];
        baseRect = [xCenter-100 yCenter-100 xCenter+100 yCenter+100];
        centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);
        
        
        Screen('FillRect', window, rectColor, centeredRect);
        
        Screen('Flip', window);
        
        %Wait between stimuli
        WaitSecs(curr_time);
        
        rectColor = [1 1 0];
        baseRect = [xCenter-100 yCenter-100 xCenter+100 yCenter+100];
        centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);
        
        Screen('FillRect', window, rectColor, centeredRect);
        Screen('Flip', window);
        
        %Time reference
        stim_presentation_time=GetSecs;
        KbQueueStart;
        
        %wait for Keystroke
        curr_time_dif=KbQueueWait;
        
        curr_estimate = curr_time_dif - stim_presentation_time;
        curr_abs_err = abs(curr_estimate-curr_time);
        
        
        if curr_abs_err > error_threshold
            dotColor = [1 0 0];  % WRONG - 10% error
            if curr_abs_err <=(0.2*curr_time)
                score = score + 0.5;
            end
        else
            dotColor = [1 1 1];  % CORRECT
            data.correct(trl)= true;
            score = score + 1;
        end
        
        
        scale_factor=500;
        
%         [ curr_abs_err    curr_abs_err*scale_factor]
        
        baseRect_err = [xCenter-curr_abs_err*scale_factor yCenter-curr_abs_err*scale_factor xCenter+curr_abs_err*scale_factor yCenter+curr_abs_err*scale_factor];

        Screen('FillOval', window, dotColor, baseRect_err );
        Screen('Flip', window);
        
        WaitSecs(0.75);
    end
    
    %% allocating relevant data to Structure
    
    data.pre_stim(trl) = curr_pre_stim;
    data.trial_time(trl) = curr_time;
    data.estimate(trl) = curr_estimate;
    data.abs_err(trl)=curr_abs_err;
    data.stim_presentation_time(trl) = stim_presentation_time;
    
    %% Trial end
end
KbQueueRelease; %Clear out the queue stuff
if ~isvisual
    % Close the audio device
    PsychPortAudio('Close', pahandle);
end
%% PsychtToolbox - closing instances and clearing buffers

Screen('CloseAll');

data.time_dist = curr_time_dist;
% Informations about session:
data.info.subject = test_subject;
data.info.modality = modality_text;
data.info.duration = interval_type;
data.info.n_trials = n_trials;
data.info.experiment = experiment;
data.info.score = score;

fprintf(1,['Congratulations! Your final score was of %f points!'],score)
%% Saving Data
if save_flag
    
    if ~exist([save_path filesep modality_text], 'dir')
        mkdir([save_path filesep modality_text]);
    end
    
    save([save_path filesep modality_text filesep modality_text '_' interval_type '_' test_subject '.mat'],'data')
end

end