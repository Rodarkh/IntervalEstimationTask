function data = intervalEstimation_BR(modality, interval_type,test_subject, save_flag)
%% Help
% modality='auditory' or 'visual';  defines what modality the task will
% run.
%% Clear the workspace
clearvars;
close all;
sca;

%% Initializing a different seed for the random number generator every session
rng('shuffle');

%% Session details are defined here
if strcmp(modality,'Visual') || strcmp(modality,'visual')
    isvisual =1;
    modality_text = 'visual';
elseif strcmp(modality,'Auditory') || strcmp(modality,'auditory')
    isvisual =0;
    modality_text = 'auditory';
else
    error('Modality unknown, the program only has Visual and Auditory modalities')
end

if strcmp(interval_type,'Long') || strcmp(interval_type,'long')
    isvisual =1;
    interval_type_text = 'long';
elseif strcmp(interval_type,'Short') || strcmp(interval_type,'short')
    isvisual =0;
    interval_type_text = 'short';
else
    error('Modality unknown, the program only has Visual and Auditory modalities')
end

n_trials = 10;


%% Time intervals definitions

short_time_def = [500 850];
long_time_def = [850 1200];


%% PsychToolbox initializations




%% Trial starts

for trl = 1:n_trials
    curr
    
    
    %% Auditory Block
    if ~isvisual
        
    end
    
    %% Visual Block
    if isvisual
        
        
    end
    

%% PsychtToolbox - closing instances and clearing buffers




%% allocating relevant data to Structure
data.estimate = zeros(n_trials,1);
data.pre_stim = zeros(n_trials,1);
data.time = zeros(n_trials,1);

data.pre_stim(trl) = curr_pre_stim;
data.time(trl) = curr_time;
data.estimate(trl) = curr_estimate;
%% Trial end
end


%% Saving Data
if save_flag
    
    if ~exist([save_path filesep modality_text], 'dir')
        mkdir([save_path filesep modality_text]);
    end
    
    save([save_path filesep modality_text filesep modality_text '_' interval_type '_' test_subject '.mat'],'data')
end

end