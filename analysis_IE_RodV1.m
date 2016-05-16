function analysis = analysis_IE_RodV1(task_version, save_flag)

%% Loading data and grouping it
%Path to a folder with all data split in folders called "auditory","visual"...
%path_folder = 'C:\Users\Rodrigo\Documents\INDP2015\Project\DATA';
analysis_folder = 'C:\Users\Rodrigo\Documents\INDP2015\Project\DATA';
path_data = [analysis_folder filesep task_version filesep , '*.mat'];

% Gets data files to analyze
d = dir(path_data);
str = {d.name};
str = sortrows({d.name}');
[s,v] = listdlg('PromptString','Select files to group up:', 'OKString', 'OK',...
    'SelectionMode','multiple',...
    'ListString', str, 'Name', 'Select a File');
names = str(s);
numFiles = size(names, 1);

counter=0;
for i=1:numFiles
    file = load([analysis_folder filesep task_version filesep names{i}]);
    
    analysis.(file.data.info.duration).estimate(:,i) = file.data.estimate;
    analysis.(file.data.info.duration).pre_stim(:,i) = file.data.pre_stim;
    analysis.(file.data.info.duration).time(:,i) = file.data.time;
    analysis.(file.data.info.duration).trial_time(:,i) = file.data.trial_time;
    analysis.(file.data.info.duration).correct(:,i) = file.data.correct;
    analysis.(file.data.info.duration).abs_err(:,i) = file.data.abs_err;
    analysis.(file.data.info.duration).stim_presentation_time(:,i) = file.data.stim_presentation_time;
    analysis.(file.data.info.duration).time_dist = file.data.time_dist;
end
%INFO
analysis.info.modality = file.data.info.modality;
duration = {'short','long'}; %use this to cycle around 
n_trials = length(file.data.estimate);
%% Performance analysis

for i=1:numFiles %subjects
    for j=2 %durations
        analysis.(duration{j}).acc(1,i) = sum(analysis.(duration{j}).correct(:,i))/n_trials;
    end
end
%% Plotting

% Saving definitions


end