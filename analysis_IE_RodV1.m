function analysis = analysis_IE_RodV1(task_version, sufix,save_flag)
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
for j=2 %durations
    n_subjects(j)=size(analysis.(file.data.info.duration).stim_presentation_time,2);
end

%% Performance analysis
for i=1:numFiles %subjects
    for j=2 %durations
        analysis.(duration{j}).acc(1,i) = sum(analysis.(duration{j}).correct(:,i))/n_trials;
        analysis.(duration{j}).error(:,i) = analysis.(duration{j}).estimate(:,i) - analysis.(duration{j}).trial_time(:,i);
    end
end

%% Binning data to possible timings
for j=2
    for i=1:n_subjects(j)
        for k=1:length(analysis.(duration{j}).time_dist)
            analysis.(duration{j}).err_bin{i,k} = analysis.(duration{j}).estimate( analysis.(duration{j}).trial_time(:,i)==analysis.(file.data.info.duration).time_dist(k) ,i);
        end
    end  
end


%% Means
for j=2 %durations
    for i=1:n_subjects(j)
        for k=1:length(analysis.(duration{j}).time_dist)
            analysis.(duration{j}).error_m(i,k) = mean(analysis.(duration{j}).err_bin{i,k});
            analysis.(duration{j}).error_std(i,k) = std(analysis.(duration{j}).err_bin{i,k});
        end
    end
end

%Population
for j=2 %durations   
    analysis.(duration{j}).population.error_m(1,:) = mean(analysis.(duration{j}).error_m,1);
    analysis.(duration{j}).population.error_se(1,:) = std(analysis.(duration{j}).error_m)/sqrt(n_subjects(j));
end

%% Fitting 
for j=2 %durations
    for i=1:n_subjects(j)       
        analysis.(duration{j}).fit_params(i,:) = polyfit(analysis.(duration{j}).time_dist,analysis.(duration{j}).error_m(i,:),1);
    end
    
    analysis.(duration{j}).population.fit_params = mean(analysis.(duration{j}).fit_params,1);
    analysis.(duration{j}).population.fit_params_se = std(analysis.(duration{j}).fit_params,1)/sqrt(n_subjects(j));
end


%% Plotting

for j=2
    plot_max(j)=max(analysis.(duration{j}).time_dist)+0.05;
    plot_min(j)=min(analysis.(duration{j}).time_dist)-0.05;
end

% Saving definitions
save_folder = 'C:\Users\Rodrigo\Documents\INDP2015\Project\ANALYSIS';
figures_folder = [save_folder filesep 'Figures' filesep task_version];

if save_flag   
    if ~exist(save_folder, 'dir')
        mkdir(save_path);
    end
        
    if ~exist(figures_folder, 'dir')
        mkdir(figures_folder);
    end  
end

%% Plot estimate vs trial time
%Individuals
for j=2
    figure
    for i=1:n_subjects(j)
        subplot(2,round(n_subjects(j)/2),i)
        hold on
        errorbar(analysis.(duration{j}).time_dist, analysis.(duration{j}).error_m(i,:), analysis.(duration{j}).error_std(i,:),'ro')
        plot(analysis.(duration{j}).time_dist,polyval(analysis.(duration{j}).fit_params(i,:),analysis.(duration{j}).time_dist),'r-','LineWidth',1.5)
        plot([plot_min(j) plot_max(j)],[ plot_min(j) plot_max(j)],'--k','LineWidth',0.5)
        axis([plot_min(j) plot_max(j) plot_min(j) plot_max(j)])
        xlabel('Time Interval(s)')
        ylabel('Estimated Time(s)')
        set(gca,'XTick',plot_min(j):0.05:plot_max(j))
    end
    if save_flag
        saveas(gcf,[figures_folder filesep duration{j} '_EstVStrial',sufix],'png')
        saveas(gcf,[figures_folder filesep duration{j} '_EstVStrial',sufix],'fig')
    end
end

% Population
for j=2
    figure
    hold on
    errorbar(analysis.(duration{j}).time_dist,analysis.(duration{j}).population.error_m,analysis.(duration{j}).population.error_se,'o')
    plot([plot_min(j) plot_max(j)],[ plot_min(j) plot_max(j)],'--k','LineWidth',0.5)
    plot(analysis.(duration{j}).time_dist,polyval(analysis.(duration{j}).population.fit_params,analysis.(duration{j}).time_dist),'b-','LineWidth',1.5)
    axis([plot_min(j) plot_max(j) plot_min(j) plot_max(j)])
    xlabel('Time Interval(s)')
    ylabel('Estimated Time(s)')
    set(gca,'XTick',plot_min(j):0.05:plot_max(j))
    if save_flag
        saveas(gcf,[figures_folder filesep duration{j} '_EstVStrial_pop', sufix],'png')
        saveas(gcf,[figures_folder filesep duration{j} '_EstVStrial_pop', sufix],'fig')
    end
end


if save_flag   
    save([save_folder filesep task_version '_'  duration{j} '_analysis' sufix '.mat'],'analysis')
end
end