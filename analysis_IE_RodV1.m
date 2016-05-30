function analysis = analysis_IE_RodV1(task_version, experiment, sufix,save_flag)
%% Loading data and grouping it
%Path to a folder with all data split in folders called "auditory","visual"...
%path_folder = 'C:\Users\Rodrigo\Documents\INDP2015\Project\DATA';

analysis_folder = 'C:\Users\Rodrigo\Documents\INDP2015\Project\DATA';
%  analysis_folder = '/Users/baylorbrangers/Desktop/Subject_Data_Bayes/Data';

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
    
    if file.data.info.experiment == experiment %Only load if its from the same experiments
        analysis.(file.data.info.duration).estimate(:,i) = file.data.estimate;
        analysis.(file.data.info.duration).pre_stim(:,i) = file.data.pre_stim;
        analysis.(file.data.info.duration).time(:,i) = file.data.time;
        analysis.(file.data.info.duration).trial_time(:,i) = file.data.trial_time;
        analysis.(file.data.info.duration).correct(:,i) = file.data.correct;
        analysis.(file.data.info.duration).abs_err(:,i) = file.data.abs_err;
        analysis.(file.data.info.duration).stim_presentation_time(:,i) = file.data.stim_presentation_time;
        analysis.(file.data.info.duration).time_dist = file.data.time_dist;
        analysis.(file.data.info.duration).score.(file.data.info.subject) = file.data.info.score;
        
    end
end

%INFO
analysis.info.modality = file.data.info.modality;
analysis.info.experiment = experiment;
duration = {'short','long'}; %use this to cycle around 
n_trials = length(file.data.estimate);
for j=2 %durations
    n_subjects(j)=size(analysis.(file.data.info.duration).stim_presentation_time,2);
end

%% Performance analysis
for j=2 %durations
    
    for i=1:numFiles %subjects
        
        analysis.(duration{j}).acc(1,i) = sum(analysis.(duration{j}).correct(:,i))/n_trials;
        analysis.(duration{j}).bias(:,i) = analysis.(duration{j}).estimate(:,i) - analysis.(duration{j}).trial_time(:,i);
   
        
        for k=1:length(analysis.(duration{j}).time_dist)
            trials_bins(i,k)= length(analysis.(duration{j}).correct( analysis.(duration{j}).trial_time(:,i)==analysis.(file.data.info.duration).time_dist(k) ,i));
            correct(i,k)= sum(analysis.(duration{j}).correct( analysis.(duration{j}).trial_time(:,i)==analysis.(file.data.info.duration).time_dist(k) ,i));
            analysis.(duration{j}).acc_bin(i,k) = correct(i,k) / trials_bins(i,k);

        end
            
        
    end
end

%% Pre stim analysis

bin_vector=[-inf 0.95 1.05 +inf];
[pre_stim_times pre_stim_bins]=histc(analysis.long.pre_stim(:,:),bin_vector);
figure
bar(pre_stim_times)

for j=2
    analysis.(duration{j}).pre_stim_b=false(n_trials,n_subjects(j),max(max(pre_stim_bins)));
    for bin=1:max(max(pre_stim_bins))
        for i=1:n_subjects(j)
            analysis.(duration{j}).pre_stim_b.bins(pre_stim_bins(:,i)==bin,i,bin) = true;
        end
    end
    
    analysis.(duration{j}).pre_stim_b.n_trials = sum(analysis.(duration{j}).pre_stim_b.bins);
    
    analysis.(duration{j}).pre_stim_b.acc = zeros(1,n_subjects(j),max(max(pre_stim_bins)));
    
    analysis.(duration{j}).pre_stim_b.bias= cell(max(max(pre_stim_bins)),n_subjects(j));
    
    for bin=1:max(max(pre_stim_bins))
        for i= 1:n_subjects(j)
            analysis.(duration{j}).pre_stim_b.acc(1,i,bin) = sum(analysis.(duration{j}).correct( analysis.(duration{j}).pre_stim_b.bins(:,i,bin),i ))/analysis.(duration{j}).pre_stim_b.n_trials(1,i,bin);
            analysis.(duration{j}).pre_stim_b.bias{bin,i} = analysis.(duration{j}).estimate(analysis.(duration{j}).pre_stim_b.bins(:,i,bin),i) - analysis.(duration{j}).trial_time(analysis.(duration{j}).pre_stim_b.bins(:,i,bin),i);
        end
    end
%     
%     for k=1:length(analysis.(duration{j}).time_dist)
%         trials_bins(i,k)= length(analysis.(duration{j}).correct( analysis.(duration{j}).trial_time(:,i)==analysis.(file.data.info.duration).time_dist(k) ,i));
%         correct(i,k)= sum(analysis.(duration{j}).correct( analysis.(duration{j}).trial_time(:,i)==analysis.(file.data.info.duration).time_dist(k) ,i));
%         analysis.(duration{j}).acc_bin(i,k) = correct(i,k) / trials_bins(i,k);
%         
%     end
    
    
end


%% Binning data to possible timings
for j=2
    for i=1:n_subjects(j)
        for k=1:length(analysis.(duration{j}).time_dist)
            analysis.(duration{j}).est_bin{i,k} = analysis.(duration{j}).estimate( analysis.(duration{j}).trial_time(:,i)==analysis.(file.data.info.duration).time_dist(k) ,i); 
        end       
    end  

end


%% Means
for j=2 %durations
    for i=1:n_subjects(j)
        for k=1:length(analysis.(duration{j}).time_dist)
            analysis.(duration{j}).est_m(i,k) = nanmean(analysis.(duration{j}).est_bin{i,k});
            analysis.(duration{j}).est_std(i,k) = nanstd(analysis.(duration{j}).est_bin{i,k});
        end
    end
    
    
    analysis.(duration{j}).biasRMS=sqrt(sum(analysis.(duration{j}).est_m(:,:).^2,2)/length(analysis.(duration{j}).time_dist));
    analysis.(duration{j}).VAR = nanmean(analysis.(duration{j}).est_std(:,:),2);
    
    
    %Population    
    analysis.(duration{j}).population.est_m(1,:) = nanmean(analysis.(duration{j}).est_m,1);
    analysis.(duration{j}).population.error_se(1,:) = nanstd(analysis.(duration{j}).est_m)/sqrt(n_subjects(j));
    
    analysis.(duration{j}).population.acc_bin_m = nanmean( analysis.(duration{j}).acc_bin,1);
    analysis.(duration{j}).population.acc_bin_se = nanstd( analysis.(duration{j}).acc_bin,1,1)/sqrt(n_subjects(j));
    
    
    analysis.(duration{j}).population.biasRMS_m= nanmean(analysis.(duration{j}).biasRMS) ;
    analysis.(duration{j}).population.biasRMS_se= nanstd(analysis.(duration{j}).biasRMS,1,1) / sqrt(n_subjects(j));
    
    analysis.(duration{j}).population.VAR_m= nanmean(analysis.(duration{j}).VAR) ;
    analysis.(duration{j}).population.VAR_se= nanstd(analysis.(duration{j}).VAR) / sqrt(n_subjects(j));
    
    
    
end


%% Fitting 
for j=2 %durations
    for i=1:n_subjects(j)       
        analysis.(duration{j}).fit_params(i,:) = polyfit(analysis.(duration{j}).time_dist,analysis.(duration{j}).est_m(i,:),1);
    end
    
    analysis.(duration{j}).population.fit_params = mean(analysis.(duration{j}).fit_params,1);
    analysis.(duration{j}).population.fit_params_se = std(analysis.(duration{j}).fit_params,1)/sqrt(n_subjects(j));
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sliding Window Analysis
for j=2 %duration
    sliding_window=20;
    for trl=1:n_trials-sliding_window
        
        
        for i=1:n_subjects(j) %subjects
            
            % Performance analysis
            analysis.(duration{j}).SW.acc(trl,i) = sum(analysis.(duration{j}).correct(trl:trl+sliding_window,i))/sliding_window;
            analysis.(duration{j}).SW.error(trl,i) = nanmean(analysis.(duration{j}).estimate(trl:trl+sliding_window,i) - analysis.(duration{j}).trial_time(trl:trl+sliding_window,i));
            analysis.(duration{j}).SW.trial_tms_beg(trl,i) = nanmean(analysis.(duration{j}).trial_time(1:trl+sliding_window,i));
            analysis.(duration{j}).SW.trial_tms_curr(trl,i) = nanmean(analysis.(duration{j}).trial_time(trl:trl+sliding_window,i));
            
            % Binning data to possible timings
            for k=1:length(analysis.(duration{j}).time_dist)
                analysis.(duration{j}).SW.est_bin{i,k,trl} = analysis.(duration{j}).estimate( analysis.(duration{j}).trial_time(trl:trl+sliding_window,i)==analysis.(file.data.info.duration).time_dist(k) ,i);
            end
            
        end
        
        % Means        
        for i=1:n_subjects(j)
            for k=1:length(analysis.(duration{j}).time_dist)
                analysis.(duration{j}).SW.est_m(i,k,trl) = nanmean(analysis.(duration{j}).SW.est_bin{i,k,trl});
                analysis.(duration{j}).SW.error_std(i,k,trl) = nanstd(analysis.(duration{j}).SW.est_bin{i,k,trl});
            end
        end
                  
    end
    
    
    
    for i=1:n_subjects(j)
        for trl=1:n_trials-sliding_window
            if sum( isnan(analysis.(duration{j}).SW.est_m(i,:,trl)) ) == 0
                analysis.(duration{j}).SW.fit_params(i,:,trl) = polyfit(analysis.(duration{j}).time_dist,analysis.(duration{j}).SW.est_m(i,:,trl),1);
            else
                non_nans = ~isnan(analysis.(duration{j}).SW.est_m(i,:,trl));
                x_values=analysis.(duration{j}).time_dist(non_nans);
                analysis.(duration{j}).SW.fit_params(i,:,trl) = polyfit(x_values,analysis.(duration{j}).SW.est_m(i,non_nans,trl),1);
            end
        end
    end
    
    %Population   
    analysis.(duration{j}).SW.pop_est_m(:,:) = squeeze(nanmean(analysis.(duration{j}).SW.est_m,1))';
    analysis.(duration{j}).SW.pop_error_se(:,:) = squeeze(nanstd(analysis.(duration{j}).SW.est_m)/sqrt(n_subjects(j)))';
    
    analysis.(duration{j}).SW.pop_fit_params = squeeze(nanmean(analysis.(duration{j}).SW.fit_params,1))';
    analysis.(duration{j}).SW.pop_fit_params_se = squeeze(nanstd(analysis.(duration{j}).SW.fit_params,1)/sqrt(n_subjects(j)))';
       
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting

for j=2
    plot_max(j)=max(analysis.(duration{j}).time_dist)+0.05;
    plot_min(j)=min(analysis.(duration{j}).time_dist)-0.05;
end

experiment_text = {'Aud to Vis'; 'Vis to Aud'};
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
    %Accuracy plot
    figure
    for i=1:n_subjects(j)
        subplot(2,round(n_subjects(j)/2),i)
        hold on
        plot(1:n_trials-sliding_window, analysis.(duration{j}).SW.acc(:,i),'g-')
        xlabel('Trial')
        ylabel('Accuracy')
        ylim([0 1])
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,[task_version '-' duration{j} ' - ' experiment_text{experiment} '-Accuracy'],'HorizontalAlignment', ...
    'center','VerticalAlignment', 'top')
    if save_flag
        saveas(gcf,[figures_folder filesep duration{j} '-' experiment_text{experiment} '_acc_sw',sufix],'png')
        saveas(gcf,[figures_folder filesep duration{j} '-' experiment_text{experiment} '_acc_sw',sufix],'fig')
    end
    
    %Accuracy per time interval
    figure
    for i=1:n_subjects(j)
        subplot(2,round(n_subjects(j)/2),i)
        hold on
        plot(analysis.(duration{j}).time_dist, analysis.(duration{j}).acc_bin(i,:),'rs-','LineWidth',2)
        xlabel('Time Interval(s)')
        ylabel('Accuracy')
        set(gca,'XTick',plot_min(j):0.05:plot_max(j))
        ylim([0 1])
        xlim([plot_min(j) plot_max(j)])
       
        
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,[task_version '-' duration{j} ' - ' experiment_text{experiment} '-Accuracy'],'HorizontalAlignment', ...
        'center','VerticalAlignment', 'top')
    if save_flag
        saveas(gcf,[figures_folder filesep duration{j} '-' experiment_text{experiment} '_acc_bin',sufix],'png')
        saveas(gcf,[figures_folder filesep duration{j} '-' experiment_text{experiment} '_acc_bin',sufix],'fig')
    end
    
    %Bias/Var
    figure
    colorss={'ro' ,'bo', 'go', 'ko', 'yo', 'co'};
    for i=1:n_subjects(j)
        hold on
        plot(analysis.(duration{j}).biasRMS(i), sqrt(analysis.(duration{j}).VAR(i)), colorss{i});
        ylim([0 1]);
        xlim([0 1]);
        
    end
    
    
    %Accuracy per time interval population
    figure
    hold on
    errorbar(analysis.(duration{j}).time_dist,analysis.(duration{j}).population.acc_bin_m,analysis.(duration{j}).population.acc_bin_se,'rs-','LineWidth',2)
    xlabel('Time Interval(s)')
    ylabel('Accuracy')
    set(gca,'XTick',plot_min(j):0.05:plot_max(j))
    ylim([0 1])
    xlim([plot_min(j) plot_max(j)])
    title([task_version '-' duration{j} ' - ' experiment_text{experiment} '-Population Accuracy'])
    if save_flag
        saveas(gcf,[figures_folder filesep duration{j} '-' experiment_text{experiment} '_acc_bin_pop',sufix],'png')
        saveas(gcf,[figures_folder filesep duration{j} '-' experiment_text{experiment} '_acc_bin_pop',sufix],'fig')
    end
    
    
    figure
    for i=1:n_subjects(j)
        subplot(2,round(n_subjects(j)/2),i)
        hold on
        errorbar(analysis.(duration{j}).time_dist, analysis.(duration{j}).est_m(i,:), analysis.(duration{j}).est_std(i,:),'ro')
        plot(analysis.(duration{j}).time_dist,polyval(analysis.(duration{j}).fit_params(i,:),analysis.(duration{j}).time_dist),'r-','LineWidth',1.5)
        plot([plot_min(j) plot_max(j)],[ plot_min(j) plot_max(j)],'--k','LineWidth',0.5)
        axis([plot_min(j) plot_max(j) plot_min(j) plot_max(j)])
        xlabel('Time Interval(s)')
        ylabel('Estimated Time(s)')
        set(gca,'XTick',plot_min(j):0.05:plot_max(j))
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,[task_version '-' duration{j} ' - ' experiment_text{experiment}],'HorizontalAlignment', ...
    'center','VerticalAlignment', 'top')

    if save_flag
        saveas(gcf,[figures_folder filesep duration{j} '-' experiment_text{experiment} '_EstVStrial',sufix],'png')
        saveas(gcf,[figures_folder filesep duration{j} '-' experiment_text{experiment} '_EstVStrial',sufix],'fig')
    end
end

% Population

for j=2
    if n_subjects(j)>1
    figure
    hold on
    errorbar(analysis.(duration{j}).time_dist,analysis.(duration{j}).population.est_m,analysis.(duration{j}).population.error_se,'o')
    plot([plot_min(j) plot_max(j)],[ plot_min(j) plot_max(j)],'--k','LineWidth',0.5)
    plot(analysis.(duration{j}).time_dist,polyval(analysis.(duration{j}).population.fit_params,analysis.(duration{j}).time_dist),'b-','LineWidth',1.5)
    axis([plot_min(j) plot_max(j) plot_min(j) plot_max(j)])
    xlabel('Time Interval(s)')
    ylabel('Estimated Time(s)')
    title([task_version '-' duration{j} '-' experiment_text{experiment} '- Population'])

    set(gca,'XTick',plot_min(j):0.05:plot_max(j))
    if save_flag
        saveas(gcf,[figures_folder filesep duration{j} '-' experiment_text{experiment} '_EstVStrial_pop', sufix],'png')
        saveas(gcf,[figures_folder filesep duration{j} '-' experiment_text{experiment} '_EstVStrial_pop', sufix],'fig')
    end
    end
end

% Plot Sliding window
for j=2
    
    % Individuals
    if n_subjects(j)>1
        figure
        for i=1:n_subjects(j)
            subplot(2,round(n_subjects(j)/2),i)
            hold on
            plot(1:n_trials-sliding_window,squeeze(analysis.(duration{j}).SW.fit_params(i,1,:)),'b-')
            plot(1:n_trials-sliding_window,squeeze(analysis.(duration{j}).SW.fit_params(i,2,:)),'r-')
            xlabel('Trial')
            ylabel('Slope / Bias')
            ylim([-1 1])
            title([task_version '-' duration{j} '-' experiment_text{experiment} ' Subject'])
            legend('Slope','Bias')
            
        end
        if save_flag
            saveas(gcf,[figures_folder filesep duration{j} '_' experiment_text{experiment} '_Slope_time_ind', sufix],'png')
            saveas(gcf,[figures_folder filesep duration{j} '_' experiment_text{experiment} '_Slope_time_ind', sufix],'fig')
        end
    end
    
    
    %Population means
    figure
    hold on
    plot(analysis.(duration{j}).SW.pop_fit_params(:,1),'b-')
    plot(analysis.(duration{j}).SW.pop_fit_params(:,2),'r-')
    xlabel('Trial')
    ylabel('Slope / Bias')
    title([task_version '-' duration{j} '-' experiment_text{experiment} '- Population'])
    %     set(gca,'XTick',plot_min(j):0.05:plot_max(j))
    if save_flag
        saveas(gcf,[figures_folder filesep duration{j} '_' experiment_text{experiment} '_Slope_time', sufix],'png')
        saveas(gcf,[figures_folder filesep duration{j} '_' experiment_text{experiment} '_Slope_time', sufix],'fig')
    end
    legend('Slope','Bias')
    

    
    %Sliding Window
   figure
    x = analysis.(duration{j}).time_dist;
    hplot = plot(x,polyval(analysis.(duration{j}).SW.pop_fit_params(1,:),x));
    hold on
    hplot2 = plot(x,analysis.(duration{j}).SW.pop_est_m(1,:),'o');
    plot([plot_min(j) plot_max(j)],[ plot_min(j) plot_max(j)],'--k','LineWidth',0.5)
    axis([plot_min(j) plot_max(j) plot_min(j) plot_max(j)])
    h = uicontrol('style','slider','units','pixel',...
                 'min',1,'max',n_trials-sliding_window,'val',1,...
                 'sliderstep',[1/(n_trials-sliding_window-1) 1], 'position',[50 10 300 20]);
    addlistener(h,'ActionEvent',@(hObject, event) makeplot(hObject,event,analysis,duration,x,j,hplot,hplot2));
          
 
end


if save_flag   
    save([save_folder filesep task_version '_'  num2str(experiment) '_analysis' sufix '.mat'],'analysis')
end
end

