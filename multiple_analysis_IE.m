function anal_mult = multiple_analysis_IE(sufix,save_flag)
%% Loading data and grouping it
%Path to a folder with all data split in folders called "auditory","visual"...
%path_folder = 'C:\Users\Rodrigo\Documents\INDP2015\Project\DATA';

analysis_folder = 'C:\Users\Rodrigo\Documents\INDP2015\Project\ANALYSIS';
figures_folder = 'C:\Users\Rodrigo\Documents\INDP2015\Project\MULT_ANALYSIS\figures';

%  analysis_folder = '/Users/baylorbrangers/Desktop/Subject_Data_Bayes/Data';

path_data = [analysis_folder filesep , '*.mat'];

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
    file = load([analysis_folder filesep names{i}]);
    if file.analysis.info.experiment==3
        if strcmp(file.analysis.info.sufix,'a')
            mdlt=2;
        elseif strcmp(file.analysis.info.sufix,'b')
            mdlt=1;
        end
    else
        if strcmp(file.analysis.info.modality,'visual')
            mdlt=2;
        else
            mdlt=1;
        end
    end
    
    anal_mult{mdlt,file.analysis.info.experiment}=file.analysis;
end
duration={'short','long'};
modality={'Auditory','Visual','VisualB'};


%% estimate
for j=2
    plot_max(j)=max(anal_mult{1,2}.(duration{j}).time_dist)+0.05;
    plot_min(j)=min(anal_mult{1,2}.(duration{j}).time_dist)-0.05;
end

figure
for j=2
    for mdlt=1:2
        subplot(1,3,mdlt)
        hold on
        errorbar(anal_mult{mdlt,1}.(duration{j}).time_dist,anal_mult{mdlt,1}.(duration{j}).population.est_m,anal_mult{mdlt,1}.(duration{j}).population.error_se,'ro')
        errorbar(anal_mult{mdlt,2}.(duration{j}).time_dist,anal_mult{mdlt,2}.(duration{j}).population.est_m,anal_mult{mdlt,2}.(duration{j}).population.error_se,'bo')
        
        plot(anal_mult{mdlt,1}.(duration{j}).time_dist,polyval(anal_mult{mdlt,1}.(duration{j}).population.fit_params,anal_mult{mdlt,1}.(duration{j}).time_dist),'r-','LineWidth',1.5)
        
        plot(anal_mult{mdlt,2}.(duration{j}).time_dist,polyval(anal_mult{mdlt,2}.(duration{j}).population.fit_params,anal_mult{mdlt,2}.(duration{j}).time_dist),'b-','LineWidth',1.5)
        
        plot([plot_min(j) plot_max(j)],[ plot_min(j) plot_max(j)],'--k','LineWidth',0.5)
        axis([plot_min(j) plot_max(j) plot_min(j) plot_max(j)])
        xlabel('Time Interval(s)')
        ylabel('Estimated Time(s)')
        title([modality{mdlt} '- Population'])
        
        set(gca,'XTick',plot_min(j):0.05:plot_max(j))
    end
    legend({'Aud2Vis','Vis2Aud'})
    
    subplot(1,3,3)
    hold on
    errorbar(anal_mult{1,3}.(duration{j}).time_dist,anal_mult{1,3}.(duration{j}).population.est_m,anal_mult{1,3}.(duration{j}).population.error_se,'go')
    errorbar(anal_mult{mdlt,3}.(duration{j}).time_dist,anal_mult{mdlt,3}.(duration{j}).population.est_m,anal_mult{mdlt,3}.(duration{j}).population.error_se,'co')
    plot(anal_mult{1,3}.(duration{j}).time_dist,polyval(anal_mult{1,3}.(duration{j}).population.fit_params,anal_mult{1,3}.(duration{j}).time_dist),'g-','LineWidth',1.5)
    plot(anal_mult{mdlt,3}.(duration{j}).time_dist,polyval(anal_mult{mdlt,3}.(duration{j}).population.fit_params,anal_mult{mdlt,3}.(duration{j}).time_dist),'c-','LineWidth',1.5)
    plot([plot_min(j) plot_max(j)],[ plot_min(j) plot_max(j)],'--k','LineWidth',0.5)
    axis([plot_min(j) plot_max(j) plot_min(j) plot_max(j)])
    xlabel('Time Interval(s)')
    ylabel('Estimated Time(s)')
    title([modality{mdlt} '- Population'])
    
    set(gca,'XTick',plot_min(j):0.05:plot_max(j))
    legend({'Vis2VisII','Vis2VisI'})
    if save_flag
        saveas(gcf,[figures_folder filesep 'mult_EstVStrial_pop', sufix],'png')
        saveas(gcf,[figures_folder filesep 'mult_EstVStrial_pop', sufix],'fig')
    end
end



%% Accuracy

figure
for j=2
    for mdlt=1:2
        subplot(1,3,mdlt)
        hold on
        errorbar(anal_mult{mdlt,1}.(duration{j}).time_dist , anal_mult{mdlt,1}.(duration{j}).population.acc_bin_m , anal_mult{mdlt,1}.(duration{j}).population.acc_bin_se,'rs-','LineWidth',2)
        
        errorbar(anal_mult{mdlt,2}.(duration{j}).time_dist , anal_mult{mdlt,2}.(duration{j}).population.acc_bin_m , anal_mult{mdlt,2}.(duration{j}).population.acc_bin_se,'bs-','LineWidth',2)

        xlabel('Time Interval(s)')
        ylabel('Accuracy')
        set(gca,'XTick',plot_min(j):0.05:plot_max(j))
        ylim([0 1])
        xlim([plot_min(j) plot_max(j)])
        title([modality{mdlt} ' - Population Accuracy'])
    end
    legend({'Aud2Vis','Vis2Aud'})
    
    subplot(1,3,3)
    hold on
    errorbar(anal_mult{1,3}.(duration{j}).time_dist , anal_mult{1,3}.(duration{j}).population.acc_bin_m , anal_mult{1,3}.(duration{j}).population.acc_bin_se,'gs-','LineWidth',2)
    
    errorbar(anal_mult{mdlt,3}.(duration{j}).time_dist , anal_mult{mdlt,3}.(duration{j}).population.acc_bin_m , anal_mult{mdlt,3}.(duration{j}).population.acc_bin_se,'cs-','LineWidth',2)
    
    xlabel('Time Interval(s)')
    ylabel('Accuracy')
    set(gca,'XTick',plot_min(j):0.05:plot_max(j))
    ylim([0 1])
    xlim([plot_min(j) plot_max(j)])
    title([modality{mdlt} ' - Population Accuracy'])
    legend({'Vis2VisII','Vis2VisI'})

end

%differene between the two
for j=2
    for mdlt=1:2
        acc_bin_m_dif(mdlt,:) = anal_mult{mdlt,1}.(duration{j}).population.acc_bin_m -  anal_mult{mdlt,2}.(duration{j}).population.acc_bin_m;
    end
end
zvalues=zscore(acc_bin_m_dif,1,2);

%% Scores organization
for j=2
    for exp=1:2
        aud_scores{exp,j}=anal_mult{1,exp}.(duration{j}).score;
        vis_scores{exp,j}=anal_mult{2,exp}.(duration{j}).score;
    end
    vis_scores{3,j}=anal_mult{2,3}.(duration{j}).score;
    vis_scores{4,j}=anal_mult{1,3}.(duration{j}).score;
end

for j=2
    for exp=1:2
        names_score_aud{exp,j} = fieldnames(aud_scores{exp,j});
        names_score_vis{exp,j} = fieldnames(vis_scores{exp,j});
    end
    names_score_vis{3,j} = fieldnames(vis_scores{3,j});
    for i= 1:length(names_score_vis{3,j})
    names_score_vis{4,j}{i,1} = [names_score_vis{3,j}{i} 'II'];
    end
end
n_subjects=zeros(2,2);
for j=2
    for exp=1:2
        n_subjects(exp,j) = length(names_score_aud{exp,j});
    end
end
aud_scores_plot=[];
vis_scores_plot=[];

for j=2
    k=1;
    for exp=1:2
        for i=1:n_subjects(exp,j)
            aud_scores_plot(j,k) = aud_scores{exp,j}.(names_score_aud{exp,j}{i,1});
            vis_scores_plot(j,k) = vis_scores{exp,j}.(names_score_vis{exp,j}{i,1});
            names_scores_plot{j,k} = names_score_aud{exp,j}{i,1};
            k=k+1;
        end
    end
    
    for i=1:length(names_score_vis{3,j})
        
        visvis_scores_plot(1,i) = vis_scores{3,j}.(names_score_vis{3,j}{i,1});   
        visvis_scores_plot(2,i) = vis_scores{4,j}.(names_score_vis{3,j}{i,1}); 
    end
      for i=1:length(names_score_vis{3,j})
    vis_names_scores_plot{j,i} = names_score_vis{3,j}{i,1};
    vis_names_scores_plot{j,i+length(names_score_vis{3,j})} = names_score_vis{4,j}{i,1};
      end
end



[max_scores subj_ind]= max(max(vis_scores_plot + aud_scores_plot));
['The winner is' (names_scores_plot(2,subj_ind)) ' with ' max_scores 'points']
figure
for j=2
    %Auditory
    subplot(3,1,1)
    bar(aud_scores_plot(j,:))
    title([modality{1}  ' - '  duration{j} '- Scores'])
    ylim([0 200])
    set(gca,'XTickLabel',names_scores_plot(j,:) )
    %Visual
    subplot(3,1,2)
    bar(vis_scores_plot(j,:))
    title([modality{2}  ' - '  duration{j} '- Scores'])
    ylim([0 200])
    set(gca,'XTickLabel',names_scores_plot(j,:) ) 
    subplot(3,1,3)
    hold on
    bar([1:length(names_score_vis{3,j})],visvis_scores_plot(1,:))
    bar([length(names_score_vis{3,j})+1:(length(names_score_vis{3,j})+length(names_score_vis{3,j}))],visvis_scores_plot(2,:))
    title([modality{2}  ' - '  duration{j} '- Scores'])
    ylim([0 200])
    set(gca,'XTick',1:(2*length(names_score_vis{3,j})) )  
    set(gca,'XTickLabel',vis_names_scores_plot(j,:) ) 
end
if save_flag
    saveas(gcf,[figures_folder filesep 'Scores',sufix],'png')
    saveas(gcf,[figures_folder filesep 'Scores',sufix],'fig')
end

%% 

% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%     text(0.5, 1,[modality{2}  ' - '  duration{j} '- Scores'],'HorizontalAlignment', ...
%         'center','VerticalAlignment', 'top')

end