function anal_mult = multiple_analysis_IE(sufix,save_flag)
%% Loading data and grouping it
%Path to a folder with all data split in folders called "auditory","visual"...
%path_folder = 'C:\Users\Rodrigo\Documents\INDP2015\Project\DATA';

analysis_folder = 'C:\Users\Rodrigo\Documents\INDP2015\Project\ANALYSIS';
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
    if strcmp(file.analysis.info.modality,'visual')
        mdlt=2;       
    else
        mdlt=1;
    end
        
    anal_mult{mdlt,file.analysis.info.experiment}=file.analysis;  
end
duration={'short','long'};

% Scores organization
for j=2
    for exp=1:2
        aud_scores{exp,j}=anal_mult{1,exp}.(duration{j}).score;
        vis_scores{exp,j}=anal_mult{2,exp}.(duration{j}).score;
    end
end

for j=2
    for exp=1:2
        names_score{exp,j} = fieldnames(aud_scores{exp,j})
    end
end


end