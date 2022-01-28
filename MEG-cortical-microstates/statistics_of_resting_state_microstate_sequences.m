function statistics_of_resting_state_microstate_sequences(subject_IDs,dirs)
% This section reproduces the results of section 3.3 and Fig 3. Note the
% analysis of non-stationarity is in a separate function
% (assessing_non_stationarity.m). The info from supplementary table 3 is
% available in cohort.stats

% Load in the results of the clustering
load('cluster_output','cohort')
    
%% Plot syntax matrix

figure('Name','Figure 3A','NumberTitle','off')
    
% Get median syntax matrix
all_syntax = [cohort.stats(1).syntax.matrix ; cohort.stats(2).syntax.matrix] ; 
P = squeeze(nanmedian(all_syntax)) ; 

% Plot
ms = microstate.individual ; 
ms.stats.syntax.matrix = P ; 
ms.plot('syntax_matrix')

    
%% Plot coverage, duration, occurrence

figure('Name','Figure 3B-D','NumberTitle','off')

% Get colour map
ms_path = microstate.functions.toolbox_path ; 
addpath(fullfile(ms_path,'+external','brewermap'))
cmap = brewermap(64,'RdBu') ;
cmap = cmap(1:32,:) ;
colormap(gca,flipud(cmap)) ;
    
% Loop over the three variables and plot
vars = {'coverage','duration','occurrence'} ; 
ylbls = {'Coverage (%)','Duration [ms]','Occurrence [s^-^1]'} ; 
for j = 1:3
        
    subplot(3,1,j)
    clear X
    for i = 1:2 % Loop over scan 1 and scan 2
        X(:,:,i) = cohort.stats(i).(vars{j}) ;
    end
    X = mean(X,3) ; % Average across scans so one value per participant
    if j==2
        X=1000*X ; % s->ms
    elseif j==1
        X = 100*X ; % decimal -> percent
    end
    plt = microstate.external.boxplot_LT(X,...
        'plotPoints',false,...
        'Box_FaceColor',cmap(end,:),...
        'Box_MedianLineColor',cmap(6,:),'Box_MedianLineWidth',1) ; 
    ylabel(ylbls{j})
    set(gca,'YGrid','on')

end        
    
end