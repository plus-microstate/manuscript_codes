function assessing_non_stationarity() 

global subject_IDs
global data_dir

%% Make surrogate data and get statistics
load('cluster_output','maps','kvec','kopt') ; 
for i = 1:length(subject_IDs)
       
    for scan = 1:2
        % load the source data
        cfg = struct ; 
        cfg.datafile = sprintf('%s/MEG-rest/sub%s-rest-%d.edf',data_dir,subject_IDs{i},scan) ; 
        source = ft_preprocessing(cfg) ; 

        % make a microstate object
        ms = microstate.individual(source.trial{1}','source',source.time{1}) ; % make microstate individual object

        % calculate stats
        ms.maps = maps{kvec == kopt} ; 
        ms = ms.cluster_alignmaps ; 
        ms = ms.stats_all ; 

        % Make some surrogates
        Nsurrogate = 100 ;  
        sur.Cs = [] ; sur.Covs = [] ; sur.Durs = [] ; sur.MDs = [] ; 
        sur.GEVs = [] ; sur.Hs = [] ; sur.Occs = [] ; 
        for surr = 1:Nsurrogate

            % Make a new microstate object but with surrogate data
            ms2 = ms ; 
            ms2.data = surrogate_iAAFT(ms.data,10,true) ; 

            % Calcualte stats
            ms2 = ms2.stats_all ; 

            % Add values to surrogate distribution
            sur.Cs(surr) = ms2.stats.complexity.complexity ; 
            sur.Covs(:,surr) = ms2.stats.coverage ; 
            sur.Durs(:,surr) = ms2.stats.duration; 
            sur.MDs(surr) = ms2.stats.mean_duration ;
            sur.GEVs(surr) = ms2.stats.gev.gev ; 
            sur.Hs(surr) = ms2.stats.hurst ; 
            sur.Occs(:,surr) = ms2.stats.occurrence ;
        end
        
        % Add original values to sur structure
        sur.C = ms.stats.complexity.complexity ; 
        sur.Cov = ms.stats.coverage ; 
        sur.Dur = ms.stats.duration ; 
        sur.MD = ms.stats.mean_duration ; 
        sur.GEV  = ms.stats.gev.gev ; 
        sur.H = ms.stats.hurst ; 
        sur.Occ = ms.stats.occurrence ; 
        
        % Z-score original data vs surrogate data
        vars = {'C','MD','GEV','H'} ; 
        for j = 1:length(vars)
                
            % Get original value and surrogate distribution
            varss = sprintf('%ss',vars{j}) ; 
            X = sur.(vars{j}) ; 
            Xs = sur.(varss) ; 
            
            % KS test to check for normal distribution
            [~,out.(vars{j}).p_sur_normal(2*i-mod(scan,2))] = kstest(zscore(Xs)) ; 
            
            % Z-score
            out.(vars{j}).Z(2*i-mod(scan,2)) = (X - mean(Xs)) / std(Xs) ; 

        end
            
        vars = {'Cov','Dur','Occ'} ; 
        for j = 1:length(vars)

            % Get original value and surrogate distribution
            varss = sprintf('%ss',vars{j}) ; 
            X = nanstd(sur.(vars{j})) ; 
            Xs = nanstd(sur.(varss)) ; 
            
            % KS test to check for normal distribution
            [~,out.(vars{j}).p_sur_normal(2*i-mod(scan,2))] = kstest(zscore(Xs)) ; 
            
            % Z-score
            out.(vars{j}).Z(2*i-mod(scan,2)) = (X - mean(Xs)) / std(Xs) ;

        end  
        
    end
   
end
       
%% Test for significant difference

% Analyse
vars = {'GEV','MD','C','H','Cov','Dur','Occ'} ;
clear X 
for i = 1:length(vars)
    % Get z-scores
    x = out.(vars{i}).Z ;
    % Average over the two scans
    x1 = reshape(x',2,length(x)/2); x1 = mean(x1) ; 
    % t-test for zero mean hypothesis
    [~,p_zscore(i)] = ttest(x1) ;
    
    % Save X as a column vector for plotting
    X(:,i) = x1 ; 
end
q_zscore = mafdr(p_zscore,'bhfdr','true') 


%% Plot Fig S3

figure('Name','Supp Fig S3','NumberTitle','off')

ms_path = microstate.functions.toolbox_path ;
addpath(fullfile(ms_path,'+external','brewermap'))
cmap = brewermap(64,'RdBu') ;
cmap = cmap(1:32,:) ; 


subplot(2,1,1)
plt = microstate.external.boxplot_LT(atanh(S0),...
        'Box_MedianFun',@(x)nanmedian(x),...
        'Box_WhiskerMinFun',@(x)max(min(x),nanmedian(x)-1.5*iqr(x)),...
        'Box_WhiskerMaxFun',@(x)min(max(x),nanmedian(x)+1.5*iqr(x)),...
        'Box_FaceColor',cmap(end,:),...
        'Box_MedianLineColor',cmap(6,:),'Box_MedianLineWidth',1,...
        'Points_Marker','.') ; 
set(gca,'YTick',atanh([0,0.5,0.75,0.9,0.99,0.999,0.9999]),'YTickLabel',[0,0.5,0.75,0.9,0.99,0.999,0.9999],'YGrid','on')
ylabel('Map similarity')
    
    
vars1 = {'GEV','MD','C','H','stdCov','stdDur','stdOcc'} ; 
plt = microstate.external.boxplot_LT(X,...
    'YLabel','z-score','XLabel',vars1,...
    'Box_FaceColor',cmap([12,12,24,12,32,24,12],:),... % Chosen these for our p-values, see fig caption in paper
    'Box_MedianFun',@(x)nanmedian(x),...
    'Box_WhiskerMinFun',@(x)max(min(x),nanmedian(x)-1.5*iqr(x)),...
    'Box_WhiskerMaxFun',@(x)min(max(x),nanmedian(x)+1.5*iqr(x)),...
    'Box_MedianLineColor',[0,0,0],'Box_MedianLineWidth',1,...
    'Points_Marker','.') ; 
set(gca,'YGrid','on')


end
