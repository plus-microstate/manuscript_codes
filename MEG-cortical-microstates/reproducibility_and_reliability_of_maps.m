function reproducibility_and_reliability_of_maps() 

global subject_IDs
global data_dir

% Load in cluster output
load('cluster_output','maps','kvec','kopt')
maps = maps{kvec == kopt} ; 
    
%% Repeat the clustering analysis multiple times with different k-means
% seeds and different GFP samples
nReps = 250 ; 
for rep = 1:nReps
    
    c = microstate.cohort ; % make microstate cohort
    for i = 1:length(subject_IDs)
        % load the source data
        cfg = struct ; 
        cfg.datafile = sprintf('%s/MEG-rest/sub%s-rest-1.edf',data_dir,subject_IDs{i}) ; 
        source = ft_preprocessing(cfg) ; 


        % get the bad samples
        artfctdef = jsondecode(fileread(sprintf('%s/MEG-rest/Artifacts/artfct-sub%s-rest-1.json',data_dir,subject_IDs{i}))) ; 
        bad_samples = [] ; 
        for mth = {'clip','jump','zscore'}
            for j = 1:size(artfctdef.(mth{1}).artifact,1)
                bad_samples = [bad_samples , (artfctdef.(mth{1}).artifact(j,1)-5):(artfctdef.(mth{1}).artifact(j,2)+5)] ; 
            end
        end
        bad_samples = unique(bad_samples) ; 

        % make a microstate object
        ms = microstate.individual(source.trial{1}','source',source.time{1}) ; % make microstate individual object
        ms = ms.add_bad_samples(bad_samples) ; 
        clear source

        % add 5000 GFP peaks from the microstate object to the cohort
        c = c.add_individuals(ms,[],5000) ; 

    end
    
    % Cluster repetition
    c = c.cluster_global(kopt,'kmeans_replicates',1) ; 
    repmaps(:,:,rep) = microstate.functions.align_maps2template(c.globalmaps,maps,'source') ; 
end
        
%% Analyse with map similarity and silhouette coefficient

% Map similarity
simf = microstate.functions.map_similarity_funhandle('source') ;

% Loop over clusters 
for i = 1:size(maps,2)
    rmapsi = squeeze(repmaps(:,i,:)) ; 

    % Cluster correlation index
    S0(:,i) = simf(rmapsi',maps(:,i)')
end

mean_map_similarity = mean(S0(:))  
sem_map_similarity = std(S0(:))/sqrt(numel(S0)) 

% Calculate silhouette coefficients
for Ci = 1:size(repmaps,2)
    rmapsi = squeeze(repmaps(:,Ci,:)) ; 
        
    d = 1-simf(rmapsi',rmapsi') ; d=d-diag(diag(d)) ;  
    a(Ci,:) = (1/(size(d,1)-1))*sum(d) ; 
        
    B = inf(size(maps,2),nReps) ; 
    for Ck = 1:size(repmaps,2)
        if Ck==Ci , continue , end
        rmapsk = squeeze(repmaps(:,Ck,:)) ; 
        d = 1-simf(rmapsi',rmapsk') ; 
        B(Ck,:) = (1/(size(d,1)-1))*sum(d) ; 
    end
    b(Ci,:) = min(B) ; 
end
s = abs(b-a)./max(a,b) ; 

mean_silhouette = mean(s(:))  
sem_silhouette = std(s(:))/sqrt(numel(s)) 
        
    
%% Plot figure S2

figure('Name','Supp Figure S2','NumberTitle','off')
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
    
    
    
subplot(2,1,2)
plt = microstate.external.boxplot_LT(s',...
        'Box_MedianFun',@(x)nanmedian(x),...
        'Box_WhiskerMinFun',@(x)max(min(x),nanmedian(x)-1.5*iqr(x)),...
        'Box_WhiskerMaxFun',@(x)min(max(x),nanmedian(x)+1.5*iqr(x)),...
        'Box_FaceColor',cmap(end,:),...
        'Box_MedianLineColor',cmap(6,:),'Box_MedianLineWidth',1,...
        'Points_Marker','.')  ; 
set(gca,'YGrid','on')
ylabel('Silhouette Coefficient')
    

end
