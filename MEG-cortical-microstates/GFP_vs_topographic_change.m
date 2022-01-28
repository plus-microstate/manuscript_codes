function GFP_vs_topographic_change() ; 

global data_dir
global subject_IDs

GMD = [] ; GFP = [] ; GFPs = [] ; 
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
    ms.data = ms.data/std(ms.data(:)) ; % global normalization
    ms = ms.add_bad_samples(bad_samples) ; 
    clear source
    
    % GFP
    ms = ms.calculate_gfp ;
    
    % GMD
    gmd = zeros(1,length(ms.time)-1) ;
    fun = microstate.functions.map_similarity_funhandle(ms.modality) ; 
    for t = 1:length(ms.time)-1
        R = fun(ms.data(t,:),ms.data(t+1,:)) ; 
        gmd(t) = 1-R.^2 ; 
    end
    
    gfp = ms.gfp(1:end-1) ; 
    gmd(isnan(gfp)) = [] ; 
    gfp(isnan(gfp)) = [] ;
    
    % Plot
    if i == 1
        figure('Name','Supp Figure S4','NumberTitle','off')
        
        subplot(2,2,1)
        plot(gfp,gmd,'.k','MarkerSize',0.5)
        xlabel('GFP')
        ylabel('Topographic change')
        title('Observed data')
        xlim([0 80])
        ylim([0 0.8])
        
        subplot(2,2,2)
        plot(log(gfp),log(gmd),'.k','MarkerSize',0.5)
        xlabel('log(GFP)')
        ylabel('log(Topographic change)')
        title('Observed data (log)')
        xlim([1 5])
        ylim([-6 0])
    
        gfps = gfp(randperm(length(gfp))) ; 
        
        subplot(2,2,3)
        plot(gfps,gmd,'.k','MarkerSize',0.5)
        xlabel('GFP')
        ylabel('Topographic change')
        title('Shuffled data')
        xlim([0 80])
        ylim([0 0.8])
        
        subplot(2,2,4)
        plot(log(gfps),log(gmd),'.k','MarkerSize',0.5)
        xlabel('log(GFP)')
        ylabel('log(Topographic change)')
        title('Shuffled data (log)')
        xlim([1 5])
        ylim([-6 0])
    end
        
    % Concatenate 
    r(i) = corr(log(gfp)',log(gmd)') ; 
    
        
end

mean_loglog_correlation = mean(r)
sem_loglog_correlation = std(r)/sqrt(length(r))