function microstate_response_to_auditory_stimulus() 

global subject_IDs
global data_dir


%% Estimate maps
% In the resting-state data, there was 1 x 8min scan from each recording
% session, and we used 5000 peaks from scan 1 for clustering. In the task
% data, there are two consecutive 5 min scans from each recording session,
% so here we will cluster from both scans from session 1, i.e. use 2500 GFP
% peaks from scan 1 and 2500 GFP peaks from scan 2. 

c = microstate.cohort ; % make microstate cohort

for i = 1:length(subject_IDs)
    
    % Both scans from session 1, i.e. scans 1-2
    for scan = 1:2
        % load the source data
        cfg = struct ; 
        cfg.datafile = sprintf('%s/MEG-auditory/sub%s-auditory-%d.edf',data_dir,subject_IDs{i},scan) ; 
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

        % add the microstate object to the cohort
        c = c.add_individuals(ms,scan,2500) ; 
    end
    clear source

end

[c,kopt,kvec,maps,gev] = c.cluster_globalkoptimum('kvec',2:40) ; 
    
% Plot Fig 5A
figure('Name','Figure 5A','NumberTitle','off') ; 
plot(kvec,gev,'k') ; 
hold on
plot(kopt,gev(kvec==kopt),'+r')
set(gca,'XTick',kvec) ; 
xlabel('k')
ylabel('GEV')
grid on
    
    
%% Backfit to full time courses across all scans

% Align maps to the resting-state maps
maps_task = maps ; kopt_task = kopt ; 
load('cluster_output','kopt','maps')
rskopt = kopt ; rsmaps = maps{kvec==rskopt} ; 
maps = maps_task ; clear maps_task ; kopt = kopt_task ; clear kopt_task ; 
[omaps,mapsim,order1] = microstate.functions.align_maps2template(maps{kvec==kopt},rsmaps,'source') ;
        
% Back-fit to the full time course
c = microstate.cohort ; 
for i = 1:length(subject_IDs) 
    for scan = 1:4
        
        % load the source data
        cfg = struct ; 
        cfg.datafile = sprintf('%s/MEG-auditory/sub%s-auditory-%d.edf',data_dir,subject_IDs{i},scan) ; 
        source = ft_preprocessing(cfg) ;
        
        % make a microstate object
        ms = microstate.individual(source.trial{sesh}','source',source.time{sesh}) ; % make microstate individual object 

        % backfit maps to data
        ms.maps = omaps ; % maps{kvec == kopt} ; 
        ms = ms.cluster_alignmaps ; 
        
        % Calculate statistics (including GEV)
        ms = ms.stats_all ; 

        % add the microstate object to the cohort
        c = c.add_individuals(ms,sprintf('scan%d_',scan)) ; 
        clear source
        
    end
end

% Collate stats for cohort
c = c.cohort_stats ; 
cohort = c ; 
save('cluster_output_task','cohort','kopt','kvec','maps','gev')

%% Plot Fig 5B & Compare some stats

% Plot Fig 5B
figure('Name','Figure 5B','NumberTitle','off') ; 


% Average each of the stats across the 4 scans, giving one value per
% participant
stats_task = microstate.functions.average_stats_across_conditions(cohort.stats) ;

% Load in the resting-state data
rs = load('cluster_output.mat') ; 
stats_rest = microstate.functions.average_stats_across_conditions(rs.cohort.stats) ;

% compare some global stats
p_stats = struct ; 
fn = {'mean_duration','hurst'} ; 
for i = 1:length(fn)
    p_stats.(fn{i}) = signrank(stats_task.(fn{i}),stats_rest.(fn{i})) ;
end
fn = {'complexity','gev'} ; 
for i = 1:length(fn)
    p_stats.(fn{i}) = signrank(stats_task.(fn{i}).(fn{i}),stats_rest.(fn{i}).(fn{i})) ;
end
p_stats % print the output to the console

%% Plot Fig 5C
    
figure('Name','Figure 5C','NumberTitle','off') ; 
lay = load('layout.mat') ; 

cohort.globalmaps = cohort.individual(1).maps ; 
cohort.plot('globalmaps',lay,'cscale',[0.5,1]) ; 


%% Compare to timing of stimuli
% Note that this is the way we did this in the manuscript, which is complex
% and long winded. +microstate has since been updated to include specific
% functionality to run this type of analysis, which is described in
% Tutorial 5 of the toolbox. It is recommended you follow tutorial 5 if
% trying to replicate this analysis in your own data, but the code from the
% manuscript is shown below. 

tpre = 0.100 ; % 100 ms 
Nlag = 90 ; % approx 350 ms at 256 Hz
statecount = zeros(1,size(omaps,2)) ;  
stimcount = zeros(Nlag,size(omaps,2)) ; 

msg = [] ; 

% Loop over each time lag
for lag = 1:Nlag
    % Can be time consuming, so print output
    repmat('\b',1,length(msg)) ; 
    msg = sprintf('lag %d of %d',lag,Nlag) ; 
    fprintf(msg) ; 
    
    % indcount just keeps track of the current scan. We have 120 scans 
    % total (30 participants, 2 sessions, 2 scans per session), so indcount
    % just counts from 1-120. 
    indcount = 0 ; 
    
    % Loop over participants
    for i = 1:length(subject_IDs)

        % Loop over scans
        for scan = 1:4
            
            % load the stimulus timings - UPDATE WHEN WRITTEN
            stim = load(sprintf('%s/%s/aud_scan%d/stimulus.mat',subject_IDs{i},scan)) ;
            indcount = indcount+1 ; 
            
            % Get post-stimulus time
            tstim = stim.time ;
            tstim = tstim+(1/256)*lag ; % note 256 used here as data sampled at 256 Hz

            lbl = cohort_task.individual(indcount).label ;
            if lag == 1 ; 
                idxpre = false(1,length(lbl)) ; 
                for stimi = 1:length(tstim) ; 
                    idxpre = idxpre | ...
                        (...
                        (cohort_task.individual(indcount).time >= (tstim(stimi)-tpre))...
                        & ...
                        (cohort_task.individual(indcount).time < tstim(stimi))...
                        ) ; 
                end

                statecount = statecount+histcounts(lbl(idxpre),(0:(size(omaps,2)))+0.5) ; 
            end

            lbl = interp1(cohort_task.individual(indcount).time,lbl,tstim,'next') ; 
            stimcount(lag,:) = stimcount(lag,:)+histcounts(lbl,(0:(size(omaps,2)))+0.5) ; 

        end
    end

end

    
% chi2 test
Nstim = sum(stimcount,2) ; 
O = stimcount ; 
E = Nstim.*(statecount/sum(statecount)) ; 
chi2 = nansum(((O-E).^2)./E,2) ;
p_stats.stim.chi2 = chi2 ; 
p = 1-chi2cdf(chi2,8) ; 
p_stats.stim.p_chi2 = p ; 
    
% find Bonferroni corrected threshold
al2 = chi2inv(1-0.05/90,8) ; 
    
% find peak
idx = find((chi2>al2) & (chi2>max(chi2)/2)) ; 
    
    
% Pearson residuals
Ox = O(idx,:) ; Ex = E(idx,:) ;  
Res = (Ox-Ex)./sqrt(Ex) ; 


%% Plot Figure 6
figure('Name','Figure 6A','NumberTitle','off') ; 
plot((1/256)*(1:Nlag),chi2,'k') ; 
xlabel('Time [s]')
ylabel('\chi^2')
xlim([0,0.35])
grid on
hold on
plot([0,0.35],[al2,al2],'r--') ; 
fill((1/256)*[idx(1) , idx(1) , idx(end), idx(end)] , [yl(1),yl(2),yl(2),yl(1)],'k','FaceAlpha',0.1,'LineStyle','none')



figure('Name','Figure 6B','NumberTitle','off')
b = bar(mean(Res)) ; 
b(1).FaceColor = [0.0196    0.4431    0.6902] ; 
grid on
set(gca,'XTick',1:9,'XTickLabel',{'Front','FrontTemp','Vis','Orb','Par','L-Temp','R-Temp','L-SensMot','R-SensMot'})
ylabel('Pearson Residual')
xtickangle(40)
hold on
xl = xlim ; 
plot(xl,[2.77,2.77],'r--')
plot(xl,-[2.77,2.77],'r--')
    


        
end
