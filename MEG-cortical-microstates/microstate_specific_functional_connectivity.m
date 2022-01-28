function microstate_specific_functional_connectivity() ;

global subject_IDs
global data_dir

% Load in the cluster maps
load('cluster_output','maps','kvec','kopt')

% Specify frequency bands
freqs = [4,8;8,13;13,30] ; 

%% Make the networks

% Loop over frequency bands
for frq = 1:size(freqs,1)
    
    % Initialize networks for this frequency
    nets = cell(2*length(subject_IDs),kopt) ; 
    
    % Loop over participants
    for i = 1:length(subject_IDs)
        
        % Loop over scans
        for scan = 1:2

            fprintf('Frq %d of %d, Participant %d of %d, Scan %d of %d\n',frq,size(freqs,1),i,length(subject_IDs),scan,2)

            % load the source data
            cfg = struct ; 
            cfg.datafile = sprintf('%s/MEG-rest/sub%s-rest-%d.edf',data_dir,subject_IDs{i},scan) ; 
            source = ft_preprocessing(cfg) ;
            
            % get the bad samples
            artfctdef = jsondecode(fileread(sprintf('%s/MEG-rest/Artifacts/artfct-sub%s-rest-%d.json',data_dir,subject_IDs{i},scan))) ; 
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

            % backfit maps to data
            ms.maps = maps{kvec == kopt} ; 
            ms = ms.cluster_alignmaps ; 

            % make networks: Note in the latest version of +microstate nets
            % have been moved to 3rd input, but in older versions nets was
            % first input
            % nets(2*i-mod(scan,2),:) = ms.networks_wpli(freqs(frq,:),true) ; 
            [~,~,nets(2*i-mod(scan,2),:)] = ms.networks_wpli(freqs(frq,:),true) ;      
            save(sprintf('nets_%d_%d_Hz',freqs(frq,1),freqs(frq,2)),'nets','-v7.3') ; 
            
            % -- Normalize networks against static for plotting --
            % make static network by assinging all points to the same
            % state
            ms.label = ones(size(ms.data,1),1) ; 

            % make networks
            [~,~,staticnet] = ms.networks_wpli(freqs(frq,:),false,[],true) ;  
            net = nets(2*i-mod(scan,2),:) ; 
            for j = 1:length(net)
                if ~isempty(net{j})
                    net{j} = (net{j}-repmat(staticnet,1,1,size(net{j},3)))./(1-net{j}) ;
                end
            end
            normnets(2*i-mod(scan,2),:) = net ;
            save(sprintf('nets_corrected_%d_%d_Hz',freqs(frq,1),freqs(frq,2)),'normnets','-v7.3') ; 
            
        end
    end
end


%% MVPA analysis
% Requires MVPA-Lite toolbox

for frq = 1:size(freqs,1)

    % Load in the networks
    load(sprintf('nets_%d_%d_Hz',freqs(frq,1),freqs(frq,2)),'nets')

    % MVPA ANALYSIS ------         
    % Get features: X is weighted degree of each node, G is the microstate.
    X = [] ; G = [] ; 
    for j = 1:size(nets,2)
        for i = 1:size(nets,1)
            n = nets{i,j} ; 
            for k = 1:size(n,3)
                nk = n(:,:,k) ; 
                if ~isempty(nk)
                    X = [X , nansum(nk,2)] ; 
                    G = [G,j] ; 
                end
            end
        end
    end

    % Make MVPA classifier
    cfg = struct ; 
    cfg.classifier = 'multiclass_lda' ; 
    cfg.sample_dimension = 2 ; 
    cfg.feature_dimension = 1 ;
    cfg.metric = 'accuracy' ; 
    [perf, result, testlabel] = mv_classify(cfg, X, G') ;

    % Get confusion matrix
    cfg.metric = 'confusion' ; 
    conf = mv_classify(cfg,X,G') ; 
        
    % Statistics
    cfg = struct ; 
    cfg.test = 'permutation' ;  
    stat = mv_statistics(cfg,result,X,G') ; 

    % Put it all into a structure
    mvpa.accuracy = perf ; 
    mvpa.confusion = conf ; 
    mvpa.result = result ; 
    mvpa.testlabel = testlabel ; 
    mvpa.stat = stat ; 
    save(sprintf('mvpa_%d_%d_Hz',freqs(frq,1),freqs(frq,2)),'mvpa')
    

    % --- PLOT FIG 4B ----
    if frq == 2 % Only plot for alpha
        % Get degree distributions for plotting
        for j = 1:max(G) ;
            ind = find((G==j)) ; 
            D(:,j) = nanmean(X(:,ind),2) ; 
        end

        % Get network edges
        load(sprintf('nets_corrected_%d_%d_Hz',freqs(frq,1),freqs(frq,2)),'normnets') ; 

        NET = cell(kopt,1) ; 
        for j = 1:size(normnets,2)

            % Get average normalized network
            normnet = zeros(size(X,1)) ; count = 0 ; 
            for i = 1:size(normnets,1)
                n = normnets{i,j} ; 
                for k = 1:size(n,3)
                    nk = n(:,:,k) ;  
                    if ~isempty(nk)
                        for m = 1:230 ; nk(m,m) = 0 ; end ;
                        normnet = normnet+nk ; count = count+1 ;
                    end
                end
            end
            normnet = normnet/count ; 

            % Threshold at 10%
            vals = sort(unique(normnet(:)),'descend') ; vals(isnan(vals)) = [] ; 
            pct = floor(0.01*length(vals)) ;
            try ; thrsh = vals(pct) ; catch ; thrsh = 0 ; end
            normnet = normnet.*(normnet > thrsh) ; 

            NET{j} = normnet ; 
        end
    
        figure('Name','Figure 4B','NumberTitle','off')
        layout = load('layout.mat') ; 
        microstate.functions.networks_plot(NET,layout) ; 
        
        fig = gcf ; 
        for j = 1:length(fig.Children)
            ax = fig.Children(j) ;
            gr = ax.Children(2) ; % Get the graph
            
            gr.NodeCData = D(:,j) ;
            gr.MarkerSize = 1.8 ; 
            gr.EdgeColor = [0,0,0] ; 
        end
        
        % Get colour map
        ms_path = microstate.functions.toolbox_path ;
        addpath(fullfile(ms_path,'+external','brewermap'))
        cmap = brewermap(64,'RdBu') ;
        cmap = cmap(1:32,:) ; 
        colormap(flipud(cmap)) ; 
        
    end
    % ---- END PLOT FIG 4B ---
    
end

%% Plot Fig 4A

figure('Name','Figure 4A','NumberTitle','off')
for frq = 1:size(freqs,1)
    subplot(1,3,frq)
    
    load(sprintf('mvpa_%d_%d_Hz',freqs(frq,1),freqs(frq,2)),'mvpa')
    
    imagesc(mvpa.confusion)
    colormap(flipud(cmap)) ; 
    xlabel('Predicted class')
    ylabel('True class')
    axis image
end

end