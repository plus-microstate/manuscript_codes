function comparing_state_estimation_methods_phase()
% Note: This code is only for a small part of section S2.3 Results from 
% the 4 activation models were given in the simulations() function. Only
% the results from the phase locking model are included here. 

% Simulate some data
rng default
N = 5 ; 
f = 1:1e-3:100 ;

% State 1: Incoherent alpha
A{1} = 1./f + exp(-(f-10).^2) ; 
Phi{1} = 2*pi*rand(N,length(f)) ; 

% State 2: Coherent alpha
A{2} = 1./f + exp(-(f-10).^2) ; 
Phi{2} = (f<8 | f>12).*A{1} + ... % incoherent background
    (f>=8 & f<=12).*(0.1*randn(5,length(f))+pi) ; % coherent alpha

% State 3: Incoherent beta
A{3} = 1./f + exp(-(f-17).^2) ; 
Phi{3} = Phi{1} ; 
      
  
%% Simulate some trials

% Some key info
K = length(A) ; % number of states 
Fs = 200 ; % Sampling freq
NTrl = 20 ; % number of trials
T = repmat(5*60*Fs,NTrl,1) ; % T for HMM-MAR toolbox

% States can last 20ms-1s
Tmin = round(0.02*Fs) ; 
Tmax = round(Fs) ; 

% Markov transitioning matrix
P = rand(K) ; P = P.*(1-eye(K)) ; P = P./sum(P,2) ; 
Pi = sum(P)/sum(sum(P)) ; % Prior probability of states

% Initialize
X = zeros(sum(T),N) ; % time series
G = nan(sum(T),1) ; % true states

for trl = 1:NTrl
    
    fprintf('Trial %d\n',trl) ; 
    
    % Generate initial state
    state = find(rand < cumsum(Pi),1) ; 
    
    idx = 0 ; 
    while idx(end)<T(trl)
        
        % Duration
        idx = idx(end)+(1:randi([Tmin,Tmax])) ; 
        idx(idx>T(trl)) = [] ; 
        
        % Fill in sequence
        fullidx = sum(T(1:(trl-1)))+idx ; 
        G(fullidx) = state ; 
        
        % Generate time series
        for fi = 1:length(f)
            X(fullidx,:) = X(fullidx,:) + ...
                A{state}(:,fi)'.*sin(2*pi*f(fi).*idx'/Fs + Phi{state}(:,fi)') ; 
        end
        
        % Get next state
        Pi = P(state,:) ; 
        state = find(rand < cumsum(Pi),1) ; 
    end
        
end

%% HMM estimation

embeddedlag = floor(60e-3*Fs/2) ; % 60 ms window, rounded to the nearest odd number of samples (then subtract time zero and divide by 2)
hmmopts = struct;
hmmopts.order = 0;
hmmopts.zeromean = 1;
hmmopts.covtype = 'full';
hmmopts.embeddedlags = -embeddedlag:embeddedlag;
hmmopts.pca = size(X,2)*2;
hmmopts.K = K;
hmmopts.Fs = Fs;
hmmopts.verbose = 1;
hmmopts.onpower = 0; 
hmmopts.standardise = 1;
hmmopts.standardise_pc = hmmopts.standardise;

[hmmest,Gamma,~,vpath] = hmmmar(X,T,hmmopts) ; 

%% Microstate estimation

% Make a +microstate cohort
coh = microstate.cohort ; 
for trl = 1:NTrl
    
    idx = sum(T(1:(trl-1))) + (1:T(trl)) ; 
    ms = microstate.individual(X(idx,:),'source',Fs) ; 
    coh = coh.add_individuals(ms,[],1); 
    
end
coh_pca = coh ; 
coh_ica = coh ; 

% Microstate
coh = coh.cluster_global(K) ;
coh = coh.cluster_globalmaps2individual ;
for trl = 1:NTrl
    coh.individual(trl).data = [] ; 
end

% PCA
coh_pca = coh_pca.cluster_global(K,'clustermethod','pca'); 
coh_pca = coh_pca.cluster_globalmaps2individual ; 
for trl = 1:NTrl
    coh_pca.individual(trl).data = [] ; 
end

% ICA
coh_ica = coh_ica.cluster_global(K,'clustermethod','ica') ; 
coh_ica = coh_ica.cluster_globalmaps2individual ; 
for trl = 1:NTrl
    coh_ica.individual(trl).data = [] ; 
end

% Split the path and data into trials
G0 = G ; X0 = X ; 
for trl = 1:NTrl
    pathtrue{trl} = G0(1:T(trl))' ; 
    G0(1:T(trl)) = [] ; 
    
    data{trl} = X0(1:T(trl),:) ; 
    X0(1:T(trl),:) = [] ; 
end

% HMM path
vpath0 = vpath ; 
for trl = 1:NTrl
    pathHMM{trl} = [nan(1,embeddedlag) , ...
        vpath0(1:(T(trl)-2*embeddedlag))' , ...
        nan(1,embeddedlag)] ; 
    
    vpath0(1:(T(trl)-2*embeddedlag)) = [] ; 
end

% Do MI
for trl = 1:NTrl
    
    ms = microstate.individual ; 
    ms.label = coh.individual(trl).label ; 
    ms = ms.stats_mutualinformation(pathtrue{trl}) ; 
    MIms(trl) = ms.stats.mutualinformation; 
    
    ms = microstate.individual ; 
    ms.label = coh_pca.individual(trl).label ; 
    ms = ms.stats_mutualinformation(pathtrue{trl}) ; 
    MIpca(trl) = ms.stats.mutualinformation; 
    
    ms = microstate.individual ; 
    ms.label = coh_ica.individual(trl).label ; 
    ms = ms.stats_mutualinformation(pathtrue{trl}) ; 
    MIica(trl) = ms.stats.mutualinformation; 
        
    ms = microstate.individual ; 
    ms.label = pathHMM{trl} ; 
    ms = ms.stats_mutualinformation(pathtrue{trl}) ; 
    MIhmm(trl) = ms.stats.mutualinformation; 
    
end

MIms
MIpca
MIica
MIhmm