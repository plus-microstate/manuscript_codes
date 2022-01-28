%% User input: List of subject IDs

% IDs of subjects to be analysed
global subject_IDs
subject_IDs = {'0001','0002'} ; % ,'0003','0004','0007','0008','0009','0014','0015','0016','0019','0021' , ...
             ...'0027','0032','0034','0035','0038','0039','0040','0041','0042','0043','0044','0045' , ...
            ...'0046','0047','0048','0049','0052','0053'} ; 
           
% Directory containing MEG data
global data_dir
data_dir = './source_meg_data' ; % Directory containing data

%% Main text sections
% % simulations() ; % Section 3.1 and Section S2.3
% % resting_state_microstate_maps() ; % Section 3.2 
% % statistics_of_resting_state_microstate_sequences() ; % Section 3.3
% % microstate_specific_functional_connectivity() ; % Section 3.4
% % microstate_response_to_auditory_stimulus() ; % Section 3.5

% Supplementary sections
% % reproducibility_and_reliability_of_maps() ; % Section S2.1
% % assessing_non_stationarity() ; % Section S2.2
% comparing_state_estimation_methods_phase(subject_IDs,dirs) ; % Section S2.3
% % GFP_vs_topographic_change() ; % Figure S4