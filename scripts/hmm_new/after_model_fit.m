clear;clc; tic;
reps = 5;
do = 1; % do = Flag to replace 0s with NaNs (in all summary measures)
        % do=0 keeps all misisng values as zeros; do=1 sets them
        % to Nans; this is important while comparing between groups!!
        
for r = 1:reps
    load(['HMMrun_rep_' num2str(r) '.mat'],'Gamma','vpath',...
          'hmm','J','T','K','n_sub'); 
      
    mean_em=zeros(J,K); % mean emissions
    for k=1:K
        mean_em(:,k)=getMean(hmm,k);
    end
    
    prob=hmm.P; % transition probabilities
    
    [FO,dign,avg_life]=summary_measures(vpath,n_sub,do);
    
    save(['Summary_measures_rep_',num2str(r) '.mat']);
end
    
    