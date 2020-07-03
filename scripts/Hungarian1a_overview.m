% we need to define the iterations where all 10 states were expressed

iters = [100 102 103 118 120 130 136 150 34 42 77 85 88 91 94]; % 15 iterations
load valid_inferences
iters=valid_inferences %(1:15);
NSTATES = 10;
% then we need to match each of these iteration's states to the other
% iteration's states.
clear mapping
mapping=struct();
mapping(:)=[];

for i=1:numel(iters)

    % load the states of i
    hmmvars_i=load(sprintf('../results_10/aroma/resta-b/HMMrun_rep_%d.mat', iters(i)));
    allw_i=[];  
    for ii=1:NSTATES
        allw_i=[allw_i hmmvars_i.hmm.state(ii).W.Mu_W']; 
    end
    
    mapping(i).i = i;
    mapping(i).iter_val = iters(i);
    
    mapping(i).assignmatrix=zeros(NSTATES, numel(iters));
    
    mapping(i).corrvals =zeros(NSTATES, numel(iters));
    mapping(i).assign_iter_vals = zeros(1, numel(iters));
    
    vpath_i = hmmvars_i.vpath;
    
    for j=setxor(1:numel(iters), i)
        
        % load the states of i
        hmmvars_j=load(sprintf('../results_10/aroma/resta-b/HMMrun_rep_%d.mat', iters(j)));
        
        vpath_j = hmmvars_j.vpath;
        allw_j=[];  
        for jj=1:NSTATES
            allw_j=[allw_j hmmvars_j.hmm.state(jj).W.Mu_W']; 
        end
        
        % bookkeeping
        mapping(i).assign_iter_vals(j) = iters(j);
        
        % do the Hungarian Algorithm to find closest match in the other
        % set; cost = 1-corr.
        
        % calculate the cost(s) associated iwth the 2 vpaths:
        JO=[];
        for iState=1:NSTATES
            for jState=1:NSTATES
                binned_i = vpath_i==iState;
                binned_j = vpath_j==jState;
                JO(iState, jState) = numel(intersect(find(binned_i),find(binned_j))) / numel(union(find(binned_i),find(binned_j)));
            end
        end
                
        
        [reassignments, costs] = munkres((1-corr(allw_i, allw_j)));
        % [reassignments, costs] = munkres((1-corr(allw_i, allw_j)) + (1-JO));
        % [reassignments, costs] = munkres((1-JO));
        
        % just check if I didn't get it wrong...
        % tmp3=sortrows([1:10; reassignments]',2);
        % alternative = tmp3(:,1);
        % reassignments = alternative;
        
        % bookeeping
        mapping(i).assignmatrix(:, j) = reassignments;
        
        % reshuffle the allw_j's
        reshuffled_all_j = allw_j(:, reassignments);
        
        % calculate correlation between allw_i and this allw_j:
        mapping(i).corrvals(:, j) = diag(corr(allw_i, reshuffled_all_j));
        
        
        mapping(i).vpaths(:, i) = vpath_i;  % we do not reshuffle the i
        
        % reshuffle the vpath according to reassignments:
        reshuffled_vpath_j = [];
        tmp3=sortrows([1:10; reassignments]',2);
        alternative = tmp3(:,1);
        for ir=1:numel(reassignments)
            reshuffled_vpath_j(vpath_j==ir) = alternative(ir);
        end
        
        mapping(i).vpaths(:, j) = reshuffled_vpath_j;
        mapping(i).unshuffled_vpaths(:, j) = vpath_j;
        
        
    end
    
    % average correlation with matched state in the other hmm inference:
    tmp=mapping(i).corrvals;
    
    
    mapping(i).avg_corr_with_others = mean(tmp);
    tmp2=mean(tmp);
    tmp2=mean(tmp2(tmp2~=0));
    mapping(i).overall_corr = tmp2;
    
end

fh=figure;lh=plot([mapping(:).overall_corr]);
set(lh,'linewidth',2, 'marker','+','color','k');
xlim([0 23]);
box off

xlabel('HMM Inference (iteration number)');
ylabel('Average correlation with other inferences');
title('Consistency of BOLD loadings of HMM inferences across Inferences');

saveas(gcf,'hmm_inferences_consistency.pdf');


save HungarianMapping.mat mapping
