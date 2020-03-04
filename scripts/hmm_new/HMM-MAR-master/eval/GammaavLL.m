function avLL = GammaavLL(hmm,Gamma,Xi,T)
% average loglikelihood for state time course
% if isfield(hmm.train,'grouping')
%     Q = length(unique(hmm.train.grouping));
% else
%     Q = 1;
% end
Q = 1; 
N = length(T); 
order = (sum(T) - size(Gamma,1))/N;
avLL = 0; K = size(Gamma,2);
jj = zeros(N,1); % reference to first time point of the segments
for in = 1:N
    jj(in) = sum(T(1:in-1)) - order*(in-1) + 1;
end
% avLL initial state
% if Q>1
%     for i = 1:Q
%         PsiDir_alphasum = psi(sum(hmm.Dir_alpha(:,i)));
%         ii = hmm.train.grouping==i;
%         for l = 1:K
%             if ~hmm.train.Pistructure(l), continue; end
%             avLL = avLL + sum(Gamma(jj(ii),l)) * (psi(hmm.Dir_alpha(l,i)) - PsiDir_alphasum);
%         end
%     end
% else
%     PsiDir_alphasum = psi(sum(hmm.Dir_alpha));
%     for l = 1:K 
%         if ~hmm.train.Pistructure(l), continue; end
%         avLL = avLL + sum(Gamma(jj,l)) * (psi(hmm.Dir_alpha(l)) - PsiDir_alphasum);
%     end
% end
PsiDir_alphasum = psi(sum(hmm.Dir_alpha));
for l = 1:K
    if ~hmm.train.Pistructure(l), continue; end
    avLL = avLL + sum(Gamma(jj,l)) * (psi(hmm.Dir_alpha(l)) - PsiDir_alphasum);
end
% avLL remaining time points
for i = 1:Q
    if Q > 1
        ii = find(hmm.train.grouping==i)';
    else
        ii = 1:length(T);
    end
    PsiDir2d_alphasum = zeros(K,1);
    for l = 1:K, PsiDir2d_alphasum(l) = psi(sum(hmm.Dir2d_alpha(l,:,i))); end
    for k = 1:K
        for l = 1:K
            if ~hmm.train.Pstructure(l,k), continue; end
            if Q==1
                avLL = avLL + sum(Xi(:,l,k)) * (psi(hmm.Dir2d_alpha(l,k))-PsiDir2d_alphasum(l));
                if isnan(avLL)
                    error(['Error computing log likelihood of the state time courses  - ' ...
                        'Out of precision?'])
                end
            else
                for n = ii
                    t = (1:T(n)-1-order) + sum(T(1:n-1)) - (order+1)*(n-1) ;
                    avLL = avLL + sum(Xi(t,l,k)) * ...
                        (psi(hmm.Dir2d_alpha(l,k,i))-PsiDir2d_alphasum(l));
                end
            end
        end
    end
end

end
