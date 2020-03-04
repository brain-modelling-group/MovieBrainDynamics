function hmm = hmmhsinit(hmm,GammaInit,T)
% Initialise variables related to the Markov chain
%
% hmm		hmm data structure
%
% OUTPUT
% hmm           hmm structure
%
% Author: Diego Vidaurre, OHBA, University of Oxford

% if isfield(hmm.train,'grouping')
%     Q = length(unique(hmm.train.grouping));
% else
%     Q = 1;
% end
Q = 1; 

% define P-priors
defhmmprior=struct('Dir2d_alpha',[],'Dir_alpha',[]);
defhmmprior.Dir_alpha = ones(1,hmm.K);
defhmmprior.Dir_alpha(~hmm.train.Pistructure) = 0;
defhmmprior.Dir2d_alpha = ones(hmm.K);
defhmmprior.Dir2d_alpha(~hmm.train.Pstructure) = 0;
for k=1:hmm.K
    defhmmprior.Dir2d_alpha(k,k) = hmm.train.DirichletDiag;
end
defhmmprior.Dir2d_alpha=hmm.train.PriorWeighting.*defhmmprior.Dir2d_alpha;
% assigning default priors for hidden states
if ~isfield(hmm,'prior')
    hmm.prior=defhmmprior;
else
    % priors not specified are set to default
    hmmpriorlist = fieldnames(defhmmprior);
    fldname = fieldnames(hmm.prior);
    misfldname = find(~ismember(hmmpriorlist,fldname));
    for i = 1:length(misfldname)
        priorval = getfield(defhmmprior,hmmpriorlist{i});
        hmm.prior = setfield(hmm.prior,hmmpriorlist{i},priorval);
    end
end

if nargin > 1 && ~isempty(GammaInit)
    hmm = hsupdate([],GammaInit,T,hmm);
else
    % Initial state
    kk = hmm.train.Pistructure;
    if Q==1
        hmm.Dir_alpha = zeros(1,hmm.K);
        hmm.Dir_alpha(kk) = 1;
        hmm.Pi = zeros(1,hmm.K);
        hmm.Pi(kk) = ones(1,sum(kk)) / sum(kk);
    else
        hmm.Dir_alpha = zeros(hmm.K,Q);
        hmm.Dir_alpha(kk,:) = 1;
        hmm.Pi = zeros(hmm.K,Q);
        hmm.Pi(kk,:) = ones(sum(kk),Q) / sum(kk);
    end
    % State transitions
    hmm.Dir2d_alpha = zeros(hmm.K,hmm.K,Q);
    hmm.P = zeros(hmm.K,hmm.K,Q);
    for i = 1:Q
        for k = 1:hmm.K
            kk = hmm.train.Pstructure(k,:);
            hmm.Dir2d_alpha(k,kk,i) = 1;
            hmm.Dir2d_alpha(k,k,i) = hmm.train.DirichletDiag;
            hmm.Dir2d_alpha(k,kk,i) = hmm.train.PriorWeighting .* hmm.Dir2d_alpha(k,kk,i);
            hmm.P(k,kk,i) = hmm.Dir2d_alpha(k,kk,i) ./ sum(hmm.Dir2d_alpha(k,kk,i));
        end
    end
end

end
