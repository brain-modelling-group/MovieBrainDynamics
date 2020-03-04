function [hmm,XW] = updateW(hmm,Gamma,residuals,XX,XXGXX,Tfactor)

K = length(hmm.state); ndim = hmm.train.ndim;
if ~isempty(hmm.state(1).W.Mu_W)
    XW = zeros(size(XX,1),ndim,K);
else
    XW = [];
end
if nargin<6, Tfactor = 1; end
if isfield(hmm.train,'B'), Q = size(hmm.train.B,2);
else Q = ndim; end
pcapred = hmm.train.pcapred>0;
if pcapred, M = hmm.train.pcapred; end

for k=1:K
    if ~hmm.train.active(k), continue; end
    setstateoptions;
    if isempty(orders) && train.zeromean, continue; end
    if strcmp(train.covtype,'diag') || strcmp(train.covtype,'full'), omega = hmm.state(k).Omega;
    else omega = hmm.Omega;
    end
    if train.uniqueAR || ndim==1 % it is assumed that order>0 and cov matrix is diagonal
        if hmm.train.pcapred>0, npred = hmm.train.pcapred;
        else npred = length(orders);
        end
        XY = zeros(npred+(~train.zeromean),1);
        XGX = zeros(npred+(~train.zeromean));
        for n=1:ndim
            ind = n:ndim:size(XX,2);
            iomegan = omega.Gam_shape / omega.Gam_rate(n);
            XGX = XGX + iomegan * XXGXX{k}(ind,ind);
            XY = XY + (iomegan * XX(:,ind)' .* repmat(Gamma(:,k)',length(ind),1)) * residuals(:,n);
        end
        if ~isempty(train.prior)
            hmm.state(k).W.S_W = inv(train.prior.iS + XGX);
            hmm.state(k).W.Mu_W = hmm.state(k).W.S_W * (XY + train.prior.iSMu); % order by 1
        else
            if train.zeromean==0 && pcapred
                regterm = diag([hmm.state(k).prior.Mean.iS; (hmm.state(k).beta.Gam_shape ./ ...
                    hmm.state(k).beta.Gam_rate) ]);  
            elseif pcapred
                regterm = diag((hmm.state(k).beta.Gam_shape ./  hmm.state(k).beta.Gam_rate));
            elseif train.zeromean==0 && ~isempty(orders)
                regterm = diag([hmm.state(k).prior.Mean.iS (hmm.state(k).alpha.Gam_shape ./ ...
                    hmm.state(k).alpha.Gam_rate) ]);
            elseif train.zeromean==0
                regterm = diag(hmm.state(k).prior.Mean.iS);
            else
                regterm = diag((hmm.state(k).alpha.Gam_shape ./  hmm.state(k).alpha.Gam_rate));
            end
            hmm.state(k).W.S_W = inv(regterm + Tfactor * XGX);
            hmm.state(k).W.Mu_W = Tfactor * hmm.state(k).W.S_W * XY; % order by 1
        end        
        for n=1:ndim
            ind = n:ndim:size(XX,2);
            XW(:,n,k) = XX(:,ind) * hmm.state(k).W.Mu_W;
        end
        
    elseif strcmp(train.covtype,'diag') || strcmp(train.covtype,'uniquediag')
        for n=1:ndim
            ndim_n = sum(S(:,n)>0);
            if ndim_n==0 && train.zeromean==1, continue; end
            regterm = [];
            if ~train.zeromean, regterm = hmm.state(k).prior.Mean.iS(n); end
            if ~isempty(orders)
                if pcapred
                    regterm = [regterm; hmm.state(k).beta.Gam_shape(:,n) ./ hmm.state(k).beta.Gam_rate(:,n)];
                else
                    alphaterm = repmat( (hmm.state(k).alpha.Gam_shape ./  hmm.state(k).alpha.Gam_rate), ndim_n, 1);
                    if ndim>1
                        regterm = [regterm; repmat(hmm.state(k).sigma.Gam_shape(S(:,n),n) ./ ...
                            hmm.state(k).sigma.Gam_rate(S(:,n),n), length(orders), 1).*alphaterm(:) ];
                    else
                        regterm = [regterm; alphaterm(:)];
                    end
                end
            end
            if isempty(regterm), regterm = 0; end
            regterm = diag(regterm);
            hmm.state(k).W.iS_W(n,Sind(:,n),Sind(:,n)) = ...
                regterm + Tfactor * (omega.Gam_shape / omega.Gam_rate(n)) * XXGXX{k}(Sind(:,n),Sind(:,n));
            hmm.state(k).W.S_W(n,Sind(:,n),Sind(:,n)) = ...
                inv(permute(hmm.state(k).W.iS_W(n,Sind(:,n),Sind(:,n)),[2 3 1]));
            hmm.state(k).W.Mu_W(Sind(:,n),n) = (( permute(hmm.state(k).W.S_W(n,Sind(:,n),Sind(:,n)),[2 3 1]) * ...
                Tfactor * (omega.Gam_shape / omega.Gam_rate(n)) * XX(:,Sind(:,n))') .* ...
                repmat(Gamma(:,k)',sum(Sind(:,n)),1)) * residuals(:,n);
        end
        XW(:,:,k) = XX * hmm.state(k).W.Mu_W;
        
    else % full or unique full - this only works if all(S(:)==1); any(S(:)~=1) is just not yet implemented 
        if pcapred
            mlW = (( XXGXX{k} \ XX') .* repmat(Gamma(:,k)',(~train.zeromean)+M,1) * residuals)';
        else
            mlW = (( XXGXX{k} \ XX') .* repmat(Gamma(:,k)',...
                (~train.zeromean)+Q*length(orders),1) * residuals)';
        end
        regterm = [];
        if ~train.zeromean, regterm = hmm.state(k).prior.Mean.iS; end % ndim by 1
        if ~isempty(orders) 
            if pcapred
                betaterm = (hmm.state(k).beta.Gam_shape ./ hmm.state(k).beta.Gam_rate)';
                regterm = [regterm; betaterm(:)];
            else
                sigmaterm = (hmm.state(k).sigma.Gam_shape ./ hmm.state(k).sigma.Gam_rate)'; 
                sigmaterm = sigmaterm(:); 
                sigmaterm = repmat(sigmaterm, length(orders), 1); % ndim*ndim*order by 1 
                alphaterm = repmat( (hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate), ...
                    length(hmm.state(k).sigma.Gam_rate(:)), 1);
                alphaterm = alphaterm(:);
                regterm = [regterm; (alphaterm .* sigmaterm)];
            end
        end
        if isempty(regterm), regterm = 0; end
        regterm = diag(regterm);
        prec = omega.Gam_shape * omega.Gam_irate;
        gram = kron(XXGXX{k}, prec);
        hmm.state(k).W.iS_W = regterm + Tfactor * gram;
        hmm.state(k).W.S_W = inv(hmm.state(k).W.iS_W);
        muW = Tfactor * hmm.state(k).W.S_W * gram * mlW(:);
        if pcapred
            hmm.state(k).W.Mu_W = reshape(muW,ndim,(~train.zeromean)+M)';
        else
            hmm.state(k).W.Mu_W = reshape(muW,ndim,~train.zeromean+Q*length(orders))';
        end
        XW(:,:,k) = XX * hmm.state(k).W.Mu_W;
    end
    
end

end