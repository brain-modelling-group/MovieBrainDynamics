function hmm = updateOmega(hmm,Gamma,Gammasum,residuals,T,XX,XXGXX,XW,Tfactor)

is_gaussian = hmm.train.order == 0; % true if Gaussian observation model is being used

K = length(hmm.state); ndim = hmm.train.ndim;
S = hmm.train.S==1; regressed = sum(S,1)>0;
if nargin<9, Tfactor = 1; end
Tres = sum(T) - length(T)*hmm.train.maxorder;
if isfield(hmm.train,'B'), Q = size(hmm.train.B,2);
else Q = ndim; end

if strcmp(hmm.train.covtype,'uniquediag') && hmm.train.uniqueAR 
    % all are AR and there's a single covariance matrix
    hmm.Omega.Gam_rate(k) = hmm.prior.Omega.Gam_rate;
    for k=1:K
        if ~isempty(XW)
            XWk = XW(:,:,k);
        else
            XWk = zeros(size(residuals));
        end
        e = (residuals - XWk).^2;
        swx2 = zeros(Tres,ndim);
        for n=1:ndim
            ind = n:ndim:size(XX,2);
            tmp = XX(:,ind) * hmm.state(k).W.S_W;
            swx2(:,n) = sum(tmp .* XX(:,ind),2);
        end
        hmm.Omega.Gam_rate = hmm.Omega.Gam_rate + ...
            0.5 * Tfactor * sum( repmat(Gamma(:,k),1,ndim) .* (e + swx2) );
    end
    hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + 0.5 * Tfactor * Tres;    
    
elseif strcmp(hmm.train.covtype,'uniquediag')
    hmm.Omega.Gam_rate(regressed) = hmm.prior.Omega.Gam_rate(regressed);
    for k=1:K
        setstateoptions;
        if ~isempty(XW)
            XWk = XW(:,:,k);
        else
            XWk = zeros(size(residuals));
        end
        e = (residuals(:,regressed) - XWk(:,regressed)).^2;
        swx2 = zeros(Tres,ndim);
        if ~isempty(hmm.state(k).W.Mu_W)
            for n=1:ndim
                if ~regressed(n), continue; end
                if ndim==1
                    tmp = XX(:,Sind(:,n)) * hmm.state(k).W.S_W;
                else
                    tmp = XX(:,Sind(:,n)) * permute(hmm.state(k).W.S_W(n,Sind(:,n),Sind(:,n)),[2 3 1]);
                end
                swx2(:,n) = sum(tmp .* XX(:,Sind(:,n)),2);
            end
        end
        hmm.Omega.Gam_rate(regressed) = hmm.Omega.Gam_rate(regressed) + ...
            0.5 * Tfactor * sum( repmat(Gamma(:,k),1,sum(regressed)) .* (e + swx2(:,regressed) ) );
    end
    hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + 0.5 * Tfactor * Tres;
    
elseif strcmp(hmm.train.covtype,'uniquefull')
    hmm.Omega.Gam_rate(regressed,regressed) = hmm.prior.Omega.Gam_rate(regressed,regressed);
    for k=1:K
        setstateoptions;
        if ~isempty(XW)
            XWk = XW(:,:,k);
        else
            XWk = zeros(size(residuals));
        end
        e = (residuals(:,regressed) - XWk(:,regressed));
        e = (e' .* repmat(Gamma(:,k)',sum(regressed),1)) * e;
        swx2 =  zeros(ndim,ndim);
        if ~isempty(hmm.state(k).W.Mu_W)  
            for n1=find(regressed)
                for n2=find(regressed)
                    if n2<n1, continue, end
                    if hmm.train.pcapred>0
                        index1 = (0:hmm.train.pcapred+(~hmm.train.zeromean)-1) * ndim + n1;
                        index2 = (0:hmm.train.pcapred+(~hmm.train.zeromean)-1) * ndim + n2;
                    else
                        index1 = (0:length(orders)*Q+(~hmm.train.zeromean)-1) * ndim + n1;
                        index2 = (0:length(orders)*Q+(~hmm.train.zeromean)-1) * ndim + n2;
                        index1 = index1(Sind(:,n1)); index2 = index2(Sind(:,n2));
                    end
                    swx2(n1,n2) = sum(sum(hmm.state(k).W.S_W(index1,index2) .* XXGXX{k}(Sind(:,n1),Sind(:,n2))));
                    swx2(n2,n1) = swx2(n1,n2);
                end
            end
        end
        hmm.Omega.Gam_rate(regressed,regressed) = hmm.Omega.Gam_rate(regressed,regressed) + Tfactor * (e + swx2(regressed,regressed));
    end
    hmm.Omega.Gam_irate(regressed,regressed) = inv(hmm.Omega.Gam_rate(regressed,regressed));
    hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tfactor * Tres;
    
else % state dependent
    for k=1:K
        if ~hmm.train.active, continue; end
        setstateoptions;
        if ~isempty(XW)
            XWk = XW(:,:,k);
        else 
            XWk = zeros(size(residuals));
        end
        if train.uniqueAR
            e = (residuals - XWk).^2;
            swx2 = zeros(Tres,ndim);
            for n=1:ndim
                ind = n:ndim:size(XX,2);
                swx2(:,n) = sum(XX(:,ind) * hmm.state(k).W.S_W .* XX(:,ind),2);
            end
            hmm.state(k).Omega.Gam_rate = hmm.state(k).prior.Omega.Gam_rate + ...
                0.5 * Tfactor * sum( repmat(Gamma(:,k),1,ndim) .* (e + swx2) );
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + 0.5 * Tfactor * Gammasum(k);
                
        elseif strcmp(train.covtype,'diag')
            e = (residuals(:,regressed) - XWk(:,regressed)).^2;
            swx2 = zeros(Tres,ndim);
            if ~isempty(hmm.state(k).W.Mu_W)  
                for n=1:ndim
                    if ~regressed(n), continue; end
                    if ndim==1
                        swx2(:,n) = sum(XX(:,Sind(:,n)) * hmm.state(k).W.S_W(Sind(:,n),Sind(:,n)) ...
                            .* XX(:,Sind(:,n)),2);
                    else
                        swx2(:,n) = sum(XX(:,Sind(:,n)) * permute(hmm.state(k).W.S_W(n,Sind(:,n),Sind(:,n)),[2 3 1]) ...
                            .* XX(:,Sind(:,n)),2);
                    end
                end
            end
            hmm.state(k).Omega.Gam_rate(regressed) = hmm.state(k).prior.Omega.Gam_rate(regressed) + ...
                0.5 * Tfactor * sum( repmat(Gamma(:,k),1,sum(regressed)) .* (e + swx2(:,regressed) ) );
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + 0.5 * Tfactor * Gammasum(k);
            
        else % full
            if is_gaussian
                if hmm.train.zeromean % If zeromean == 1 for Gaussian model, then XWk is zero and we don't need to subtract at all
                    e = residuals(:,regressed);
                else % If zeromean == 0 for Gaussian model, then XWk has all the same rows, and bsxfun() is a fair bit faster
                    e = bsxfun(@minus,residuals(:,regressed),XWk(1,regressed));
                end
            else
                e = (residuals(:,regressed) - XWk(:,regressed));
            end
            e = (bsxfun(@times,e,Gamma(:,k)))' * e;
            swx2 =  zeros(ndim,ndim);
            if ~isempty(hmm.state(k).W.Mu_W)  
                orders = formorders(train.order,train.orderoffset,train.timelag,train.exptimelag);
                if isempty(orders)
                    swx2 = hmm.state(k).W.S_W * XXGXX{k};
                else
                    for n1=find(regressed)
                        for n2=find(regressed)
                            if n2<n1, continue, end
                            if hmm.train.pcapred>0
                                index1 = (0:hmm.train.pcapred+(~train.zeromean)-1) * ndim + n1;
                                index2 = (0:hmm.train.pcapred+(~train.zeromean)-1) * ndim + n2;
                            else
                                index1 = (0:length(orders)*Q+(~train.zeromean)-1) * ndim + n1;
                                index2 = (0:length(orders)*Q+(~train.zeromean)-1) * ndim + n2;
                                index1 = index1(Sind(:,n1)); index2 = index2(Sind(:,n2));
                            end
                            swx2(n1,n2) = sum(sum(hmm.state(k).W.S_W(index1,index2) .* XXGXX{k}(Sind(:,n1),Sind(:,n2))));
                            swx2(n2,n1) = swx2(n1,n2);
                        end
                    end
                end
            end
            if hmm.train.firsteigv
                [V,D] = svd(Tfactor * e);
                hmm.state(k).Omega.Gam_rate = sum(diag(D)) * (V(:,1) * V(:,1)') + hmm.state(k).prior.Omega.Gam_rate;
            else
                hmm.state(k).Omega.Gam_rate(regressed,regressed) = hmm.state(k).prior.Omega.Gam_rate(regressed,regressed) + ...
                    Tfactor * (e + swx2(regressed,regressed));
            end
            hmm.state(k).Omega.Gam_irate(regressed,regressed) = inv(hmm.state(k).Omega.Gam_rate(regressed,regressed));
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + Tfactor * Gammasum(k);
        end
    end
    
end

end
