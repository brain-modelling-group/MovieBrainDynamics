function [fo,dign,avg_life] = summary_measures( st_path,n_sub,do )
%this function calculates empoirical probability for every subject
% Inputs: 
%         st_path = 1 x n_sub*n_time vector as outputted from fit_hmm
%         function
%           n_sub = Number of subjects
%              do = Flag to replace 0s with NaNs
%                   do=0 keeps all misisng values as zeros; do=1 sets them
%                   to Nans
% Output:
%        fo = Fractional occupancies (in %) calculated for
%             every subject  in each state. Size = n_sub X n_state
%      dign = Diagonal probabilites (probability to remain in a given state)
%             calculated for every subject  in each state. 
%             Size = n_sub X n_state
%  avg_life = Time spent in a state when visited (on average)
%             Size = n_sub X n_state

S=max(st_path); % Number of states 
n_time=length(st_path)/n_sub;

%%%%%%%%%% FO
fo=zeros(n_sub,S);
for sub=1:n_sub
    for s=1:S
    fo(sub,s)=nnz(st_path((sub-1)*n_time+1:sub*n_time)==s);
    end
end
fo=100*fo/n_time;

%%%%%%%%%%%%%%%% Emp. Prob
subj_emp=zeros(S,S,n_sub);
for sub=1:n_sub
    emp=zeros(S,S);
    stat=squeeze(st_path((sub-1)*n_time+1:sub*n_time));
    for j=1:S
        for k=1:S
            for q=2:length(stat)
                if stat(q)==k && stat(q-1)==j
                    emp(j,k)=emp(j,k)+1;
                end
            end
        end
    end
    for r=1:S
        emp(r,:)=emp(r,:)/sum(emp(r,:)); % to make sure every row adds to 1
        % i.e., total probability = 1
    end
    emp(isnan(emp))=0; % removing any NaNs due to division by zero cases
    % from above step
    
   % if a row is all zero, that means the subject never visited that state.
   % Here the diagonal prob is set to 1, just to make sure that evry row
   % adds upto 1
   % Uncomment the lines below if you want to set diagonals to 1 in the
   % case of states not visited
%     row_has_all_zeros = ~any(emp, 2);
%     indices = find(row_has_all_zeros);
%     for k=1:length(indices)
%         emp(indices(k),indices(k))=1;
%     end
    subj_emp(:,:,sub)=emp;
end
%%%%%%%%%%%%%%%%%%% Diagonal Prob.
dign=zeros(size(fo));
for sub=1:n_sub
    % keyboard;
    %tmp=squeeze(subj_emp(sub,:,:));
    tmp=squeeze(subj_emp(:,:,sub));
    dign(sub,:)=diag(tmp);
end
if do
    for k=1:S
        tmp=dign(:,k); 
        loc=find(tmp==0); %cnt_t1(k)=nnz(loc); % find subjects with missing
%         cnt_p(k)=nnz(loc>41);cnt_c(k)=nnz(loc<42);
        dign(loc,k)=NaN;fo(loc,k)=NaN;
    end
end
%%%%%%%%%%%%%%%%%%%%%% 
life_time=cell(n_sub,S); avg_life=zeros(n_sub,S);
for sub=1:n_sub
    sub_st=st_path((sub-1)*n_time+1:sub*n_time);
    for k=1:S
        loc=(sub_st==k);
        try
            cnt=diff([0;find(diff(loc));numel(loc)]);
        catch
            keyboard;
        end
        if loc(1)==1
            life_time{sub,k}=cnt(1:2:end);  
        else
            life_time{sub,k}=cnt(2:2:end);
        end
        avg_life(sub,k)=mean(life_time{sub,k});
    end
end

for k=1:S
    tmp=avg_life(:,k);
    loc=find(isnan(tmp)); 
    if do
%         cnt_t2(k)=nnz(loc); % no of subjects thrown away if put NaN
        avg_life(loc,k)=NaN;fo(loc,k)=NaN;
    else
        avg_life(loc,k)=0;
    end
end
