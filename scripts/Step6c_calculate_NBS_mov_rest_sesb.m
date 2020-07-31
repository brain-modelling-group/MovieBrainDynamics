
clear all
close all
% do NBS on the shit
addpath('NBSdirected');



load('outr.mat');
ANALYSIS='all';

load valid_inferences_all.mat
RUN=valid_inferences(1); % you might change this depending on which HMM inference you wish to check out.


% RUN=10;

COMP_A = 'resta';
COMP_B = 'mov1a';

run=sprintf('run%d',RUN); 
% subs_to_use_real = [2     3     7     8     9    10    11    12    14    16    17    18    19    20];



mats_a = outr.(run).aroma.(ANALYSIS).(COMP_A).emp(:, :, :);
mats_b = outr.(run).aroma.(ANALYSIS).(COMP_B).emp(:, :, :);



C=cat(3,mats_a,mats_b);

nsubs=size(mats_a,3);

X=[ones(nsubs,1) zeros(nsubs,1); zeros(nsubs,1) ones(nsubs,1)];





%Any significant results are stored in variable called con_mat

%Connectivity matrices (regions x regions x subjects)
% C=randn(10,10,6);
% C(1,2,1:3)=C(1,2,1:3)+10;
% C(1,3,1:3)=C(1,3,1:3)+10;
% C(1,4,1:3)=C(1,4,1:3)+10;
% C(4,3,1:3)=C(4,3,1:3)+10;
% C(4,1,1:3)=C(4,3,1:3)+10;

%Total number of permutations to generate
GLM.perms=5000;
%Design matrix
GLM.X=X;
 %Contrast
 GLM.contrast=[1 -1]; 
 %Type of test
 GLM.test='ttest'; % 'ttest' or 'ftest'
 %Exchange block for repeated measures
 %GLM.exchange=[]; STAT
 STATS.size='Extent'; %'Intensity' or 'Extent'
 %Threshold
 STATS.thresh=3.1; 
 %Significance (usually 0.05)
 STATS.alpha=0.05; 

 %NO NEED TO CHANGE ANYTHING BEYOND THIS POINT
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 bgl=0; 
 nDisp=100; 
 
 %Number of nodes
 N=size(C,1);
 ind_uplo=union(find(triu(ones(N,N),1)),find(tril(ones(N,N),-1)));
 
 GLM.y=zeros(size(GLM.X,1),length(ind_uplo));
 for i=1:size(C,3)
     tmp=C(:,:,i);
     GLM.y(i,:)=tmp(ind_uplo);
 end

%Precompute test stat
STATS.test_stat=NBSglm(GLM); 

%Number of edges
J=length(ind_uplo); 

%Determine whether test statistics have been precomputed and determine
%index of edges exceeding the primary threshold
if ~isempty(STATS.test_stat)
    %Precomputed test statistics
    ind=ind_uplo(STATS.test_stat(1,:)>STATS.thresh); 
    %Number of permutations
    K=size(STATS.test_stat,1)-1; 
else
    %Never get to this case
end

%Size of a component measured using extent or intensity? 
Intensity=0;
if strcmp(STATS.size,'Intensity')    
    %If size measure using intensity, create an N x N matrix cotaining the 
    %test statistic for each edge minus the test statistic threshold
    %(primary threshold)
    Intensity=1; 
    %Compute a test statistic matrix
    test_stat_mat=zeros(N,N); 
    if ~isempty(STATS.test_stat)
        %Precomputed
        test_stat_mat(ind_uplo)=STATS.test_stat(1,:)-STATS.thresh;
%       test_stat_mat=(test_stat_mat+test_stat_mat');
    else
        %Never reach this case. 
    end
end
    
adj=spalloc(N,N,length(ind));
adj(ind)=1; 
%Only consider components comprising more than one node, equivalent to at
%least one edge
if bgl==1
    [a,sz]=components((adj+adj')/2); 
else
    [a,sz]=get_components((adj+adj')/2); 
end
ind_sz=find(sz>1);
sz_links=zeros(1,length(ind_sz));
max_sz=0; 
for i=1:length(ind_sz)
    nodes=find(ind_sz(i)==a);
    if Intensity
        %Measure size as intensity
        sz_links(i)=sum(sum(adj(nodes,nodes).*test_stat_mat(nodes,nodes))); %/2;
    else
        %Measure size as extent
        sz_links(i)=sum(sum(adj(nodes,nodes))); %/2;
    end
    adj(nodes,nodes)=adj(nodes,nodes)*(i+1);
    if max_sz<sz_links(i)
        max_sz=sz_links(i);
    end
end

%Subtract one to remove edges not part of a component
%Although one is also subtracted from edges comprising a component, this is 
%compensated by the (i+1) above
adj(~~adj)=adj(~~adj)-1;

%Repeat above for each permutation
%Empirical null distribution of maximum component size
null_dist=zeros(K,1); 
str1='| Permutation | Max Size | Max Size | Lowest  |';
str2='|             |  Random  |  Actual  | p-value |';
try tmp=get(H,'string'); set(H,'string',[{str1};{str2};tmp]); drawnow;
catch;  fprintf([str1,'\n',str2,'\n']); end 
p_approx=0;
%Store what is already displayed in the listbox
try pre_str=get(H,'string'); catch; end
new_str={};
%First row of test_stat is the observed test statistics, so start at the
%second row
for i=2:K+1
    if ~isempty(STATS.test_stat)
        %Precomputed test statistics 
        ind=ind_uplo(STATS.test_stat(i,:)>STATS.thresh); 
    else

    end
    if Intensity 
        %Compute a test statistic matrix
        test_stat_mat=zeros(N,N); 
        if ~isempty(STATS.test_stat)
            test_stat_mat(ind_uplo)=STATS.test_stat(i,:)-STATS.thresh;
        else

        end    
    end
    adj_perm=spalloc(N,N,length(ind));
    adj_perm(ind)=1;
    if bgl==1
        [a,sz]=components((adj_perm+adj_perm')/2); 
    else
        [a,sz]=get_components((adj_perm+adj_perm')/2); 
    end
    ind_sz=find(sz>1);
    max_sz_perm=0; 
    for j=1:length(ind_sz)
        nodes=find(ind_sz(j)==a);
        if Intensity
            tmp=sum(sum(adj_perm(nodes,nodes).*test_stat_mat(nodes,nodes))); %/2;
        else
            tmp=sum(sum(adj_perm(nodes,nodes))); %/2;
        end
        if tmp>max_sz_perm
            max_sz_perm=full(tmp);
        end   
    end
    null_dist(i-1)=max_sz_perm; 
    if max_sz_perm>=max_sz
        p_approx=p_approx+1;
    end
   % str=sprintf('|   %5d/%5d |     %4d |     %4d |   %0.3f |',...
   %v1.1.2 Changed to %6.0f to %6.1f to allow fractional component sizes
   %that arise when component size is measured with intensity. 
   str=sprintf('| %5d/%5d |   %6.1f |   %6.1f |  %0.4f |',...
            i-1,K,max_sz_perm,max_sz,p_approx/(i-1));
        %Display no mare than nDisp most recent permutations
        new_str=[str,{new_str{1:min(nDisp,length(new_str))}}]';
        try set(H,'string',[new_str;pre_str]); drawnow; 
        catch;  fprintf([str,'\n']); end 
end
str1='| Permutation | Max Size | Max Size | Lowest  |';
str2='|             |  Random  |  Actual  | p-value |';
try tmp=get(H,'string'); set(H,'string',[{str1};{str2};tmp]); drawnow;
catch;  fprintf([str1,'\n',str2,'\n']); end 


%Determine components satisfying alpha significance threshold
n_cnt=0; 
for i=1:length(sz_links)
    tmp=sum(null_dist>=sz_links(i))/K;
    if tmp<=STATS.alpha
        n_cnt=n_cnt+1;
        ind=find(adj==i);
        con_mat{n_cnt}=spalloc(N,N,length(ind)*2);
        con_mat{n_cnt}(ind)=1; 
        con_mat{n_cnt}=con_mat{n_cnt};
        pval(n_cnt)=tmp;
    end
end
if n_cnt==0
    pval=[]; con_mat=[]; 
end


if numel(con_mat)>0
    m=full(con_mat{1});
    if numel(con_mat) > 1
        for i=2:numel(con_mat)
            m = m + con_mat{i};
        end
    end
else
    m = zeros(10, 10);
end

save Step6c_calculate_NBS_mov_rest_matrix_sesb.mat m


