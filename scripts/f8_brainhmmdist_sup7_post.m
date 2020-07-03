
% to be run after f8_ script...

% we have questions, the N-by-4 answers:

questions_dist_matrix = s(1).questions_dist_matrix;
% do the adjacency stuff, as per Luke Chang's suggestion:
questions_text={'bored','enjoy','emotional','audio'};

qD=zeros(max(size(questions_dist_matrix)),max(size(questions_dist_matrix)));
questions_dist_matrix(:,end+1)=zeros(size(questions_dist_matrix,1),1);

for i=1:size(questions_dist_matrix,1)
    for j=1:i
        qD(j,i) = questions_dist_matrix(i,j);
        qD(i,j) = questions_dist_matrix(i,j);
    end
end


% now we have the qD

% do the magic -- figure out the coords in 2D that preserves as much as
% possible the distances.
[Y,eigvals] = cmdscale(qD);


% inspect whether 2D or 3D is enough:
format short g
[eigvals eigvals/max(abs(eigvals))]



answers = {};
for i=1:size(questions,1)
    ans_txt = sprintf('%d,',questions(i,:));
    ans_txt(end)=[];
    answers{end+1} = [sprintf('P%2d  ',i) ans_txt];

end


% make the plot -- and also plot what were the answers:
figure;plot(Y(:,1),Y(:,2),'.','color','k','marker','o','MarkerFaceColor','k')
for i=1:numel(answers)
    text(Y(i,1)+0.05,Y(i,2),answers{i});
end
% title('questions: Bored?, Enjoyed?, Emotional?, Audio?');
 title('Interview Answers');

 

%%
% make the plot with the other distance matrices?
hmm_dist_matrix = s(1).hmm_dist_matrix;




qD=zeros(max(size(hmm_dist_matrix)),max(size(hmm_dist_matrix)));
hmm_dist_matrix(:,end+1)=zeros(size(hmm_dist_matrix,1),1);

for i=1:size(hmm_dist_matrix,1)
    for j=1:i
        qD(j,i) = hmm_dist_matrix(i,j);
        qD(i,j) = hmm_dist_matrix(i,j);
    end
end
[Y,eigvals] = cmdscale(qD);


% inspect whether 2D or 3D is enough:
format short g
[eigvals eigvals/max(abs(eigvals))]


subs_inds_text={};
for i=1:size(hmm_dist_matrix,1)
    subs_inds_text{end+1} = sprintf('S:%d',i);
end
% make the plot -- and also plot what were the answers:
figure;plot(Y(:,1),Y(:,2),'.')
for i=1:numel(subs_inds_text)
    text(Y(i,1)+0.01,Y(i,2),subs_inds_text{i});
end

title('2D Distances (HMM FO)');

% make the plot in 3D -- and also plot what were the answers:
% figure;plot3(Y(:,1),Y(:,2), Y(:, 3),'o')

% text(Y(:,1)+0.2,Y(:,2),Y(:, 3), answers)

%%
% from that figure, we can see the following:
Engaged = [1 11 15 16 3 18 4 10 14 8 9];

NotEngaged = [6 7 2 5 17 12 ];

% DidntLike = [12];

% Unknown = [13];


% now -- compare FO, for each state:
hmm_fo_Engaged = hmm_fo(Engaged,:);
hmm_fo_NotEngaged = hmm_fo(NotEngaged,:);
% hmm_fo_DidntLike = hmm_fo(DidntLike,:);


% make that into butterfly plot
figure('color','w');
for i=1:10
    subplot(1,10,i);
    
    plot(ones(size(hmm_fo_Engaged,1),1)-0.5, hmm_fo_Engaged(:,i),'ks');
    hold on;
    plot(ones(size(hmm_fo_NotEngaged,1),1)+0.5, hmm_fo_NotEngaged(:,i),'rs');
    
    % plot(ones(size(hmm_fo_DidntLike,1),1)+0.5, hmm_fo_DidntLike(:,i),'rx');

    xlim([0, 2]);
    ylim([0, 25]);
    
    xlabel(i);
    
    if i>1
        set(gca,'yticklabel',[],'ytick',[]);
    end
    set(gca,'xticklabel',[],'xtick',[]);
    
    [H,P,CI,STATS] = ttest2(hmm_fo_Engaged(:,i), [ hmm_fo_NotEngaged(:,i); hmm_fo_DidntLike(:,i)]);
    
    if H==1
        text(1,26,sprintf('P=%.2f',P),'horizontalalignment','center');
    end
    
end











%% Quick and DIRTY NBS:

addpath('../../hmm-mar-matlab/4_transition_probabilities/')
addpath('../hmm-necs-scripts/')
addpath(genpath('/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/4_transition_probabilities/NBSdirected/directed'));

Engaged = [11 15 16 3 18 4 10 14 8 9];
NotEngaged = [6 7 2 5 17  ];

mats_a = hmm_emps(:,:,Engaged);
mats_b = hmm_emps(:,:,NotEngaged);



C=cat(3,mats_a,mats_b);

nsubs=size(mats_a,3);

X=[ones(numel(Engaged),1) zeros(numel(Engaged),1); zeros(numel(NotEngaged),1) ones(numel(NotEngaged),1)];
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
 GLM.contrast=-1*[1 -1]; 
 %Type of test
 GLM.test='ftest'; % 'ttest' or 'ftest'
 %Exchange block for repeated measures
 %GLM.exchange=[]; STAT
 STATS.size='Extent'; %'Intensity' or 'Extent'
 %Threshold
 STATS.thresh=2.1; 
 %Significance (usually 0.05)
 STATS.alpha=0.15; 

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








