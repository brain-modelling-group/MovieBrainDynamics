% f7_cca

% getting the emp's and the emission prob

% clear all;close all;clc;

% load ../hmm-necs-scripts/outr.mat
cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript

addpath(genpath('../../hmm-mar-matlab/hmm_new/'))

% obtain/decode the data we need:


% make a 'sub list' - who did what
% only the subs that are present... oh
% man.

% re-run this twice - first WITH - second WITHOUT the 'bad guys'.

NSTATES=[8 10 24];





slist=struct();
slist.mov1a = [1:21];
slist.mov2a = [1:16 18:21];
slist.resta = [1:21];
slist.mov1b = [2:5 7:14 16:20];
slist.restb = [2:5 7:14 16:20];

% incorporate bad subjects -- that are OVERALL bad:
bads = [4 13 15]; bads = [];
fns=fieldnames(slist);
for i=1:numel(fns)
    for j=1:numel(bads)
        slist.(fns{i})(slist.(fns{i})==bads(j))=[];
    end
end

% then incorporate the occasioal bads:
slist.mov2a(slist.mov2a==11)=[];
slist.mov2a=sort([slist.mov2a 15]);
slist.resta(slist.resta==5)=[];
slist.resta(slist.restb==6)=[];


NULL_FIRST = 5;

J=14;  % the number of networks!


analyses = {...
    {'mov1a'},'mov1a---2'; ... % only mov1
    {'mov2a'},'mov2a---2'; ... % only mov2 (repeat on d1 from 11min onwards)
    {'mov1b'},'mov1b---2'; ... % only mov2 day2
    {'mov1a','mov1b'},'mov1a-1b---2'; ... % mov1 day and and day 2
    {'mov1a','mov2a'},'mov1a-2a---2'; ... % mov1 day1 + repeat on day 1
    {'resta'},'resta---2'; ...  % only resta
    {'restb'},'restb---2'; ...  % only rest day 2
    {'resta','restb'},'resta-b---2'; ...  % rest a and rest day 2
    {'mov1a','mov1b','resta','restb'},'all---2'; ...  % all movie and rest day 1 aand day 2 but not repeats
    {'mov1a','mov2a','mov1b','resta','restb'},'allPlus---2'; ...  % ALL functional data.
    {'mov1a','resta'},'sessiona---2'; ...
    };



% run hmm for normal, aroma, and aroma-gsr extracted timeseries.



for i_NSTATES=2 %1:numel(NSTATES)
    this_NSTATES=NSTATES(i_NSTATES);
    
    covregtypes={'normal','aroma','aroma-gsr'};
    covregpres={'','aroma-','aroma-gsr-'};
    
    for i_covregtype=2 % 1:3
        covregtype = covregtypes{i_covregtype};
        covregpre  = covregpres{i_covregtype};
        
        for i_analysis = 1 %6 % 1:size(analyses,1)
            
            
            scan_names = analyses{i_analysis, 1};
            %             save_dir =  ['results_new2_' num2str(this_NSTATES) '/' covregtype '/' analyses{i_analysis, 2}];
            %
            %             if exist(save_dir,'dir')
            %                 rmdir(save_dir,'s');
            %                 mkdir(save_dir);
            %             else
            %                 mkdir(save_dir);
            %             end
            
            
            % scan_names = {'mov1a','mov1b'};
            % scan_names= {'mov1a'};
            
            subs_to_use = [1:21]; % magic to make sure that IF you specified more than one functional task,
            % subjects that do NOT have done them ALL are excluded.Otherwise
            % hmm gets confused, I think.
            
            nsubs=21;
            for i=1:numel(scan_names)
                subs_to_use = intersect(subs_to_use, slist.(scan_names{i}));
                nsubs=numel(subs_to_use);
            end
            
            alld={};
            T=[];% whos%
            
            for i_scan_name = 1:numel(scan_names)
                scan = scan_names{i_scan_name};
                % scan='mov1';
                % d=dir(['extracted_timeseries/*' scan '*r01.txt']);
                d=[];
                for i=1:nsubs
                    d = [d dir(['../extracted_timeseries/' covregpre 'ts-' scan '*' sprintf('s%.2d',subs_to_use(i)) '*r01.txt'])];
                end
                
                
                switch scan
                    case {'mov1a','mov1b'}
                        nscans=535;
                    case {'mov2a'}
                        nscans=250;
                    case {'resta','restb'}
                        nscans=220;
                end
                
                
                scans_sum=0;
                % alld={};
                
                for i=1:numel(d);
                    
                    
                    scans_sum = scans_sum + nscans - NULL_FIRST;
                    % subj = regexp(fname,'s[0-9]{2}','match'); subj=subj{1};
                    
                    alld{end+1} = [];
                    for j=1:14
                        
                        
                        fname=regexprep(['../extracted_timeseries/' d(i).name],'r01.txt',sprintf('r%.2d.txt',j));
                        
                        
                        disp(fname)
                        v=load(fname);
                        
                        % SOOOOoooo.... subject s03, restb, is ONLY 219
                        % volumes.
                        tmp_sub=str2num(regexprep(fname,'.*-s([0-9]{2}).*','$1'));
                        if tmp_sub == 3 && strcmp(scan,'restb')
                            v(end+1)=v(end);
                        end
                        
                        disp(length(v));
                        % discard first 5 volumes (!!!)
                        v(1:NULL_FIRST) = [];
                        
                        
                        dv=detrend(v,'constant');
                        
                        % NOT smooth yet...
                        % fdv=filtfilt([1,1,1]/3,1,dv);
                        
                        sfdv=dv/std(dv);
                        
                        
                        alld{end} = [alld{end} sfdv];
                    end
                    
                    
                end
                
                T(end+1) = scans_sum;
                
                % so - count how many??
                
                
                
            end
            
            
            % do another iteration later - to figure out IF this makes a
            % huge impact...
            % T=sum(T);  % concat to one big T.., no differentiation between sessions?
            
            
            
            % merge -- this is the data going into hmm-mar
            dat=[];
            for i=1:numel(alld)
                dat=[dat; alld{i}];
            end
            
            
            
            % smooth here:
            ndat=dat;
            for i=1:14
                ndat(:,i) = smoothdata(dat(:,i));
            end
            dat=ndat;
            
            % [T, J] = size(dat);
            
            
            % run the hmm.
            %             % addpath(genpath('hmm_new'))
            %
            %
            %             tic;
            %             % addpath(genpath('HMM-MAR-master/'));
            %
            %             % load your data -- single subject? - there is an option to specify segment
            %             % lengths -- where is that? Are all segments equal length?
            %             % dat = rand(100,3); % random data jst as an example
            %             % n_time points X n_networks format
            %
            %             % [T,J] = size(dat); % T = number of time points in concatenated data
            %             % J = number of networks if method is same as ours
            %
            %             n_sub = numel(subs_to_use); % no. of subjects(used for concatenation)
            %             K = this_NSTATES; % no. states
            %             reps = 15; % to run it multiple times (saves all the results
            %             % as a seperate mat file)
            %
            %             TR = 2.20;  % TR of rest data
            %
            %             % options for model estimation; I have set the usual choices here;
            %             % no need to change anything in geenral
            %             options = struct();
            %             options.K = K; % number of states
            %             options.order = 0; % no autoregressive components
            %             options.zeromean = 0; % model the mean
            %             options.covtype = 'full'; % full covariance matrix
            %             options.Fs = 1/TR;
            %             options.verbose = 1;
            %             options.standardise = 1;
            %             options.inittype = 'HMM-MAR';
            %             options.cyc = 500;
            %             options.initcyc = 10;
            %             options.initrep = 3;
            %
            %             % run the HMM multiple times
            %             for r = 1:reps
            %                 disp(['RUN ' num2str(r)]);
            %                 [hmm, Gamma, ~, vpath, ~, ~, fe] = hmmmar(dat,T,options);
            %                 save([save_dir '/HMMrun_rep_' num2str(r) '.mat'],'Gamma','vpath',...
            %                     'hmm','T','J','K','n_sub','fe');
            %             end
            %
            %
            %             % do the summary measures:
            %
            %             %     % clear;clc; tic;
            %             %     reps = 5;
            %             %     do = 1; % do = Flag to replace 0s with NaNs (in all summary measures)
            %             %         % do=0 keeps all misisng values as zeros; do=1 sets them
            %             %         % to Nans; this is important while comparing between groups!!
            %             %
            %             %     for r = 1:reps
            %             %         load([save_dir '/HMMrun_rep_' num2str(r) '.mat'],'Gamma','vpath',...
            % %             %               'hmm','J','T','K','n_sub');
            %             %
            %             %         mean_em=zeros(J,K); % mean emissions
            %             %         for k=1:K
            %             %             mean_em(:,k)=getMean(hmm,k);
            %             %         end
            %             %
            %             %         prob=hmm.P; % transition probabilities
            %             %
            %             %         [FO,dign,avg_life]=summary_measures(vpath,n_sub,do);
            %             %
            %             %         save([save_dir '/Summary_measures_rep_',num2str(r) '.mat']);
            %             %     end
            
            
        end
        
        
        
    end
    
    
end



% so that gives us dat -- now load in the hmm from our favorite HMM run:
% also -- assuming NO bads:, those are: [4 13 15] -- we take care of that in our CCA
dat=dat;
vars = load('../results_new2_10/aroma/all---2/HMMrun_rep_13.mat');
hmm=vars.hmm;
%%

[Gamma,Xi] = hmmdecode(dat,T,vars.hmm,0);
[vpath] = hmmdecode(dat,T,vars.hmm,1);


%%
% let's plot it:
% (see f7_cca_sup1)


%%




% so we could - theoretically - check things out...

% use our hmm to decode everything in mov1a - - then in second iteration
% also the discareded guys from mov1a.

% we then decide how many the states & and question answers - we are going
% to be using.

% then we check EACH PERMUTATION of states and questions and repeat the
% CCA, to see if we can get any significant result. Hopefully things will
% make sense. Otherwise, we are doomed.


% once that is done, plot the loadings/etc like discussed earlier.




% make the CCA, matrix A - 1 per subj.
% we will use all subs in this ST, right..

% hmm=[];
% part='mov1a';
% analysis='all';
% preprocessing='aroma';



% we'd need to re-calculate the emp and st from vpath:
% emps = outr.run13.(preprocessing).(analysis).(part).emp(:,:,:);
% st=outr.run13.(preprocessing).(analysis).(part).vpath';
% emps = emp-per-subject, code is somewhere used previously to generate
% f6_transitions...

% so we'd need to re-calculate smps and st, from code that we obtained
% somewhere below:

selected_subs=1:21;
nvol=530;
st_path=vpath(1:530*21);

NSTATES=10;
S=NSTATES; % max(st_path); % Number of states -- should be 12!

% create the subject-specific bpaths:
sub_specific_vpath=[];

% ... do all the magick
subj_emp = zeros(S,S,numel(selected_subs));
emps=subj_emp;




% n_time=length(st_path)/n_sub;
n_sub = numel(selected_subs);
n_time = nvol;

%%

for sub=1:n_sub
    emp=zeros(S,S);
    stat=squeeze(st_path((sub-1)*n_time+1:sub*n_time));
    sub_specific_vpath(:,end+1) = stat;
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
    row_has_all_zeros = ~any(emp, 2);
    indices = find(row_has_all_zeros);
    for k=1:length(indices)
        emp(indices(k),indices(k))=1;
    end
    subj_emp(:,:,sub)=emp;
end


emps=subj_emp;

% st=reshape(st_path,n_sub,[]); % reshape to n_sub x n_time format
st=reshape(vpath,[],n_sub)';
S=10; % Nstates = 10.
n_sub=size(st,1);
n_time=size(st,2);


fo=zeros(n_sub,S);

for kk=1:n_sub
    stat=st(kk,:); % getting state path seq of individual subject
    for k=1:S
        fo(kk,k)=sum(stat==k); % counting numbe rof time points spent in each state
    end
end
%fo=100*fo./176%(sum(fo,2));% normailzing to get the % time spent in each state
fo=100*fo./n_time;


%% SO we got all the data we need - now it's a matter of figuring out what to do with it...
SELSUBS = [1:3 5:12 14 16:21]; % since 4, 13 and 15 are bad... -- this for mov1a

% SELSUBS(end-1) = [];

% we got the questions...
m=load('../../questionnaires/m.mat');
m=m.m;
% m=m(m(:,20)==1,:); -- do NOT select, for now...

questions=m(:,1:4);
questions=questions(SELSUBS,:); % qtitles = 'bored','enjoyed','emitional','audio' = 1:4; 6:9 is session b


% get the (full) transition matrices:
hmm_emps = emps(:,:,SELSUBS);

hmm_fo = fo(SELSUBS,:);

OMIT_ZEROS_FROM_TRANSITION_MATRICES = 1;
OMIT_TRACE_FROM_TRANSITION_MATRICES = 1;
USE_ONLY_TRACE_FROM_TRANSITION_MATRICES = 0;
SKIP_ZERO_IN_PTEST = 0;
SAVEIT = 1;
QUESTION_INDICES = [1:4];

HMM_FEATURE = 'TRANSITIONS'; % or: 'TRANSITIONS';
% HMM_FEATURE = 'TRANSITIONS';

% we only look at STATE TRANSITIONS here
% we can look also at STATE OCCUPANCY (using fractional occupancy)
% or we can look at STATE ROBUSTNESS (i.e. chance of staying in a state)

hmm_dist_matrix = [];
questions_dist_matrix = [];
hmm_dist_vec = [];
questions_dist_vec = [];

%  normalized squared euclidean distance
sq_euclid = @(u, v) 0.5 * var(u-v) / (var(u) + var(v)); % use the L2 Norm.
sq_euclid_option2 = @(u, v) pdist([reshape(u,1,[]); reshape(v,1,[])],'euclid');
correlation_option = @(u, v) pdist([reshape(u,1,[]); reshape(v,1,[])],'correlation');


get_upper_inds = @(x) find(triu(reshape(1:(x^2),x,[]),1)>0);

get_lower_inds = @(x) find(tril(reshape(1:(x^2),x,[]),-1)>0);

get_tr_inds = @(x) diag(reshape(1:(x^2),x,[]));

ORIG_SELSUBS = SELSUBS;
NPERMS = 5000;

s=struct(); % for saving our permutation results...

for iPerm=1:NPERMS
    
    s(iPerm).PermutedSubs = randperm(length(ORIG_SELSUBS));
    
    if iPerm == 1  % make sure that our first permutation is the original one...
        s(iPerm).PermutedSubs = 1:length(ORIG_SELSUBS);
    end
    
    % shuffle the questions according to that:
    permutedQuestions = questions(s(iPerm).PermutedSubs,:);
    
    questions = questions(:,QUESTION_INDICES);
    for i=1:numel(SELSUBS)
        for j=1:i
            if j ~= i
                
                hmm_emps_i = hmm_emps(:,:,i);hmm_emps_i=hmm_emps_i(:);
                hmm_emps_j = hmm_emps(:,:,j);hmm_emps_j=hmm_emps_j(:);
                
                % find region where it's > 0 (otherwise maybe not fair
                % calculating distance?
                sel_vals_i = find(hmm_emps_i>0);
                sel_vals_j = find(hmm_emps_j>0);
                
                if OMIT_ZEROS_FROM_TRANSITION_MATRICES == 1
                    select_these_vals = intersect(sel_vals_i, sel_vals_j);
                else
                    select_these_vals = 1:100;
                end
                
                if OMIT_TRACE_FROM_TRANSITION_MATRICES == 1
                    trace_indices = [1    12    23    34    45    56    67    78    89   100];
                    select_these_vals=setxor(select_these_vals, trace_indices);
                end
                
                if USE_ONLY_TRACE_FROM_TRANSITION_MATRICES == 1  % the other 2 options must be 0
                    trace_indices = [1    12    23    34    45    56    67    78    89   100];
                    select_these_vals = trace_indices;
                    
                end
                
                hmm_emps_i=hmm_emps_i(select_these_vals);
                hmm_emps_j=hmm_emps_j(select_these_vals);
                
                
                
                
                
                hmm_fo_i = hmm_fo(i,:);
                hmm_fo_j = hmm_fo(j,:);
                
                questions_i = permutedQuestions(i,:);
                questions_j = permutedQuestions(j,:);
                
                
                
                
                % keyboard;
                %
                
                if strcmp(HMM_FEATURE,'FO')
                    this_u = hmm_fo_i;
                    this_v = hmm_fo_j;
                elseif strcmp(HMM_FEATURE,'TRANSITIONS')
                    this_u = hmm_emps_i;
                    this_v = hmm_emps_j;
                end
                
                
                hmm_dist_matrix(i, j) = correlation_option(this_u, this_v);
                hmm_dist_vec(end+1) = hmm_dist_matrix(i, j);
                % if hmm_dist_matrix(i, j) == 0
                %     keyboard;
                % end
                questions_dist_matrix(i, j) = sq_euclid_option2(questions_i, questions_j);
                questions_dist_vec(end+1) = questions_dist_matrix(i, j);
                
            end
        end
    end
    
    
    
    if SKIP_ZERO_IN_PTEST == 1
        
        inds_hmm = find(hmm_dist_vec>0);
        inds_questions = find(questions_dist_vec>0);
        
        pval_inds = intersect(inds_hmm, inds_questions);
    else
        pval_inds = 1:length(hmm_dist_vec);
    end
    
    
    
    % figure out whether there's a correlation:
    matrix_to_test = [hmm_dist_vec(pval_inds)', questions_dist_vec(pval_inds)'];
    [RHO,PVAL] = corr(matrix_to_test);
    
    if iPerm == 1
        s(iPerm).hmm_dist_matrix = hmm_dist_matrix;
        s(iPerm).hmm_dist_vec = hmm_dist_vec;
        s(iPerm).questions_dist_matrix = questions_dist_matrix;
        s(iPerm).questions_dist_vec = questions_dist_vec;
        s(iPerm).pval_inds = pval_inds;
        s(iPerm).matrix_to_test = matrix_to_test;
    end
    s(iPerm).RHO = RHO;
    s(iPerm).PVAL = PVAL;
    
    
end



fh=figure('color','w');
set(fh,'position',[ 265         200        1094         600]);
subplot(2,2,1);
imagesc(s(1).hmm_dist_matrix);
colormap hot
colorbar;
title('hmm state transition distances');
xlabel('participant');
ylabel('participant');

subplot(2,2,2);
imagesc(s(1).questions_dist_matrix);
colormap hot
colorbar;
title('interview distances');
xlabel('participant');
ylabel('participant');


subplot(2,2,3)
plot(s(1).matrix_to_test(:,1), s(1).matrix_to_test(:,2),'bO');
xlabel('hmm distances');
ylabel('interview distances');
title(sprintf('rho = %.3f, pval = %.3f',s(1).RHO(1,2),s(1).PVAL(1,2)));
% ylim([0, 6]);
H=lsline;
set(H,'color','r');


subplot(2,2,4)
% plot histogram of ALL rho values - from 2 onwards.
all_rho = [];
for i=1:numel(s)
    all_rho(end+1) = s(i).RHO(1,2);
end
histogram(all_rho(2:end));
% what are the lower and upper quantile 
l_05 = quantile(all_rho,0.05);
l_95 = quantile(all_rho,0.95);
line(l_05*[1 1],get(gca, 'ylim'),'color',0.4*[1 1 1]);
line(l_95*[1 1],get(gca, 'ylim'),'color',0.4*[1 1 1]);
line(s(1).RHO(1,2)*[1 1], get(gca,'ylim'),'color','r');
xlabel('rho');
ylabel('permutation count');

% cle


if SAVEIT == 1
    set(fh,'paperunits','centimeters');
    set(fh,'papersize',[25 20],'paperposition',[0 0 20 15]);
    
    ftarget_file = 'f8_brainhmmdist_sup6.jpg';
    
    ftarget1=['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' ftarget_file];
    ftarget2=['/home/johan/ldrive/Lab_LucaC/10_Johan/1_HMMMovie/' ftarget_file];
    
    saveas(fh,ftarget1);
    saveas(fh,ftarget2);
    
    saveas(fh,[ftarget1(1:end-3) 'pdf']);
    
    
    cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript
end




%%
%
%
%
%
% so here's the logic to determine which features to select, etc.
%
%
%fo
%
%
%
% % first we go through the STATES, 3 at a time:
%
%
% NVARS = 2;
% NPROBS = 2;
%
% % hm, rather arbitrary choice of numbers. But, well.
% SELSUBS = [1:3 5:12 14 16:21]; % bads = [4 13 15];
%
%
% % list of all the transitions:
% all_trans = [nchoosek(1:10,2);nchoosek(10:-1:1,2)];
%
%
% fo = [1, 2, 3, 4, 6, 7];
% TR_PROBS = [[6, 4]; [3, 2]];
%
% % EMP_STATES = [4, 6];
% % TR_PROBS = [[6, 4]];
% C = nchoosek(1:10,NVARS);
%
% list_two_combos=nchoosek(1:90,NPROBS);
%
%
% for Li=1:size(list_two_combos,1)
%
%     TR_PROBS = all_trans(list_two_combos(Li,:),:);
%
%
%     for Ci=1:size(C,1)
%
%
%         EMP_STATES = C(Ci,:);
%         % TR_PROBS = [];
%
%
%         m=[];
%         for i=1:n_sub
%
%             s_emp =  emps(:,:,i);
%             s_vp = st(i,:);
%
%             for ii=1:numel(EMP_STATES)
%                 m(i,ii) = sum(s_vp==EMP_STATES(ii))/numel(s_vp);
%             end
%
%             for ii=1:size(TR_PROBS,1)
%                 m(i,numel(EMP_STATES)+ii) = s_emp(TR_PROBS(ii,1), TR_PROBS(ii,2));
%             end
%
%         end
%         hmm=m;
%         % select SUBS here
%         hmm=hmm(SELSUBS,:);
%
%
%         % figure;imagesc(hmm);colorbar;
%
%
%
%
%
%         % how - get the answers to questionnaires.
%         % for overview what M is composed of -- see the XLS sheet.
%         % m columns 1:4 are: bored, enjoy, emotional, audioQ, mov1a.
%         % m columns 6:9 are same as for mov1b
%         % column 20 is whether we're using the scroes or not for 'all' analysis.
%         m=load('../../questionnaires/m.mat');
%         m=m.m;
%         % m=m(m(:,20)==1,:); -- do NOT select, for now...
%
%         questions=m(:,1:3);
%         questions=questions(SELSUBS,:);
%         % and here..
%
%
%
%         % figure;imagesc(questions);colorbar;
%
%
%
%         % use this in order to find our TRUE result.
%         [A,B,r,U,V,stats] = canoncorr(hmm,questions);
%
%
%         % let's permute - around 5000 times. Store all the Permuted Results into
%         % big matrices.
%         NUMPERMS=5000;
%         p=struct();
%         for perm=1:NUMPERMS
%
%             perm_hmm = hmm(randperm(size(hmm,1)),:);
%
%             [pA,pB,pr,pU,pV,pstats] = canoncorr(perm_hmm,questions);
%
%             p(perm).perm_hmm = perm_hmm;
%             p(perm).A = pA;
%             p(perm).B = pB;
%             p(perm).r = pr;
%             p(perm).U = pU;
%             p(perm).V = pV;
%             p(perm).stats = pstats;
%
%
%             if pr < 0
%                 keyboard;
%             end
%
%         end
%
%         % allstats = 0;
%         statsWilks_1 = [];
%         for perm=1:NUMPERMS
%             statsWilks_1(end+1) = p(perm).stats.Wilks(1);
%         end
%
%         corrs_1 = [];
%         for perm=1:NUMPERMS
%             corrs_1(end+1) = p(perm).r(1);
%         end
%
%         fh=figure;
%
%
%         ah=subplot(2,2,1);
%
%         hist(corrs_1,100)
%
%
%
%         ts= sprintf('%d, ',EMP_STATES);
%         for TRi=1:size(TR_PROBS,1)
%             ts=[ts sprintf(' %d -> %d  ',TR_PROBS(TRi,1),TR_PROBS(TRi,2))];
%         end
%         % ts = [ts sprintf('%.3f', r(1))];
%         % ts = [ts sprintf('\nWilks=%.3f ', stats.Wilks(1))];
%         ts = [ts sprintf('\np = %.3f ', stats.p(1))];
%         title(ts);
%         line(r(1)*[1 1], get(gca,'ylim'),'color','r');
%         xlabel('correlation value');
%         ylabel('count (out of 5000)');
%         q_p05 = quantile(corrs_1,0.05);
%         q_p95 = quantile(corrs_1,0.95);
%
%         q_p01 = quantile(corrs_1,0.01);
%         q_p99 = quantile(corrs_1,0.99);
%
%
%         line(q_p05*[1 1], get(gca,'ylim'),'color','k');
%         line(q_p95*[1 1], get(gca,'ylim'),'color','k');
%
%         line(q_p01*[1 1], get(gca,'ylim'),'color','k');
%         line(q_p99*[1 1], get(gca,'ylim'),'color','k');
%
%
%         % next up -- loadings, a la Zalesky style - just take the first
%         % canonical variate. We need to plot A(:,1) and B(:,1)
%         % and calculate hmm correlated with all U, and quesiton correlated with
%         % V - check diagonals.
%         a=subplot(2,2,3);
%         bar(A(:,1));
%         corrvals_hmm = corr(hmm,U(:,1));
%         for icorrvals=1:numel(corrvals_hmm)
%             th=text(icorrvals,A(icorrvals,1),sprintf('%.2f',corrvals_hmm(icorrvals)));
%             set(th,'verticalalignment','middle','horizontalalignment','center');
%         end
%         box off
%         title('state coeff & loadings');
%
%         a=subplot(2,2,4);
%         bar(B(:,1));
%         corrvals_questions = corr(questions,V(:,1));
%         for icorrvals=1:numel(corrvals_questions)
%             th=text(icorrvals,B(icorrvals,1),sprintf('%.2f',corrvals_questions(icorrvals)));
%             set(th,'verticalalignment','middle','horizontalalignment','center');
%         end
%
%         box off
%         title('questions coeff & loadings');
%
%
%         a=subplot(2,2,2);
%         plot(U(:,1),V(:,1),'bO');
%         xlabel('states');
%         ylabel('answers');
%         ts=[];
%         ts = [ts sprintf('r = %.3f', r(1))];
%         title(ts);
%         % fit a straight line, too, as best as you can, pls:
%         % P=polyfit(U(:,1),V(:,1),1);
%         % this_xl=get(gca,'xlim');
%         % line(this_xl,P(1)+P(2)*this_xl,'color','r');
%         tr=refline;
%         set(tr,'color','r');
%
%
%         set(fh,'paperunits','centimeters');
%         set(fh,'papersize',[15 10],'paperposition',[0 0 15 10]);
%
%         ftarget_file = sprintf('f7_cca_sup2_%dcombo%d%d.jpg', NVARS, Ci,Li);
%         ftarget1=['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' ftarget_file];
%         ftarget2=['/home/johan/ldrive/Lab_LucaC/10_Johan/1_HMMMovie/' ftarget_file];
%
%         saveas(fh, ftarget1);
%         saveas(fh, ftarget2);
%
%         %     if r(1) > q_p95 || r(1) < q_p05
%         %         ftarget_file = sprintf('f7_cca_sup2_%dsig_combo%d.jpg', NVARS, Ci);
%         %         ftarget1=['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' ftarget_file];
%         %         ftarget2=['/home/johan/ldrive/Lab_LucaC/10_Johan/1_HMMMovie/' ftarget_file];
%         %         saveas(fh, ftarget1);
%         %         saveas(fh, ftarget2);
%         %     end
%         %
%         %       if r(1) > q_p99 || r(1) < q_p01
%         %         ftarget_file = sprintf('f7_cca_sup2_%d99sig_combo%d.jpg', NVARS, Ci);
%         %         ftarget1=['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' ftarget_file];
%         %         ftarget2=['/home/johan/ldrive/Lab_LucaC/10_Johan/1_HMMMovie/' ftarget_file];
%         %         saveas(fh, ftarget1);
%         %         saveas(fh, ftarget2);
%         %     end
%
%
%         if stats.p(1) < 0.05
%             ftarget_file = sprintf('f7_cca_sup2_%dp05sig_combo%d%d.jpg', NVARS, Ci,Li);
%             ftarget1=['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' ftarget_file];
%             ftarget2=['/home/johan/ldrive/Lab_LucaC/10_Johan/1_HMMMovie/' ftarget_file];
%             saveas(fh, ftarget1);
%             saveas(fh, ftarget2);
%         end
%
%
%         close(fh);
%
%         cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript
%
%
%     end
%
% end



% some code to save...

%
% set(newbigf,'paperunits','centimeters');
% set(newbigf,'papersize',[25 20],'paperposition',[0 0 25 20]);
%
% fsource = 'f6_transitions_sup1.jpg';
% ftarget_file = fsource;
% saveas(newbigf,fsource);
%
%
% ftarget1=['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' ftarget_file];
% ftarget2=['/home/johan/ldrive/Lab_LucaC/10_Johan/1_HMMMovie/' ftarget_file];
%
% copyfile(fsource,ftarget1);
% copyfile(fsource,ftarget2);
%
% delete(fsource);
%
%
% cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript
%
%
%
%
