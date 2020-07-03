% so we alter the following function so as to work with OUR data.
% the emp_prob should work with our slists
% we also already know the nvol(s)
% so instead of subj_emp, we should get like a struct (maybe) with
% different subj_emps -- that depend on analysis, preprocessing strategy,
% and RUN.

% this function recovers from all of the analyses and options the subject's
% transition matrix.

% I update it to also get the subject's specific (and session specific)
% viterbi path.





clear all;close all;clc;

ANALYSIS='all';
NSTATES=10;

this_cwd = pwd;
load(['valid_inferences_' ANALYSIS '.mat']);



% cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/hmm-necs-scripts/

outr=struct();


analyses = {...
    {'mov1a'},'mov1a'; ... % only mov1
    {'mov1b'},'mov1b'; ... % only mov2 day2
    {'mov1a','mov1b'},'mov1a-1b'; ... % mov1 day and and day 2
    {'resta'},'resta'; ...  % only resta
    {'restb'},'restb'; ...  % only rest day 2
    {'resta','restb'},'resta-b'; ...  % rest a and rest day 2
    {'mov1a','mov1b','resta','restb'},'all'; ...  % all movie and rest day 1 aand day 2 but not repeats
    };


% deal with the subjects...
slist=struct();
slist.mov1a = [1:21];
slist.mov2a = [1:16 18:21];
slist.resta = [1:21];
slist.mov1b = [2:5 7:14 16:20];
slist.restb = [2:5 7:14 16:20];

% incorporate bad subjects -- that are OVERALL bad:
bads = [4 13 15];
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



preprocessings={'aroma'};


for i_summarymeasures=valid_inferences
    SUMMARY_MEASURES_FILE=i_summarymeasures;
    
    for i_preprocessing=1:numel(preprocessings)
        preprocessing_dir=preprocessings{i_preprocessing};
        
        
        i_analysis = find(strcmp(analyses(:, 2), ANALYSIS));
        
        analysis=regexprep(analyses{i_analysis,2},'-','_');
        scan_names = analyses{i_analysis, 1};
        
        subs_to_use = [1:21]; % magic to make sure that IF you specified more than one functional task,
        % subjects that do NOT have done them ALL are excluded.Otherwise
        % hmm gets confused, I think.
        
        nsubs=21;
        for i=1:numel(scan_names)
            subs_to_use = intersect(subs_to_use, slist.(scan_names{i}));
            nsubs=numel(subs_to_use);
        end
        
        selected_subs=subs_to_use;
        disp(['analysis = ' analysis]);
        outr.(sprintf('run%d',i_summarymeasures)).(regexprep(preprocessing_dir,'-','_')).(analysis).subs=selected_subs;
        nsubs=numel(selected_subs);
        
        
        
        load([this_cwd '/../results_10/' preprocessing_dir '/' analyses{i_analysis,2} '/HMMrun_rep_' num2str(SUMMARY_MEASURES_FILE) '.mat']); % this will get our vpath.
        
        % so now we have the 'vpath'. Figure out below what the ONSET
        % is where we should start looking?
        
        my_onset=0;
        for j_subpart=1:numel(analyses{i_analysis,1})
            
            
            analysis=regexprep(analyses{i_analysis,2},'-','_');
            % now we have vpath.
            this_subpart=analyses{i_analysis,1}{j_subpart};
            
            % IF - in the future - we're going to include >1 subject..
            % and figure out how to do that with HMM-MAR...
            outr.(sprintf('run%d',i_summarymeasures)).(regexprep(preprocessing_dir,'-','_')).(analysis).(this_subpart).nsub = selected_subs;
            
            disp(['preproc_dir = ' regexprep(preprocessing_dir,'-','_')]);
            disp(['analysis = ' analysis]);
            disp(['this_subpart = ' this_subpart]);
            
            switch this_subpart
                case {'mov1a','mov1b'}
                    nvol=530; % it IS actually 215 now; not 530!
                case {'resta','restb'}
                    nvol=215;
                    
            end
            
            
            
            % st_path is what we select it to be, now:
            st_path=vpath(my_onset+(1:(nvol*numel(selected_subs))));
            
            NSTATES=10;
            S=NSTATES; % max(st_path); % Number of states -- should be 12!
            
            
            % create the subject-specific bpaths:
            sub_specific_vpath=[];
            
            
            % ... do all the magick
            subj_emp = zeros(S,S,numel(selected_subs));
            
            
            
            
            % n_time=length(st_path)/n_sub;
            n_sub = numel(selected_subs);
            n_time = nvol;
            
            
            
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
            
            %emp+
            
            
            outr.(sprintf('run%d',i_summarymeasures)).(regexprep(preprocessing_dir,'-','_')).(analysis).(this_subpart).vpath=sub_specific_vpath;
            outr.(sprintf('run%d',i_summarymeasures)).(regexprep(preprocessing_dir,'-','_')).(analysis).(this_subpart).emp=subj_emp;
            
            
            
            % add to onset (i.e. to be able to move to the next task?
            my_onset = my_onset + nvol * numel(selected_subs);
            
            
            
            
        end
    end
end



save outr.mat outr
%
% function [fo] = fract_occ(st_path,n_sub)
% %This function calculates the fractional occupancy in different states for
% %each subject
% % Inputs:
% %         st_path = Viterbi-decoded state path sequence (concatenated), as
% %                   outputted from fit_hmm. Size =1 X n_sub*n_time
% %           n_sub = Number of subjects (as per the data used in the process
% %                   of concatenation. Size 1 X 1
% % Ouptputs     fo = % Fractional Occupancy for each subject. Size n_sub X S
% %                   where S is the number of hidden states
% %
% S=max(st_path);
% st=reshape(st_path,n_sub,[]); % reshape to n_sub x n_time format
% fo=zeros(n_sub,S);
%
% for kk=1:n_sub
%     stat=st(kk,:); % getting state path seq of individual subject
%     for k=1:S
%         fo(kk,k)=sum(stat==k); % counting numbe rof time points spent in each state
%     end
% end
% %fo=100*fo./176%(sum(fo,2));% normailzing to get the % time spent in each state
% n_time=length(st_path)/n_sub %Changed by Luca
% fo=100*fo./n_time
% end

%
%
%
%
%
%
%                 function [ subj_emp ] = emp_prob( st_path,n_sub )
%                 %this function calculates empoirical probability for every subject
%                 % Inputs:
%                 %         st_path = 1 x n_sub*n_time vector as outputted from fit_hmm
%                 %         function
%                 %           n_sub = Number of subjects
%                 % Output:
%                 %        subj_emp = Empirical transition probabilitites calculated for
%                 %        every subject. Size = SxSxn_sub
%

