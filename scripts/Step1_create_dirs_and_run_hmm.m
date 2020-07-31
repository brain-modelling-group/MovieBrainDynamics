
% this script does the following
% (1) set up the results folder structure
%
% results
% ├── aroma
% │   ├── all
% │   ├── mov1a
% │   ├── mov1b
% │   ├── mov1a-1b
% │   ├── resta
% │   └── restb
% │   ├── resta-b
% ├── normal
% │   ├── all
% │   ├── mov1a
% │   ├── mov1b
% │   ├── mov1a-1b
% │   ├── resta
% │   └── restb
% │   ├── resta-b
% ├── aroma-gsr
% │   ├── all
% │   ├── mov1a
% │   ├── mov1b
% │   ├── mov1a-1b
% │   ├── resta
% │   └── restb
% │   ├── resta-b
%
% aroma = use fmriprep's automated aroma preprocessing (as outlined in the
%         paper) as the to deal with motion. After AROMA, covariate
%         regression step only includes regressor from the white matter,
%         regressor from the CSF ('CSF and 'WhiteMatter' from the
%         fMRIPREP's .csv outout), and the boxcar regressors from applying
%         the threshold of 0.4 on to the 'FramewiseDisplacement' .csv
%         output of fMRIPREP.
% 
% normal = use the 6 motion parameters ('X', 'Y', 'Z', 'RotX', 'RotY',
%          'RotZ') from the fMRIPREP's output to deal with motion. In
%          addition, in the single covariate regression, also includes
%          regressors from White Matter and CSF ('CSF and 'WhiteMatter'
%          from the fMRIPREP .csv output), as well as the boxcar regressors
%          as described above.
%
% aroma-hsr = exactly as with aroma, only in the covariate regression, the
%             Global Signal Regression is also taken into account. The
%             'GlobalSignal' regressor is included from the fMRIPREP's .csv
%             output.
%
% [normal and aroma-hsr are available upon request]
%
% in all cases, covariate regression was performed with the 'clean_img'
% function of the image module within the python nilearn package
% (https://nilearn.github.io/modules/generated/nilearn.image.clean_img.html)
%
%
% (2) read in the time-series data. The time-series data have been exported
% as plain text files. Each of the 14 Brain Networks has one text file
% encoding the average BOLD signal from that network. The time series have
% been obtained by simple multiplication of a binary mask of the network to
% the fMRI images after they have been preprocessed with fMRIPREP and after
% covariate regression.
%
%
% (3) Concatenate the time-series data of all participants and all sessions
% into a matrix of 14 (networks) by N (total timepoints) as an input to the
% HMM Analysisis. The total timepoints depends upon the selected scans put
% into the analysis.
%
%
% (4) Run the HMM Analysis and save the results into a .mat file.
%
%
% Different analysis are run for:
% - Number of brain states (8, 10, 12, 24)
% - For each brain state, 15 iterations are run
% - This is repeated for different conditions
%


% matlab might complain about pca; since matlab does not have any modules
% structure, everything is in the same path.
% rmpath(fileparts(which('pca')))

%
%
% this is for looping over total amount of brain nstates
%
%
NSTATES=[10]; % we run only 10 states...
HMMREPS = 15; % re-run the HMM analysis 15 times.



% run the hmm.
addpath(genpath('hmm_new')) % load in the hmm toolbox + additional functions




analyses = {...
    {'mov1a'},'mov1a'; ... % only day 1 movie
    {'mov1b'},'mov1b'; ... % only day 2 movie
    {'mov1a','mov1b'},'mov1a-1b'; ... 
    {'resta'},'resta'; ...  % only dat 1 rest
    {'restb'},'restb'; ...  % only day 2 rest
    {'resta','restb'},'resta-b'; ...  
    {'mov1a','mov1b','resta','restb'},'all'; ...  % all
    };

WHICH_ANALYSES = [1 7]; % 

TRUNCATE_VOLS_TO_220 = 0; % matlab does not have True and False, so 0 and 1 is the way to go.


%
%
% this keeps score which subject did which scan, and who should be excluded
% because of unusable data. This is assessed from a visual inspection of
% the fMRIPREP data. Excluded subjects are not enumerated in the list.
%
%

slist=struct();
slist.mov1a = [1:21]; % all subjects did movie1a and rest1a
slist.resta = [1:21];
slist.mov1b = [2:5 7:14 16:20]; % these are all subjects who returned 
slist.restb = [2:5 7:14 16:20];


% update with bad subjects -- motion was too high here. 
bads = [4 13 15];
fns=fieldnames(slist);
for i=1:numel(fns)
    for j=1:numel(bads)
        slist.(fns{i})(slist.(fns{i})==bads(j))=[];
    end
end

% then incorporate the occasioal bads (except for this one, the other scans
% look OK). 
slist.resta(slist.resta==5)=[];
slist.resta(slist.restb==6)=[];


%
%
% This is to specify removal of the first 5 scans
%
%

NULL_FIRST = 5;

%
%
% To say that there are 14 brain networks
%
%


J=14;  % the number of networks!



%
%
% Run the main loop which sorts through everything, creates directories and
% files, runs the HMM, and saves the results to the .mat file.
%
%



for i_NSTATES=1:numel(NSTATES)
    this_NSTATES=NSTATES(i_NSTATES);
    
    % the names of the time-series data are called in a systematic way;
    % this is to build the filename which should be loaded.
    covregtypes={'normal','aroma','aroma-gsr'};
    covregpres={'','aroma-','aroma-gsr-'};

    
    for i_covregtype=2 % we only run AROMA, to run all change to [1:3]
        covregtype = covregtypes{i_covregtype};
        covregpre  = covregpres{i_covregtype};
        
        for i_analysis = WHICH_ANALYSES % optionally run ALL analysis; change to [1:6]
            
            
            scan_names = analyses{i_analysis, 1};
            save_dir =  ['../results_' num2str(this_NSTATES) '/' covregtype '/' analyses{i_analysis, 2}];
            
            % use with caution, since we potentially remove directory which
            % already exists.
            if exist(save_dir,'dir')
                rmdir(save_dir,'s');
                mkdir(save_dir);
            else
                mkdir(save_dir);
            end
            
            
            % figure out which subjects need to be included. You can
            % include more subjects IF you do different analysis; for
            % example in an HMM run on mov1a + resta; you might have 18
            % total. for all, you have 14. For a subjects to be included
            % into the HMM analysis, the person needs to have done ALL the
            % scans (mov1a, resta, mov1b, restb) for the analysis 'all',
            % and should also have usable data.
            subs_to_use = [1:21]; % elements from here will be removed to figure out which subjects will be used.
            
            nsubs=21;
            for i=1:numel(scan_names)
                subs_to_use = intersect(subs_to_use, slist.(scan_names{i}));
                nsubs=numel(subs_to_use);
            end
            
            
            
            alld={}; % container for loading the timeseries
            T=[];    % container for number of time points
            
            
            for i_scan_name = 1:numel(scan_names)
                scan = scan_names{i_scan_name};

                % d=dir(['extracted_timeseries/*' scan '*r01.txt']);
                
                
                d=[]; % this build up a file list of subjects; the file being of the first roi/network.
                for i=1:nsubs
                    d = [d dir(['../data/extracted_timeseries/' covregpre 'ts-' scan '*' sprintf('s%.2d',subs_to_use(i)) '*r01.txt'])];
                end
                
                % this is to later use to set the number of scans properly.
                switch scan
                    case {'mov1a','mov1b'}
                        nscans=535;
                        
                        if TRUNCATE_VOLS_TO_220
                            scans=220;
                        end
                        
                    case {'resta','restb'}
                        nscans=220;
                end
                
                
                scans_sum=0;
                
                for i=1:numel(d)
                    
                    
                    scans_sum = scans_sum + nscans - NULL_FIRST;
                    % subj = regexp(fname,'s[0-9]{2}','match'); subj=subj{1};
                    
                    alld{end+1} = [];
                    for j=1:14
                        
                        
                        fname=regexprep(['../data/extracted_timeseries/' d(i).name],'r01.txt',sprintf('r%.2d.txt',j));
                        
                        
                        disp(fname)    % double-check what's being loaded
                        v=load(fname); % load the timeseries data
                        
                        if TRUNCATE_VOLS_TO_220
                            % truncate vols to these indices (I added 530/2
                            % but you can leave it out to select a
                            % different part of the movie
                            cut_sel_inds = (1:220) + 530/2;
                            v=v(cut_sel_inds);
                        end
                        
                        disp(length(v));
                        % discard first 5 volumes (!!!)
                        
                        
                        
                       
                        v(1:NULL_FIRST) = [];
                        
                        
                        
                        
                        dv=detrend(v,'constant'); % remove the mean
                        
                        
                        sfdv=dv/std(dv); % scale with standard deviation
                        
                        
                        alld{end} = [alld{end} sfdv]; % append the data
                    end
                    
                    
                end
                
                T(end+1) = scans_sum; % add the number of scans to the list
                
                
                
                
            end
            
            
            
            % concatenate  all loaded data -- this is the data going into
            % hmm-mar
            dat=[];
            for i=1:numel(alld)
                dat=[dat; alld{i}];
            end
            
            % smooth the loaded data
            ndat=dat;
            for i=1:14
                ndat(:,i) = smoothdata(dat(:,i));
            end
            dat=ndat;
            

            
            
            tic;
            
            n_sub = numel(subs_to_use); % no. of subjects(used for concatenation)
            K = this_NSTATES; % no. states
            reps = HMMREPS; % to run it multiple times (saves all the results
            % as a seperate mat file)
            
            TR = 2.20;  % TR of timeseries
            
            % options for model estimation; I have set the usual choices here
            % for explanation see the HMM-MAR Wiki page at:
            % https://github.com/OHBA-analysis/HMM-MAR/wiki
            options = struct();
            options.K = K; % number of states
            options.order = 0; % no autoregressive components
            options.zeromean = 0; % model the mean
            options.covtype = 'full'; % full covariance matrix
            options.Fs = 1/TR;
            options.verbose = 1;
            options.standardise = 1;
            options.inittype = 'HMM-MAR';
            options.cyc = 500;
            options.initcyc = 10;
            options.initrep = 3;
            
            % run the HMM multiple times - and save the results to .mat
            % files
            for r = 1:reps
                disp(['RUN ' num2str(r)]);
                [hmm, Gamma, ~, vpath, ~, ~, fe] = hmmmar(dat,T,options);
                save([save_dir '/HMMrun_rep_' num2str(r) '.mat'],'Gamma','vpath',...
                    'hmm','T','J','K','n_sub','fe');
                
                % keyboard;
                % we also save the data now...
                save([save_dir '/HMMrun_rep_' num2str(r) '_data.mat'], 'dat');
                
                
                % calculate summary measures for this HMM and save those
                % too to the disk:
                mean_em=zeros(J,K); % mean emissions
                for k=1:K
                    mean_em(:,k)=getMean(hmm,k);
                end
                
                prob=hmm.P; % transition probabilities
                
                
                do = 1; % do = Flag to replace 0s with NaNs (in all summary measures)
                % do=0 keeps all misisng values as zeros; do=1 sets them
                % to Nans; this is important while comparing between groups!!
                [FO,dign,avg_life]=summary_measures(vpath,n_sub,do);
                
                save([save_dir '/Summary_measures_rep_',num2str(r) '.mat'],...
                    'mean_em','prob','FO','dign','avg_life');
                
                
            end
            
            
        end
        
        
    end
    
    
end

