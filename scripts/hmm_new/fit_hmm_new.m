clear;clc;tic;
addpath(genpath('HMM-MAR-master/'));

% load your data -- single subject? - there is an option to specify segment
% lengths -- where is that? Are all segments equal length?
dat = rand(100,3); % random data jst as an example
% n_time points X n_networks format

[T,J] = size(dat); % T = number of time points in concatenated data
                 % J = number of networks if method is same as ours
                 
n_sub = 21; % no. of subjects(used for concatenation)
K = 12; % no. states
reps = 5; % to run it multiple times (saves all the results
          % as a seperate mat file)

TR = 2.20;  % TR of rest data

% options for model estimation; I have set the usual choices here; 
% no need to change anything in geenral 
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

% run the HMM multiple times
for r = 1:reps
    disp(['RUN ' num2str(r)]);
    [hmm, Gamma, ~, vpath] = hmmmar(dat,T,options);
    save(['HMMrun_rep_' num2str(r) '.mat'],'Gamma','vpath',...
          'hmm','T','J','K','n_sub');   
end
