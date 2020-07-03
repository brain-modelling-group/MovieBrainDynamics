addpath('/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/revisions/R12/HMM-MAR');
addpath(genpath('/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/revisions/R12/HMM-MAR/'))

preprocessing='aroma';
analysis='all';
run='run1241';

% grab the HMM stuff:
hmmstruct=load(sprintf('../results_10/aroma/resta-b/HMMrun_rep_%s.mat',run(4:end)));
hmm=hmmstruct.hmm;

hmm.train.Gamma_constraint=[];

Gamma=hmmstruct.Gamma;
T=hmmstruct.T;

hmmdata=load(sprintf('../results_10/aroma/resta-b/HMMrun_rep_%s_data.mat',run(4:end)));
dat=hmmdata.dat;

% load HMM.mat
% load Gamma.mat
% load T.mat
% load data.mat

% this will give us hmm, T and data. Gamma is optional! So we won't use it.
% No idea how to load in Xi. Yet.



[Gamma, Xi] = hmmdecode(dat, T, hmm, 0);
vpath = hmmdecode(dat, T, hmm, 1);

% ...HMM-MAR/wiki/User-Guide
fe = hmmfe(dat, T, hmm);

% Then according to website:

hmm.train.Gamma_constraint=[];


[Gamma, Xi] = hmmdecode(dat, T, hmm, 0);
vpath = hmmdecode(dat, T, hmm, 1);

% ...HMM-MAR/wiki/User-Guide
fe = hmmfe(dat, T, hmm);

% it has a SECOND OUTPUT ARGUMENT !?
% Downloaded the new version from github...
% addpath('../../hmm_new/HMM-MAR-master/')

% we need to set this, I think.
hmm.train.lowrank=0;

[fe, likelihood_of_the_model_per_time_point] = hmmfe(dat, T, hmm);

% got error hmm.train.lowrank... what's the defualt? -- estimate hmm again
% with new toolbox...
% re-shape according to our sequence. Which was MA, MB, RA, RB; all subs
% concatenated... Check Step1....

% make for each sub (!) their own likelihood trace...

sub_points_to_extract  =[];
for i_sub=1:14
    
    points_to_extract = [];
    points_to_extract = [points_to_extract (1:215) + (i_sub-1)*215];  % mov a
    points_to_extract = [points_to_extract (1:215) + (i_sub-1)*215 + 215*14]; % mov b
    %points_to_extract = [points_to_extract (1:215) + (i_sub-1)*215 + (215 + 215)*14]; % mov b
    %points_to_extract = [points_to_extract (1:215) + (i_sub-1)*215 + (215 + 215 + 215)*14]; % mov b
    % points_to_extract = points_to_extract + (i_sub-1) * (215+215+530+530);
    
    sub_points_to_extract(:, end+1) = points_to_extract;
end

% make for each sub the fe
likelihood_per_sub = [];
for i_sub=1:14
    
    likelihood_per_sub(:, end+1) = likelihood_of_the_model_per_time_point(sub_points_to_extract(:, i_sub));
end

% make it for mova, b, resta  and rest b
subject_means = mean(likelihood_per_sub);
grand_mean = mean(likelihood_per_sub(:));
sample_means = mean(likelihood_per_sub,2);
sample_sds = std(likelihood_per_sub,[],2);
grand_sd = std(likelihood_per_sub(:));


t_13_0975 = 2.160;
ci_grand_mean = grand_mean + [1, -1] * t_13_0975 * grand_sd / sqrt(14);

new_data = likelihood_per_sub - subject_means + grand_mean;


% mean_ll_per_sub = [];
% median_ll_per_sub = [];
% 
% for i_sub = [1:14]
%     figure(i_sub);
%     
%     subplot(1,4,1);
%     poffset=0;
%     plot(new_data((1:530)+poffset,i_sub));
%     ylims = get(gca,'ylim');
%     
%     % calculate the CI here:
%     
%     title(sprintf('%s, %s', 'mova',sprintf('sub%d',i_sub)));
%     mean_ll_per_sub(i_sub, 1) = mean(new_data((1:530)+poffset,i_sub));
%     median_ll_per_sub(i_sub, 1) = median(new_data((1:530)+poffset,i_sub));
%     
%     subplot(1,4,2);
%     poffset=530;
%     plot(new_data((1:530)+poffset,i_sub));
%     set(gca,'ylim', ylims);
%     title(sprintf('%s', 'movb'));
%     mean_ll_per_sub(i_sub, 2) = mean(new_data((1:530)+poffset,i_sub));
%     median_ll_per_sub(i_sub, 2) = median(new_data((1:530)+poffset,i_sub));
%     
%     subplot(1,4,3);
%     poffset=530*2;
%     plot(new_data((1:215)+poffset,i_sub));
%     set(gca,'ylim', ylims);
%     title(sprintf('%s', 'resta'));
%     mean_ll_per_sub(i_sub, 3) = mean(new_data((1:215)+poffset,i_sub));
%     median_ll_per_sub(i_sub, 3) = median(new_data((1:215)+poffset,i_sub));
%     
%     
%     subplot(1,4,4);
%     poffset=530*2+215;
%     plot(new_data((1:215)+poffset,i_sub));
%     set(gca,'ylim', ylims);
%     title(sprintf('%s', 'restb'));
%     mean_ll_per_sub(i_sub, 4) = mean(new_data((1:215)+poffset,i_sub));
%     median_ll_per_sub(i_sub, 4) = median(new_data((1:215)+poffset,i_sub));
%     
% end
% 
% % ttests...
% % mov1 vs rest1
% [H,P,CI,STATS] = ttest2(new_data(:, 1), mean_ll_per_sub(:, 3));
% [H,P,CI,STATS] = ttest2(new_data(:, 1), mean_ll_per_sub(:, 3));


% mov2 vs rest2
% [H,P,CI,STATS] = ttest2(mean_ll_per_sub(:, 2), mean_ll_per_sub(:, 4));
% [H,P,CI,STATS] = ttest2(median_ll_per_sub(:, 2), mean_ll_per_sub(:, 4));




%
%
% n_sub = 14; %numel(subs_to_use); % no. of subjects(used for concatenation)
% K = 10;  %this_NSTATES; % no. states
% reps = 15; % to run it multiple times (saves all the results
% % as a seperate mat file)
%
% TR = 2.20;  % TR of timeseries
%
% % options for model estimation; I have set the usual choices here
% % for explanation see the HMM-MAR Wiki page at:
% % https://github.com/OHBA-analysis/HMM-MAR/wiki
% options = struct();
% make the fig, with averages...
% close(15);
figure(15);
inds={1:215, 215+(1:215)}; %, 215+215+(1:215), 215+215+215+(1:215)};
seq = [1 2];
titles={'resta','restb'}; %,'resta','restb'};

for si=1:2
    i=seq(si);
    ph(si) = subplot(1,2,si);
    % ci_data = mean(new_data(inds{i},:),2) +  std(new_data(inds{i},:),[], 2) * t_13_0975 / sqrt(14) * sqrt(numel(inds{i})/(numel(inds{i})-1)) * [1, -1];
    % ci_data = mean(new_data(inds{i},:),2) +  std(new_data(inds{i},:),[], 2) * t_13_0975 / sqrt(14) * sqrt(4/3) * [1, -1];
    ci_data = mean(new_data(inds{i},:),2) +  std(new_data(inds{i},:),[], 2) * t_13_0975 / sqrt(14) * sqrt(2/1) * [1, -1];
    
    plot(ci_data(:, 1),'color',[0.7 0.7 0.7]);hold on;
    plot(ci_data(:, 2),'color',[0.7 0.7 0.7]);
    
    plot(mean(new_data(inds{i},:),2),'k', 'linewidth',2);hold on;
    ylim([-22, -8]);
    ylims = get(gca,'ylim');
    title(sprintf('%s', titles{i}));
    xlim([1 numel(inds{i})]);
end



% % rescaling operation according to length...
% % we just change the 'position' parameter (length and offset) according to
% % the xlim
% 
% min_x = get(ph(1), 'position') * [1 0 0 0]';
% max_x = get(ph(2), 'position') * [1 0 1 0]';
% total_els = 215+215; %+215+215;
% norm_els = total_els / 4;
% 
% spacing = get(ph(2), 'position') * [1 0 0 0]' - get(ph(1), 'position') * [1 0 1 0]';


% new_offset = min_x;
% for i_ph = 1:2
%     
%     
%     old_x_offset = get(ph(i_ph), 'position') * [1 0 0 0]';
%     old_x_size   = get(ph(i_ph), 'position') * [0 0 1 0]';
%     
%     this_els = max(get(ph(i_ph), 'xlim'));
%     reshape_factor = this_els / norm_els;
%     
%     new_x_size = old_x_size * reshape_factor;
%     old_ph_position = get(ph(i_ph), 'position');
%     set(ph(i_ph), 'position', [new_offset, old_ph_position(2),new_x_size , old_ph_position(4)]);
%     
%     
%     new_offset = new_offset + spacing + new_x_size;
%     
% end
% options.K = K; % number of states
% options.order = 0; % no autoregressive components
% options.zeromean = 0; % model the mean
% options.covtype = 'full'; % full covariance matrix
% options.Fs = 1/TR;
% options.verbose = 1;
% options.standardise = 1;
% options.inittype = 'HMM-MAR';
% options.cyc = 500;
% options.initcyc = 10;
% options.initrep = 3;
%
% options = checkoptions(options)
%
%
% [hmm, Gamma, ~, vpath, ~, ~, fe] = hmmmar(dat,T,options);

fsource = ['Check_FE_' '8_minutes_run' '.pdf'];

saveas(gcf,fsource);




% for full comprehensiveness --> do it as a bargraph...
% grab the data - again...


operation = @median;
% inds={1:530, 530+(1:530), 530+530+(1:215), 530+530+215+(1:215)};

inds={1:215, 215+(1:215)}; % 215+215+(1:215), 215+215+215+(1:215)};
datmatrix = zeros(14, 2);
reorder = [1 2]; %2 4 1 3];
for i=1:2
    for j=1:14
        datmatrix(j, reorder(i)) = operation(likelihood_per_sub(inds{i}, j));
    end
end

grand_mean = mean(datmatrix(:));
subj_means = mean(datmatrix(:), 2);

new_data_2 = datmatrix - subject_means' + grand_mean;

% ci = mean(new_data_2)' + (t_13_0975 * std(new_data_2) / sqrt(14) * sqrt(4/3))' * [1, -1];
ci = mean(new_data_2)' + (t_13_0975 * std(new_data_2) / sqrt(14) * sqrt(2/1))' * [1, -1];

figure;errorbar(mean(new_data_2), (t_13_0975 * std(new_data_2) / sqrt(14) * sqrt(4/3))')
ylim([-22, -8]);
xlim([0 5]);

% set(gca,'xlim',[0, 5])
% set(gca,'ylim',[-20, -10])


