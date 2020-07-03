


addpath(genpath(pwd));
clear all;close all;clc;
% addpath(genpath('../hmm-necs-scripts/'))


load outr.mat

ORDER_STATES=[1 2 3 4 5 6 7 8 9 10];
NSTATES=ORDER_STATES;


RUN = 10;  % select which run to use
ANALYSIS = 'all'; % which analysis (usually 'all')
parts =  {'mov1a', 'resta'}; % which part of that analysis - compare mov with rest, session A
preprocessing='aroma';




analysis = ANALYSIS;
% parts =  {'mov1a','resta','mov1b','restb'};
% parts =  {'mov1a','mov1b'};
% parts =  {'restb','mov1b'};


% which is the optimal mapping?
% load HungarianMapping.mat

% mcorrs=[mapping(1:15).overall_corr];
% figure;plot([mapping(1:15).overall_corr]);
% best_run = mapping(find(mcorrs==max(mcorrs))).iter_val;

run=RUN;


allfo=[];
dwellt={};
for i_parts=1:numel(parts)
    part=parts{i_parts};
    
    % st=reshape(st_path,n_sub,[]); % reshape to n_sub x n_time format
    st=outr.(sprintf('run%d', RUN)).(preprocessing).(analysis).(part).vpath';
    
    st = st(:, :); % fix the subject issue... --> this will get rid of subject 5, shich had relatively bad restign state a.
    % this is to keep the analysis between rest and movie consistent in
    % who hs actually ioncluded in the anaysis...
    % st is the state path of everyone (now a 14-by-timepoints matrix).
    
    % select the right subject(s)...
    
    
    
    S=10; % Nstates = 12.
    n_sub=size(st,1);
    n_time=size(st,2);
    
    
    fo=zeros(n_sub,S);
    
    for kk=1:n_sub
        stat=st(kk,:); % getting state patih seq of individual subject
        for k=1:S
            fo(kk,k)=sum(stat==k); % counting number of time points spent in each state
            
            dwellt{kk,k,i_parts} = grab_collection_of_ones(stat==k);
            
            pauses_inbetween{kk, k, i_parts} = grab_collection_inbetweens(stat==k);
            
        end
    end
    %fo=100*fo./176%(sum(fo,2));% normailzing to get the % time spent in each state
    fo=100*fo./n_time;
    
    % because it == 0, but we don't wish to mess up the violinplot:
    % fo(:,5) = rand(n_sub,1)*0.001;
    
    allfo(:,:,i_parts)=fo;
end

% figure;
%boxplot(allfo);


% doing the stats:
m=allfo(:,:,1) - allfo(:,:,2);
r=struct();r(:)=[];for i=1:size(m,2);[H,P,CI,STATS]=ttest(m(:,i));r(end+1).H=H;r(end).P=P;r(end).CI=CI;r(end).STATS=STATS;end


addpath(genpath('Violinplot-Matlab-master'));


% fh=figure;
% violinplot(fo, repmat({'mov1a'},1,n_sub));


% make the stupid figure - 1
% for printing the network -- in COLORS...
cmat=[[230, 25, 75];...
    [60, 180, 75];...
    [255, 225, 25];...
    [0, 130, 200];...
    [245, 130, 48];...
    [70, 240, 240];...
    [240, 50, 230];...
    [250, 190, 190];...
    [0, 128, 128];...
    [230, 190, 255];...
    [170, 110, 40];...
    [255, 250, 200];...
    [128, 0, 0];...
    [170, 255, 195];...
    [0, 0, 128];...
    [128, 128, 128];...
    [255, 255, 255];...
    [0, 0, 0]]/255;



fo_for_boxplot= [];
for i=1:numel(parts)
    fo_for_boxplot=[fo_for_boxplot allfo(:,ORDER_STATES,i)];
end




labels = {};
positions=[];
box_colors=[];
for i_parts=1:numel(parts)
    positions = [positions i_parts:numel(parts):numel(parts)*S];
    part=parts{i_parts};
    for i_S=ORDER_STATES
        labels{end+1} = sprintf('S%d',i_S);
        box_colors = [box_colors; cmat(i_S,:)];
    end
    
end

box_colors=[];
for j=ORDER_STATES
    for i=1:numel(parts)
        box_colors(end+1,:)=cmat(j,:);
    end
end
set(0,'defaulttextInterpreter','none')
set(groot, 'defaultAxesTickLabelInterpreter','none'); set(groot, 'defaultLegendInterpreter','none');

fh=figure('color','w');
subplot(2,1,1); box off;

boxplot(fo_for_boxplot,'labels',labels,'positions',positions,'colors',box_colors,'labelverbosity','minor','widths',0.5);

tags={'Box','Upper Adjacent Value','Lower Adjacent Value','Upper Whisker','Lower Whisker','Median'};
le=[];for i=1:numel(tags)
    le=[le; findobj(gcf,'type','line','tag',tags{i})];
end
set(le,'linewidth',1.5)


for i=1:10
    for j=1:14
        % for k=1:2
        
        % begin of line
        byl=fo_for_boxplot(j,i);
        eyl=fo_for_boxplot(j,i+10);
        
        bxl=(i-1)*2 + 1.35;
        exl=(i-1)*2 + 1.65;
        
        if byl ~= 0 && eyl ~= 0
            
            lh=line([bxl, exl],[byl,eyl]);
            set(lh,'color',[0.6 0.6 0.6],'marker','.') %cmat(i,:));
        end
        
    end
end

% xlabel('state');
ylabel('occupancy (%)');
ylim([0 90]);
% title('all');
set(gca,'xlim',[0 numel(parts) * 10 + 1]);
box off;
set(gca,'xtick',[],'xticklabel',{});


% our stats...
for i=1:numel(r)
    if ~ isnan(r(i).H)
        if r(i).P<0.005
            text(i*2-0.5, max(get(gca,'ylim')),'*');
        end
    end
end


% % figure;
% % violinplot(fo_for_boxplot,labels);
% TR=2.2;
% % make a violinplot, then...
% violinstruct=struct(); % here we make it!!
% for i=ORDER_STATES
%     for j=1:numel(parts)
%      violinstruct.(['s' num2str(i) '_' parts{j}])=[];
%     end
% end
%
%
% % in order to do some states -- % reduce dwellt to matrix like above:
% testmatp=[];
% for i=1:size(pauses_inbetween,1)
%     for j=1:size(pauses_inbetween,2)
%         for k=1:size(pauses_inbetween,3)
%             this_vals= pauses_inbetween{i,j,k};
%             if numel(this_vals)==0 || sum(this_vals) == 0
%                 testmatp(i, j, k) = NaN;
%             else
%                 testmatp(i,j,k) = mean(this_vals);
%             end
%
%         end
%     end
% end
%
%
% % doing the stats:
% mp=testmatp(:,:,1) - testmatp(:,:,2);
% rp=struct();rp(:)=[];for i=1:size(mp,2);[H,P,CI,STATS]=ttest(mp(:,i));rp(end+1).H=H;rp(end).P=P;rp(end).CI=CI;rp(end).STATS=STATS;end




%
% for i=1:size(pauses_inbetween,1) % subjs
%     for j=1:size(pauses_inbetween,2)  % states
%         for k=1:size(pauses_inbetween,3)  % rest/mov?
%
%             for p=1:numel(pauses_inbetween{i,j,k})
%                 violinstruct.(['s' num2str(j) '_' parts{k}])(end+1) = pauses_inbetween{i,j,k}(p) * TR;
%
%             end
%
%         end
%     end
% end
%
%
% subplot(3,1,3); box off;
% vs=violinplot(violinstruct);
% ylabel('time between visits (s)');
% set(gca,'xticklabel',regexprep(get(gca,'xticklabel'),'_',' '));
% set(gca,'xticklabel',regexprep(get(gca,'xticklabel'),'mov1','mov'));

%
% for i=1:numel(vs)
%     vs(i).ViolinColor = box_colors(i,:);
% end

% keyboard


%
%
% set(gca,'xlim',[0 numel(parts) * 10 + 1]);
% xtickangle(45)
% ylim([0 350]);
% % our stats...
% for i=1:numel(rp)
%     if rp(i).H
%         text(i*2-0.5, max(get(gca,'ylim')),'*');
%
%     end
% end




%
%
% TR=2.2;
% % make a violinplot, then...
% violinstruct=struct();
% for i=ORDER_STATES
%     for j=1:numel(parts)
%      violinstruct.(['s' num2str(i) '_' parts{j}])=[];
%     end
% end


% in order to do some states -- % reduce dwellt to matrix like above:
testmat=[];
for i=1:size(dwellt,1)
    for j=1:size(dwellt,2)
        for k=1:size(dwellt,3)
            this_vals= dwellt{i,j,k};
            if numel(this_vals)==0
                testmat(i, j, k) = NaN;
            else
                testmat(i,j,k) = mean(this_vals);
            end
            
        end
    end
end

% 
% % doing the stats:
mdt=testmat(:,:,1) - testmat(:,:,2);
rdt=struct();rdt(:)=[];for i=1:size(mdt,2);[H,P,CI,STATS]=ttest(mdt(:,i));rdt(end+1).H=H;rdt(end).P=P;rdt(end).CI=CI;rdt(end).STATS=STATS;end

%
% for i=1:size(dwellt,1) % subjs
%     for j=1:size(dwellt,2)  % states
%         for k=1:size(dwellt,3)  % rest/mov?
%
%             for p=1 %:numel(dwellt{i,j,k})
%                 violinstruct.(['s' num2str(j) '_' parts{k}])(end+1) = mean(dwellt{i,j,k}(p) * TR);
%
%             end
%
%         end
%     end
% end

TR = 2.2;

% obtain for MOVIE -- av. dwell time PER SUBJECT PER STATE
% a cell,
av_dwell_per_sub_and_state=[];
% av_dwell_per_sub_and_state_rest=[];
av_dwell_per_sub_and_state_mov=[];
av_dwell_per_sub_and_state_rest=[];
dwell_for_boxplot=[];
for i=1:size(dwellt, 1) % subjs
    for j=1:size(dwellt, 2)
        
        
        k=1;
        if numel(dwellt{i,j,k}) == 0
            av_dwell_per_sub_and_state_rest(i,j) = NaN;
        else
            av_dwell_per_sub_and_state_rest(i,j) = mean(dwellt{i,j,k}) * TR;
        end
        
        k=2;
        if numel(dwellt{i,j,k}) == 0
            av_dwell_per_sub_and_state_mov(i,j) = NaN;
        else
            av_dwell_per_sub_and_state_mov(i,j) = mean(dwellt{i,j,k}) * TR;
        end
        
    end
end
tmp=[av_dwell_per_sub_and_state_rest av_dwell_per_sub_and_state_mov];
dwell_for_boxplot = tmp(:,[1 11 2 12 3 13 4 14 5 15 6 16 7 17 8 18 9 19 10 20]);



subplot(2,1,2); box off;

labels={'rest s1','mov s1','rest s2','mov s2','rest s3','mov s3','rest s4','mov s4','rest s5','mov s5','rest s6','mov s6','rest s7','mov s7','rest s8','mov s8','rest s9','mov s9','rest s10','mov s10'}; % ','rest

boxplot(dwell_for_boxplot,'labels',labels,'positions',[1:20],'colors',box_colors,'labelverbosity','minor');

tags={'Box','Upper Adjacent Value','Lower Adjacent Value','Upper Whisker','Lower Whisker','Median'};
le=[];for i=1:numel(tags)
le=[le; findobj(gcf,'type','line','tag',tags{i})];
end
set(le,'linewidth',1.5)

% xlabel('state');
ylabel('dwell time (s)');
% title('all');
set(gca,'xlim',[0 numel(parts) * 10 + 1]);
box off;
% set(gca,'xtick',[],'xticklabel',{});



% let's make some lines

for i=1:10
    for j=1:14
        % for k=1:2
   
        % begin of line
        byl=dwell_for_boxplot(j,(i-1)*2+1);
        eyl=dwell_for_boxplot(j,(i-1)*2+2);
        
        bxl=(i-1)*2 + 1.35;
        exl=(i-1)*2 + 1.65;
        
        lh=line([bxl, exl],[byl,eyl]);
        set(lh,'color',[0.6 0.6 0.6],'marker','.'); % cmat(i,:));
        
    end
end

%
% subplot(3,1,2); box off;
% vs=violinplot(violinstruct);
% ylabel('dwell time (s)');


set(gca,'xlim',[0 numel(parts) * 10 + 1]);
xtickangle(45)
ylim([0 50]);
% set(gca,'xtick',[],'xticklabel',{});
% 
% our stats...
for i=1:numel(rdt)
    if ~isnan(rdt(i).H)
        if rdt(i).P < 0.005
            text(i*2-0.5, max(get(gca,'ylim')),'*');
        end
    end
end
box off


set(fh,'paperunits','centimeters');
set(fh,'papersize',[22 15],'paperposition',[0 0 22 15]);

% saveas(fh, 'occdwell_ses1.jpg');
% fsource = 'f5_occdwell_new1_sesa.jpg';
fsource = ['../figures/Fig4-' preprocessing '-' ANALYSIS '-occdwell_' sprintf('run%d',run) '.pdf'];

% fsource2 = 'f5_occdwell_new1_sesa.pdf';
% ftarget_file=fsource;
saveas(fh,fsource);


% create the Table...
 % Matlab tables are definitely not pandas dataframes:
Participant = {};
for i=1:size(dwell_for_boxplot,1)
    Participant{end+1} = sprintf('Participant_%d',i); 
end
Participant=Participant';

% the FO has columns:
s='';
for i=1:NSTATES
    for j=1:numel(parts)
        s=[s sprintf('FO_%s_S%d = fo_for_boxplot(:,%d)\n',parts{j}, i, i+(j-1)*10)];
    end
end

% the DT also has columns, but different order than the FO's...
% s='';
for i=1:NSTATES
    for j=1:numel(parts)
        s=[s sprintf('DT_%s_S%d = dwell_for_boxplot(:,%d)\n',parts{j}, i, 2*i-1+j-1)];
    end
end
eval(s)

s='';
s=[s sprintf('myTable = table(')];
s=[s sprintf('Participant,')];

mytypes={'FO','DT'};
for k=1:numel(mytypes)
    for i=1:NSTATES
        for j=1:numel(parts)
            s=[s sprintf('%s_%s_S%d,',mytypes{k},parts{j}, i)];
        end
    end
end
s(end)=[]; % remove comma
s=[s sprintf(');')];

eval(s)

% fsource = ['../figures/Fig4-' preprocessing '-' ANALYSIS' '-occdwell_' run '.pdf'];
fsource_xls = ['../figures/Fig4-' preprocessing '-' ANALYSIS '-occdwell_' sprintf('run%d',run) '_values.xls'];

writetable(myTable,fsource_xls);






% saveas(fh,fsource2);
% 
% 
% 
% 
% 
% ftarget1=['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' ftarget_file];
% ftarget2=['/home/johan/ldrive/Lab_LucaC/10_Johan/1_HMMMovie/' ftarget_file];
% 
% print('-djpeg','-r600','/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/Drawings_Inkscape/F5_FO_DT_sa');
% 
% copyfile(fsource,ftarget1);
% copyfile(fsource,ftarget2);
% 
% saveas(fh, ['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' fsource2]);
% saveas(fh, ['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' fsource]);
% 
% 
% delete(fsource);


% cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript





%
%
%
% for i=1:12;
%     % [p,c,stat]=ttest2(ADHD(:,i),CTLS(:,i))
%     % so it'd be possible to calculate the fo as a matrix for two
%     % situations -- and then perform some stats on them.
%     % pair-wise stats...
%
% end