% f4_state_ephys_annots

% addpath('/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript');
addpath(genpath('hmm_new'));
% addpath('/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/hmm-necs-scripts');




% so clean up the following code a little bit -- then allow selection of N
% states - then allow determination of task, too.

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
    [0, 0, 0]];


%
%
% all this to extract subs_to_use -- which subject we will use in this
%
%
%

% analysis (i.e., 'all').
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


analyses = {...
    {'mov1a'},'mov1a'; ... % only mov1
    {'mov2a'},'mov2a'; ... % only mov2 (repeat on d1 from 11min onwards)
    {'mov1a','mov1b'},'mov1a-1b'; ... % mov1 day and and day 2
    {'resta'},'resta'; ...  % only resta
    {'restb'},'restb'; ...  % only rest day 2
    {'resta','restb'},'resta-b'; ...  % rest a and rest day 2
    {'mov1a','mov1b','resta','restb'},'all'; ...  % all movie and rest day 1 aand day 2 but not repeats
    };


scan_names = analyses{7, 1};

analysis = analyses{7, 2}(1:end-4);
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




%
%
%  some fuctions that would be useful later on:
%
%
%




% so there are 2 figures to be made

mydice = @(v1, v2) 2*sum(v1.*v2)/(sum(v1)+sum(v2));
myshuffle = @(v) v(randperm(numel(v)));
% we also use, for our plotting, my_create_fig...




%
%
% load the data generated while we did fig. 2: it's 'sd', and has many many
% fields.
%
%
%
%




% the first is the subject-specific figure - just (simply) showing the full
% data available for each subject - that's how it is also stored (and
% recoverable)

load sd.mat




% get 'all' data (i.e. mov1a+1b+resta+restb, skipping mov2a) - preprocessed
% with aroma, hmm run no. 13.

load valid_inferences_all.mat
this_inference = valid_inferences(1);

% matlab dynamic fieldnames written up as an sprinted string...
s=sd(1).aroma.all.(sprintf('sm%d',this_inference));



%

NPERMS=500;
% whichsessions = {'resta','restb','mov1a','mov1b'};
% whichsessions={};
STATES=[1 2 3 4 5 6 7 8 9 10];
bigfs=[];

ds=struct();


sd2=struct();

anno_mat = load('../data/annotations/annobigmat.mat');
annobigmat = anno_mat.bigmat(6:end,:); % getting rid of 1st 5 volumes



% we set our thing to mov1a + mov1b - and concatenate?

XLIM = [0.055  0.145];

% MODALITIES={'hr','eye','gsr'};
highlows = {'high','low'};  % for >20% and <20% -- for now..
ANNOTATIONS = {'L','CP','F+','F-','S+','S-'};
ANNOTATIONS_for_struct = {'L','CP','Fpos','Fneg','Spos','Sneg'};
ANNOTATIONStitles = {'language','changepoint','faces +','faces -','scenes +','scenes -'};
ANNOTATIONStitles=ANNOTATIONS;

bigf=figure;
bigfs(end+1) = bigf; % save figure handle...
set(bigf,'color','w');

fcols=numel(ANNOTATIONS);
frows=numel(STATES);

spacingx=0.02;
spacingy=0.02;

% borders of our image:
bigspacingxl = 0.20;
bigspacingyl = 0.10;

bigspacingxu = 0.10;
bigspacingyu = 0.15;

factorx = (1-bigspacingxl-bigspacingxu-(fcols-1)*spacingx)/fcols;
factory = (1-bigspacingyl-bigspacingyu-(frows-1)*spacingy)/frows;


new_ahs=[];
my_titles={};
for iState=1:numel(STATES)
    state=STATES(iState);
    
    
    % we ALSO (!) loop over the kind of electrophysiological modality
    for iAnnotation = 1:numel(ANNOTATIONS)
        annotation = ANNOTATIONS{iAnnotation};
        
        
        
        my_titles{iState,iAnnotation} = sprintf('%s',ANNOTATIONStitles{iAnnotation});
        
        % now figure out which subjects have this modality:
        collected={};
        
        % tmp=[];for i=1:21;tmp(end+1) = 0;if isfield(s,sprintf('s%d',i));if isfield(s.(sprintf('s%d',i)),whichsess); tmp(end)=isfield(s.(sprintf('s%d',i)).(whichsess),sprintf('%smat',modality));end;end;end;subs=find(tmp);
        
        subs_for_this = subs_to_use;
        
        
        % get the state vectors - for each sub: - put into matrix
        statemat=[];
        annotmat=[];
        excluded_subs=[];
        included_subs=[];
        
        
        
        for isubj=1:numel(subs_for_this)
            if sum(isnan(s.(sprintf('s%d',subs_for_this(isubj))).('mov1a').hmmmat(:,state)))==0 && sum(isnan(s.(sprintf('s%d',subs_for_this(isubj))).('mov1b').hmmmat(:,state)))==0
                % disp('go');
                included_subs(end+1) = subs_for_this(isubj);
                statemat(end+1,:) = [s.(sprintf('s%d',subs_for_this(isubj))).('mov1a').hmmmat(:,state); s.(sprintf('s%d',subs_for_this(isubj))).('mov1b').hmmmat(:,state)];
                annotmat(end+1,:) = [annobigmat(:,iAnnotation); annobigmat(:,iAnnotation) ];
            else
                excluded_subs(end+1) = subs_for_this(isubj);
            end
            
        end
        
        
        
        
        ah=[];
        im=[];
        
        old_xp = 0; % old_position(1);
        old_yp = 0; % old_position(2);
        old_xs = 1; % old_position(3);
        old_ys = 1; % old_position(4);
        
        x_offset = 0;
        
        x_offset = x_offset + (factorx+spacingx)*(iAnnotation-1);
        
        % x_offset = x_offset + factorx + spacingx;
        fj = iState;
        
        new_xp = old_xp * factorx + x_offset + bigspacingxl;
        new_yp = 1-(old_yp * factory + (fj) * (factory + spacingy) + bigspacingyl);
        new_xs = old_xs * factorx;
        new_ys = old_ys * factory;
        
        new_position = [new_xp new_yp new_xs new_ys];
        
        
        % ah(end+1) = axes('parent',bigf,'position',new_position);
        
        % now - perform all of the permutation magic...
        
        
        % set(gca,'xlim',XLIM);
        ylims=[0 660];
        xlims=[0.055 0.125]-0.01;
        centroids = xlims(1):diff(xlims)*0.01:xlims(end);
        xlims=[0.055 0.125]-0.01;
        centerit = 1;
        if iState==1
            print_n=0;
        else
            print_n=0;
        end
        
        NPERMS=500;
        centerit = 1;
        [fh, h, p_05, p_95, p_real,  pseudo_z_score, vd_mean, vd_std] = my_permutations_3(statemat,annotmat, included_subs,xlims,ylims,centroids,centerit,print_n,NPERMS);
        
        sd2(1).(sprintf('S%d',iState)).(ANNOTATIONS_for_struct{iAnnotation}).p_05 = p_05;
        sd2(1).(sprintf('S%d',iState)).(ANNOTATIONS_for_struct{iAnnotation}).p_95 = p_95;
        sd2(1).(sprintf('S%d',iState)).(ANNOTATIONS_for_struct{iAnnotation}).p_real = p_real;
        sd2(1).(sprintf('S%d',iState)).(ANNOTATIONS_for_struct{iAnnotation}).pseudo_z_score = pseudo_z_score;
        sd2(1).(sprintf('S%d',iState)).(ANNOTATIONS_for_struct{iAnnotation}).vd_mean = vd_mean;
        sd2(1).(sprintf('S%d',iState)).(ANNOTATIONS_for_struct{iAnnotation}).vd_std = vd_std;
        
        
        to_copy=get(fh,'children');
        new_ax=copyobj(to_copy,bigf);
        set(new_ax,'position',new_position);
        new_ahs(iState,iAnnotation) = new_ax;
        
        close(fh);
        
        
        
    end
    
end


% here I can do whatever..
for i=1:size(new_ahs,1)
    for j=1:size(new_ahs,2)
        
        %if i~=size(new_ahs,1)
        %    set(new_ahs(i,j),'xticklabel',{},'xtick',[]);
        %end
        
        set(new_ahs(i,j),'fontsize',6);
        
        if j~=1
            set(new_ahs(i,j),'yticklabel',{},'ytick',[]);
        end
        
        if i==1
            title(new_ahs(i,j),my_titles{i,j});
        end
        
        if j==1
            
            thxl=get(new_ahs(i,j),'xlim');
            thyl=get(new_ahs(i,j),'ylim');
            state=STATES(i);
            
            myXPOS=thxl(1) - diff(thxl) * 0.25;
            myYPOS=mean(thyl)+diff(thyl)*0.2;
            
            th3=text(myXPOS-diff(get(new_ahs(i,j),'xlim'))*0.0005,myYPOS,sprintf('State %d',state),'fontweight','bold','verticalalignment','top','horizontalalignment','right','parent',new_ahs(i,j));
            th3=text(myXPOS+diff(get(new_ahs(i,j),'xlim'))*0.0005,myYPOS,sprintf('State %d',state),'fontweight','bold','verticalalignment','top','horizontalalignment','right','parent',new_ahs(i,j));
            th3=text(myXPOS,myYPOS-diff(get(new_ahs(i,j),'ylim'))*0.0005,sprintf('State %d',state),'fontweight','bold','verticalalignment','top','horizontalalignment','right','parent',new_ahs(i,j));
            th3=text(myXPOS,myYPOS+diff(get(new_ahs(i,j),'xlim'))*0.0005,sprintf('State %d',state),'fontweight','bold','verticalalignment','top','horizontalalignment','right','parent',new_ahs(i,j));
            th3=text(myXPOS,myYPOS+diff(get(new_ahs(i,j),'xlim'))*0.0005,sprintf('State %d',state),'fontweight','bold','verticalalignment','top','horizontalalignment','right','parent',new_ahs(i,j),'color',cmat(state,:)/255);
            %th=text(thxl(1),thyl(2),sprintf(' Network %d',state),'fontweight','bold','Color',cmat(state,:)/255);
            
            % th2=text(new_ahs(i,j),thxl(1),thyl(2),sprintf('\n\n%s',network_descriptions{state}),'verticalalignment','top');
            
        end
        
        
    end
end
anno_ax = axes('parent',bigf,'position',[0 0 1 1]);
set(anno_ax,'visible','off');
th=text(0.5, 0.99,sprintf('HMM States ~ Annotations'));
set(th,'horizontalalignment','center','verticalalignment','top');




%%
% and here should come -- the annotations, in very similar fashion as above
% - but then different.
% 
% 
% set(newbigf,'paperunits','centimeters');
% set(newbigf,'papersize',[30 25],'paperposition',[0 0 30 25]);
% 
% fsource = 'f11_radio_annots_electrophys.jpg';
% ftarget_file = fsource;
% saveas(newbigf,fsource);
% 
% 
% ftarget1=['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' ftarget_file];
% % ftarget2=['/home/johan/ldrive/Lab_LucaC/10_Johan/1_HMMMovie/' ftarget_file];
% 
% copyfile(fsource,ftarget1);
% % copyfile(fsource,ftarget2);
% 
% delete(fsource);


% cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript




print('-djpeg','-r300','../figures/states-annotations.jpg');



% I do it with a string because I prefer to no longer use the Table
% constructor with Matlab. It usage gives me a headache.

s=sprintf('Z-scores of Brain State Paths ~ Annotations; mov session A and B are concatenated\n');
s=[s sprintf('State:\t')];
for i=1:10
    s=[s sprintf('S%d\t',i)];
end
s=[s sprintf('\n')];

for j=1:numel(ANNOTATIONS)
    this_annot = ANNOTATIONS{j};
    this_annot_struct = ANNOTATIONS_for_struct{j}; % matlab doesn't support any string for dyn fieldnames, unlike Python. a '+' is not allowed.
    s=[s sprintf('%s\t',this_annot)];
    for i=1:10
        s=[s sprintf('%.2f\t', sd2.(sprintf('S%d',i)).(sprintf(this_annot_struct)).pseudo_z_score)];
    end
    s=[s sprintf('\n')];
end
s(end)=[];

fprintf(s);
fprintf('\n\n');




