% f6_transitions
% re-do it for sessionB in a separate script? - which involves Step6a and
% Step6b AGAIN.



% cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/revisions_rerun_hmm_10_times_only_rest/revisions/scripts

% addpath(genpath(pwd));
clear all;close all;clc;
addpath(genpath(pwd));
% addpath('../../hmm-mar-matlab/4_transition_probabilities/')

% addpath /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/hmm-necs-scripts
% addpath /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/revisions_rerun_hmm_10_times/revisions/R14


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


WHICHSTATES=[1 2 3 4 5 6 7 8 9 10];
% WHICHSTATES=[1 2 3 4 6 7 ];


% obtained from calculate_nbs_stuff.m in hmm-mar-matlab/hmm-necs-scripts folder






NodeLabels={};for i=1:numel(WHICHSTATES); NodeLabels{i}=num2str(WHICHSTATES(i));end



% hope to not have to use this...
REMOVE_STATE = [];
load outr.mat


preprocessing='aroma';
analysis = 'all';
% parts =  {'mov1a','mov1b','resta','restb'};
parts =  {'mov1a','resta'};
RUN=10;
NSTATES=10;
 %ANALYSIS=analysis;



run=sprintf('run%d', RUN);
% subs_to_use_real = [2     3     7     8     9    10    11    12    14    16    17    18    19    20];

% prev_subs = outr.(run).aroma.(analysis).subs; new_inds=[];
% for i=1:numel(prev_subs)
% 	if any(prev_subs(i) == subs_to_use_real)
%         new_inds(end+1) = i;
%     end
% end

m=load('Step6a_calculate_NBS_mov_rest_matrix_sesa.mat');
nbs_results.(parts{1}) = m.m;
nbs_titles.(parts{1})='NBS Movie > Rest session A';

m=load('Step6c_calculate_NBS_mov_rest_matrix_sesb.mat');
nbs_results.(parts{2}) = m.m;
nbs_titles.(parts{2})='NBS Rest > Movie Session A';


savedir = ['../results_10/' preprocessing '/' regexprep(analysis, '_', '-') '/'];

dns=0.20;




matrix_fhs = [];
bigfs=[];

allfo=[];
collected_handles=[];
firstone=1;
for i_parts=1:numel(parts)
    part=parts{i_parts};
    
    
    
    % part=parts{i_parts};
    
    % st=reshape(st_path,n_sub,[]); % reshape to n_sub x n_time format
    st=outr.(run).(preprocessing).(analysis).(part).vpath';
    st=st(:, :);
    
    
    % should we do something here regarding 0 values of the transition for
    % (some) subjects, maybe??
    %
    %
    % % this should probably be.. checked!
    %
    %
    
    mtx = mean(outr.(run).(preprocessing).(analysis).(part).emp(:,:,:),3);

    
    pmtx=mtx;
    mtx(logical(eye(size(pmtx))))=0; %set the diagonal to zero
    ptransp = 1-logical(eye(size(pmtx)));
    fh=figure('color','w');
    ah=axes;
    myim = imagesc(pmtx,'AlphaData',ptransp);
    
    set(gca,'visible','off')
    set(gca,'clim',[0 0.09]);
    
    th=text(mean(get(gca,'xlim')),min(get(gca,'ylim')),'To State','verticalalignment','bottom','horizontalalignment','center','fontsize',12);
    th=text(min(get(gca,'xlim')),mean(get(gca,'ylim')),'From State','verticalalignment','bottom','horizontalalignment','center','fontsize',12);
    set(th,'rotation',90);
    th=text(min(get(gca,'xlim')),min(get(gca,'ylim')),part,'verticalalignment','bottom','horizontalalignment','right','fontsize',12);
    
    % cb=colorbar('SouthOutside');
    % cb=colorbar('South');
    % title(part);
    
    matrix_fhs(end+1) = fh;
    
    
    % make the excel (ugh, excel..)
    for ii=1:NSTATES
        s=sprintf('ToState%d=pmtx(:,%d)',ii,ii);
        eval(s);
    end
    
    s=sprintf('x = {');
    for ii=1:NSTATES
        s=[s sprintf('\''FromState%d\'';',ii)];
    end
    s(end)=[];
    s=[s sprintf('}')];
    eval(s);
    
    s=sprintf('myTable=table(x,');
    for ii=1:NSTATES
        s=[s sprintf('ToState%d,',ii)];
    end
    s(end)=[];
    s=[s sprintf(')')];
    eval(s);
    

    % savefig_title...
    savefig_i = sprintf('../figures/Fig5_NBSFigure_StateTransitions_%s_data.xls',part);
    writetable(myTable, savefig_i);
    
    

    
    
    mtx=mtx(WHICHSTATES,:);
    mtx=mtx(:,WHICHSTATES);
    
    
    % remove state 5...
    % mtx(REMOVE_STATE,:)=[];
    % mtx(:,REMOVE_STATE)=[];
    
    
    raw_mtx=mtx;
    
    
    wts=[];for i=1:10;wts(i) = sum(st(:)==i);end;wts=wts/sum(wts);
    my_wts = wts*10*10;
    my_wts=my_wts(WHICHSTATES);
    % my_wts(REMOVE_STATE)=[];
    
    
    
    
    %mtx=mCTLS_Trans_prob;
    % dns=1;% density. % of connections that needs to be seen 0.1 = top 10% connections are displayed
    %hf=figure;
    %hf.Color='w';
    %hf.Position=[100 100 1000 500];
    
    
    
    mtx(logical(eye(size(mtx))))=0; %set the diagonal to zero
    [id,od,deg] = degrees_dir(mtx>=0.001);  % median of the transition values --> ONLY THE TRANSITION VALUES
                                            % this will effectively
                                            % dis-regard the trace, which
                                            % contains only 0's anyway
                                            % it'll also disregard nonsense
                                            % transition values.
    
                                            
                                            
                                            
    
    %tmp=threshold_proportional(mtx,dns);
    % G=digraph(tmp);
    
    
    
    % clear;clc; close all; tic;
    
    % con_avg_prb=rand(12,12); % the trans. prob. mtx that needs to be visulaized
    % avearge transitions in controls
    
    % dns=0.15;% density.
    % prcntage of connections that needs to be seen
    % 0.1 = top 10% connections are displayed
    
    hf=figure;hf.Color='w';
    bigfs(end+1) = hf;
    
    
    hf.Position=[100 100 1000 500];
    
    tmp=threshold_proportional(raw_mtx,dns);
    % thresholding trans_prob mtx
    
    G=digraph(tmp);
    
    
    % make the excel (ugh, excel..)
    for ii=1:NSTATES
        s=sprintf('ToState%d=tmp(:,%d)',ii,ii);
        eval(s);
    end
    
    s=sprintf('x = {');
    for ii=1:NSTATES
        s=[s sprintf('\''FromState%d\'';',ii)];
    end
    s(end)=[];
    s=[s sprintf('}')];
    eval(s);
    
    s=sprintf('myTable=table(x,');
    for ii=1:NSTATES
        s=[s sprintf('ToState%d,',ii)];
    end
    s(end)=[];
    s=[s sprintf(')')];
    eval(s);
    

    % savefig_title...
    savefig_ii = sprintf('../figures/Fig5_NBSFigure_max20percent_%s_data.xls',part);
    writetable(myTable, savefig_ii);
    
    
    
    
    w=G.Edges.Weight; w=w/max(w); w=w*5; % scaling the edges
    % you might need to vary 5 depending on your figure. Increase if edges are
    % thin and decrease if they appear very thick.
    LWidths=w;
    h=plot(G,'LineWidth',LWidths,'Layout','layered');
    % h=plot(G,'LineWidth',LWidths,'Layout','layered','XData',my_xdata,'YData',my_ydata,'ZData',my_zdata);
    
    title(sprintf('%s', regexprep(part,'1','')),'Interpreter','none');
    axis off; h.ArrowSize=10; % changes the arrow size
    % again might need to vary arrow size
    whos t
    h.NodeColor=cmat(WHICHSTATES,:); % set node color to green
    h.NodeLabel=NodeLabels;
    
    % setting nodes 1,2 &5 to cyan color
    % this might hep to distinguish between sattes that have
    % significantly different Fo values between groups
    % highlight(h,[1 2 5],'NodeColor','cy');
    
    % h.MarkerSize=4*(1:12); % changing node size
    h.MarkerSize=0.25*round(my_wts)+4; % changing node size
    % since our trans_prob is 12x12, we have 12 nodes and I'm just
    % making them increase in size from 1 to 12.
    % in practice, you can chnage this according toavg state FO values;
    % so that states occupied most appear larger
    
    % for customising the plot further, look at the properties of object 'h'
    
    collected_handles(end+1) = h;
    % set(h,'XData',my_xdata,'YData',my_ydata,'ZData',my_zdata);
    
    
    
    
    % do it again!! for NBS-processed graphs:
    % we have results from NBS (see above):
    hf=figure;hf.Color='w';
    bigfs(end+1) = hf;
    
    
    hf.Position=[100 100 1000 500];
    
    % tmp=threshold_proportional(raw_mtx,dns);
    % thresholding trans_prob mtx
    
    G=digraph(nbs_results.(part));
    
    
    
    
    tmpnbs=nbs_results.(part);
    % make the excel (ugh, excel..)
    for ii=1:NSTATES
        s=sprintf('ToState%d=tmpnbs(:,%d)',ii,ii);
        eval(s);
    end
    
    s=sprintf('x = {');
    for ii=1:NSTATES
        s=[s sprintf('\''FromState%d\'';',ii)];
    end
    s(end)=[];
    s=[s sprintf('}')];
    eval(s);
    
    s=sprintf('myTable=table(x,');
    for ii=1:NSTATES
        s=[s sprintf('ToState%d,',ii)];
    end
    s(end)=[];
    s=[s sprintf(')')];
    eval(s);
    

    % savefig_title...
    savefig_iii = sprintf('../figures/Fig5_NBSFigure_nbsresult_%s_data.xls',part);
    writetable(myTable, savefig_iii);
    
    
    
    w=G.Edges.Weight; w=w/max(w); w=w*1; % scaling the edges
    % you might need to vary 5 depending on your figure. Increase if edges are
    % thin and decrease if they appear very thick.
    if numel(w)>0
        LWidths=w;
        h=plot(G,'LineWidth',LWidths,'Layout','layered');
    else
        h=plot(G);
    end
    % h=plot(G,'LineWidth',LWidths,'Layout','layered','XData',my_xdata,'YData',my_ydata,'ZData',my_zdata);
    
    
    title(nbs_titles.(part),'Interpreter','none');
    axis off; h.ArrowSize=10; % changes the arrow size
    % again might need to vary arrow size
    whos t
    h.NodeColor=cmat(WHICHSTATES,:); % set node color to green
    h.NodeLabel=NodeLabels;
    
    % setting nodes 1,2 &5 to cyan color
    % this might hep to distinguish between sattes that have
    % significantly different Fo values between groups
    % highlight(h,[1 2 5],'NodeColor','cy');
    
    % h.MarkerSize=4*(1:12); % changing node size
    h.MarkerSize=0.25*round(my_wts)+4; % changing node size
    % since our trans_prob is 12x12, we have 12 nodes and I'm just
    % making them increase in size from 1 to 12.
    % in practice, you can chnage this according toavg state FO values;
    % so that states occupied most appear larger
    
    % for customising the plot further, look at the properties of object 'h'
   %  else
     %    h=plot(G);
   %  end
    
    collected_handles(end+1) = h;
    % set(h,'XData',my_xdata,'YData',my_ydata,'ZData',my_zdata);

    % keyboard;
    
    
    
    
end


% % die a 1000 deaths for this ugliness which is Matlab-inherent. -- set each
% % to the first one...
% properties_hs={'XData','YData','ZData'};
% for j=1:numel(properties_hs)
%     this_property=properties_hs{j};
%     for i=2:numel(collected_handles)
%         set(collected_handles(i),this_property,get(collected_handles(1),this_property));
%     end
% end
% % layout(h,'layered')
% %h2.XData=h.XData;h2.YData=h.YData;h2.ZData=h.ZData;










%% 
% combine all the 4 figures into ONE figure:

all_figs = [matrix_fhs bigfs];
all_figs = all_figs([1 4 2 5 3 6]);
all_figs=all_figs([1 5 2 3 4 6]);
bigfs=all_figs;

fcols=3;
frows=2;
inums=[2 1];

spacingx=0.02;
spacingy=0.02;

% borders of our image:
bigspacingxl = 0.05;
bigspacingyl = 0.05;

bigspacingxu = 0.05;
bigspacingyu = 0.05;

factorx = (1-bigspacingxl-bigspacingxu-(fcols-1)*spacingx)/fcols;
factory = (1-bigspacingyl-bigspacingyu-(frows-1)*spacingy)/frows;

newbigf=figure('color','w');
set(newbigf,'position',[50 50 1440 900]);
% combine into one UBER figure -- of 2-by-2
% inums=[2 1];  % so ugly, thanks matlab...
for ii=1:numel(inums)
    i=inums(ii);
    for j=1:fcols
        
        this_fig=bigfs((i-1)*fcols+j);
        to_copy=get(this_fig,'children');
        
        for k=1:numel(to_copy)
            
            old_pos = get(to_copy(k),'position');
            
            
            old_xp = old_pos(1);
            old_yp = old_pos(2);
            old_xs = old_pos(3);
            old_ys = old_pos(4);
            
            
            x_offset = (factorx+spacingx)*(j-1);
            
            % x_offset = x_offset + factorx + spacingx;
            
            
            new_xp = old_xp * factorx + x_offset + bigspacingxl;
            new_yp = old_yp * factory + (ii-1) * (factory + spacingy) + bigspacingyl;
            new_xs = old_xs * factorx;
            new_ys = old_ys * factory;
            
            new_position = [new_xp new_yp new_xs new_ys];
            
            
            newobj=copyobj(to_copy(k),newbigf);
            set(newobj,'position',new_position);
            
            
            
        end
    end
end


% add in our colorbar... in the middle, preferably
all_chs=get(newbigf,'children');
p1=get(all_chs(3),'Position');
p4=get(all_chs(6),'position');
ypos=p4(2)+p4(4)+0.05;
xpos=p1(1);
xsize = p1(3);
ysize = 0.05;
cb_ax = axes('parent',newbigf,'position',[xpos ypos xsize ysize]);

cb=colorbar('South');
set(cb_ax,'visible','off');
set(cb_ax,'clim',[0 0.09]);


% close(bigfs);


% and here should come -- the annotations, in very similar fashion as above
% - but then different.

analysis = 'mov1a-1b';
set(newbigf,'paperunits','centimeters');
set(newbigf,'papersize',[25 20],'paperposition',[0 0 25 20]);
% 
fsource = ['../figures/Fig5_Dunamics_' run '.pdf'];
% ftarget_file = fsource;
saveas(newbigf,fsource);
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

for i=1:6
    
    if i==1 || i==2
        papersize=6;
    else
        papersize=6;
    end
   
    set(i,'paperunits','centimeters')
    set(i,'papersize',papersize * [1 1]);
    set(i,'paperposition', papersize * [0 0 1 1]);
    saveas(i,[sprintf('tmpfig_sesA_%d_%s.jpg',i, run)]);
    saveas(i,[sprintf('tmpfig_sesA_%d_%s.pdf',i, run)]);
end








%
%
% clear;clc; close all; tic;
%
% con_avg_prb=rand(12,12); % the trans. prob. mtx that needs to be visulaized
% % avearge transitions in controls
%
% dns=0.1;% density.
% % prcntage of connections that needs to be seen
% % 0.1 = top 10% connections are displayed
%
% hf=figure;hf.Color='w';
% hf.Position=[100 100 1000 500];
%
% tmp=threshold_proportional(con_avg_prb,dns);
% % thresholding trans_prob mtx
% G=digraph(tmp);
% w=G.Edges.Weight; w=w/max(w); w=w*5; % scaling the edges
% % you might need to vary 5 depending on your figure. Increase if edges are
% % thin and decrease if they appear very thick.
% LWidths=w;
% h=plot(G,'LineWidth',LWidths,'Layout','layered');
% title('Avg Trans_Prob - Controls','Interpreter','none');
% axis off; h.ArrowSize=17; % changes the arrow size
% % again might need to vary arrow size
%
% h.NodeColor='g'; % set node color to green
%
% % setting nodes 1,2 &5 to cyan color
% % this might hep to distinguish between sattes that have
% % significantly different Fo values between groups
% highlight(h,[1 2 5],'NodeColor','cy');
%
% h.MarkerSize=4*(1:12); % changing node size
% % since our trans_prob is 12x12, we have 12 nodes and I'm just
% % making them increase in size from 1 to 12.
% % in practice, you can chnage this according toavg state FO values;
% % so that states occupied most appear larger
%
% % for customising the plot further, look at the properties of object 'h'
%
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % now doing the same for patients
% pat_avg_prb=rand(12,12); % avearge transitions in patients
%
% dns=0.1;% density. % of connections that needs to be seen 0.1 = top 10% connections are displayed
% hf=figure;hf.Color='w';
% hf.Position=[100 100 1000 500];
%
% tmp=threshold_proportional(pat_avg_prb,dns);
% G=digraph(tmp);
% w=G.Edges.Weight; w=w/max(w); w=w*5; % scaling the edges
% LWidths=w;
% h2=plot(G,'LineWidth',LWidths,'Layout','layered');
% title('Avg Trans_Prob - Patients','Interpreter','none');
% axis off; h2.ArrowSize=17; % changes the arrow size
%
% h2.NodeColor='g'; % set node color to green
% highlight(h2,[1 2 5],'NodeColor','cy'); % setting nodes 1,2 &5 to cyan color
%
% h2.MarkerSize=4*(1:12); % changing node size
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % This is the key part resettting the x,y,z co-ordinates of second plot.
% % Just copying those from the fist plot, so that state 1 in both plots
% % appear at the same palce, for easy comaprison.
% h2.XData=h.XData;h2.YData=h.YData;h2.ZData=h.ZData;
% toc;