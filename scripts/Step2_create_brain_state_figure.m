base_cwd = pwd;
addpath(base_cwd);
% we need spm12 to read in images
addpath('/home/arasha/Documents/toolboxes/spm12/'); 


load bw.mat

addpath(genpath(pwd));
CLOSE_FIGS=0;

NSTATES=10;

VISIBILITY='on';

preprocessings={'normal','aroma','aroma-gsr'};



% descriptors of the network - only valid for NSTATES = 10, and permutation
% 13 (in our case)
network_descriptions = {...
    sprintf('High EXEC, SENS, LAN,\nlow SAL (ANS)');...
    sprintf('High DMN, SAL, low EXEC');...
    sprintf('High DMN (dDMN),\nhigh SENS (AUD), high LAN');...
    sprintf('High visual (hVIS)');...
    sprintf('Mean brain networks activity');...
    sprintf('High DMN and LAN,\nlow SENS (VSN)');...
    sprintf('High VIS (pVIS),\ngeneral networks suppression');...
    sprintf('High DMN, SAL (ASN),\nhigh SENS (AUD), low LAN');...
    sprintf('Low brain networks activity');...
    sprintf('High SAL and EXEC (BGN),\nlow SENS (pVIS and hVIS)');...
};...

network_descriptions = {'','','','','','','','','',''};

% these colors are for coding the networks. They should be maximally
% discernable. 
% See: https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
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

% network_descriptions = {'','','','','','','','','','','',''};
close all;


% types of analyses
analyses = {...
    {'mov1a'},'mov1a'; ... % only mov1
    {'mov1b'},'mov1b'; ... % only mov2 day2
    {'mov1a','mov1b'},'mov1a-1b'; ... % mov1 day and and day 2
    {'resta'},'resta'; ...  % only resta
    {'restb'},'restb'; ...  % only rest day 2
    {'resta','restb'},'resta-b'; ...  % rest a and rest day 2
    {'mov1a','mov1b','resta','restb'},'all'; ...  % all movie and rest day 1 aand day 2 but not repeats
    };


fbgcolor = [1 1 1]; % [0.95 0.95 0.95 ];
BOUTLINE = 1; % brain-outline


for i_summarymeasures=1:15 % set to [1:15] if you wish to see everything
    SUMMARY_MEASURES_FILE=i_summarymeasures;
    
    for i_preprocessing=2 % set to [1:3] if you wish to make all figures
        preprocessing_dir=preprocessings{i_preprocessing};
        
        for i_analyses=7  % set to [1:7] if you wish to make this figure for all analyses
            
            
            % change-dir into the correct analysis directory
            analysis_dir=analyses{i_analyses, 2};
            
            cd([base_cwd '/../results_' num2str(NSTATES)]);
            cd(preprocessing_dir)
            cd(analysis_dir)
            
            
            

            cmap_activities=colormap(bw);
            
            
            % we have 1 file for each of the spatial masks
            roi_base = [base_cwd '/../data/rois/'];
            
            roi_files={};
            roi_vol={};
            for i=1:14
                roi_files{i} = [roi_base 's' sprintf('%.2d',i) '-to-2009c-in-2009c-dims.nii'];
                roi_vol{i} = spm_vol(roi_files{i});
                roi_data{i} = spm_read_vols(roi_vol{i});
                
                
            end
            
            % the standard t1 scan
            t1_file = [roi_base '2009c-brain.nii'];
            t1_vol = spm_vol(t1_file);
            t1_data = spm_read_vols(t1_vol); %-- it's a 193x229x193
            
            % a slightly dilated t1 scan (to make outline of the brain)
            t1_outline_file = [roi_base '2009c-brain-binned-dil1.nii'];
            t1_outline_vol = spm_vol(t1_outline_file);
            t1_outline_data = spm_read_vols(t1_outline_vol);
            
            
            % check Summary_measures_1
            load(['HMMrun_rep_' num2str(SUMMARY_MEASURES_FILE) '.mat']);
            load(['Summary_measures_rep_' num2str(SUMMARY_MEASURES_FILE) '.mat']);
            NSTATES=10;
            
            % get the coeffs of all the states this yields a roi-by-state
            % matrix. each col is a representation of which roi's/networks
            % are represented in that state. we go through states later and
            % make nice T1 figs.
            allw=[];
            for i=1:NSTATES
                allw=[allw hmm.state(i).W.Mu_W']; 
            end
            
            
            %coronal_slice=100;
            sl=100;
            thr_fig=0.3;    % theshold to binarize the brainmask (has 
                            % continuous values for the T1 scan
                            % -- the brain masks for the functional data
                            % only have 0/1 values.
            clims=[-0.5 0.5]; % lim for normalized projection of the brain 
                              % state  onto the brain network
            
            
            fh=[];  % to store figure handles. We basically make 10 
                    % separate figures (one for each state). In each figure
                    % there are 2 + 14 (T1 + Outline + 14 brain network)
                    % images drawn over each other with transparency set
                    % properly. Then one big figure is made, and for each
                    % of the 10 figures, the contents are copied into the
                    % big figure.
            for state=1:NSTATES
                
                
                fh(state)=figure('position',[400,400,1000,200],'color',fbgcolor,'visible',VISIBILITY);
                set(fh(state),'units','normalized')
                
                
                if BOUTLINE
                    ah0=axes('parent',fh(state),'position',[0.20 0 0.23 1]);
                    t1_outline=0.5 * squeeze(t1_outline_data(sl,25:220,02:173))';
                    im1=imagesc(ah0,1-t1_outline);
                    set(ah0,'clim',[0 1]);
                    % set(im1,'alphadata',t1_outline>0);
                    axis xy
                    colormap(ah0,'gray')
                    set(ah0,'visible','off');
                end

                
                ah1=axes('parent',fh(state),'position',[0.20 0 0.23 1]);
                t1data=squeeze(t1_data(sl,25:220,02:173))';
                im1=imagesc(ah1,t1data);
                set(im1,'alphadata',t1data>0);
                axis xy
                colormap(ah1,'bone')
                set(ah1,'visible','off');
                
                
                ah_state=axes('parent',fh(state),'position',[0.45 0 0.53 0.1]);
                set(ah_state,'xlim',[0 14],'ylim',[0 1])
                colormap(ah_state,cmap_activities);
                set(ah_state,'visible','off');
                set(ah_state,'clim',clims);
                
                
                % to alter the sequence of the networks so things are
                % consistent with Kottaram et al. Paper.
                translate_mat = [4 9 8 1 11 14 6 5 2 10 7 12 3 13];
                phs=[];
                for i=1:14
                    
                    this_w = allw(i, state);
                    ah2=axes('parent',fh(state),'position',[0.20 0 0.23 1]);
                    roidata=squeeze(roi_data{i}(sl,25:220,02:173))';
                    im2=imagesc(ah2,(roidata>thr_fig)*this_w);
                    axis xy
                    colormap(ah2,cmap_activities);
                    set(ah2,'clim',clims);
                    set(im2,'alphadata',roidata>thr_fig);
                    set(ah2,'visible','off');
                    colormap(ah2);
                    % colormap(im2,cmap_activities);
                    
                    
                    patch_coords_x = [0 1 1 0] + (translate_mat(i)-1);
                    patch_coords_y = [0 0 1 1];
                    phs(translate_mat(i)) = patch(ah_state,patch_coords_x,patch_coords_y,this_w);
                    
                end
                
                % some magic to make the title of the brain network more
                % legible.
                pdiff=0.06;
                set(phs(3),'XData', get(phs(3),'XData') - [0 pdiff pdiff 0]')
                set(phs(4),'XData', get(phs(4),'XData') + [pdiff 0 0 pdiff]')
                set(phs(5),'XData', get(phs(5),'XData') - [0 pdiff pdiff 0]')
                set(phs(6),'XData', get(phs(6),'XData') + [pdiff 0 0 pdiff]')
                set(phs(8),'XData', get(phs(8),'XData') - [0 pdiff pdiff 0]')
                set(phs(9),'XData', get(phs(9),'XData') + [pdiff 0 0 pdiff]')
                set(phs(13),'XData', get(phs(13),'XData') - [0 pdiff pdiff 0]')
                set(phs(14),'XData', get(phs(14),'XData') + [pdiff 0 0 pdiff]')
                
                
                % brain outline?
                if BOUTLINE
                    % brain outline - to make the networks more pronounced if
                    % they're white and on the edge of the brain.
                    ah1axoutl=axes('parent',fh(state),'position',[0.0 0 0.20 1]);
                    t1data=0.5 * squeeze(t1_outline_data(:,:,80))';
                    im1=imagesc(ah1axoutl,1-t1data);
                    set(ah1axoutl,'clim',[0 1]);
                    % set(im1,'alphadata',t1data>0);
                    axis xy
                    colormap(ah1axoutl,'bone')
                    set(ah1axoutl,'visible','off');
                end
                
                % make an axial section?
                ah1ax=axes('parent',fh(state),'position',[0.0 0 0.20 1]);
                t1data=squeeze(t1_data(:,:,80))';
                im1=imagesc(ah1ax,t1data);
                set(im1,'alphadata',t1data>0);
                axis xy
                colormap(ah1ax,'bone')
                set(ah1ax,'visible','off');
                
                
                % do for each ROI the following:
                translate_mat = [4 9 8 1 11 14 6 5 2 10 7 12 3 13];
                phs=[];
                for i=1:14
                    
                    this_w = allw(i, state);
                    ah2ax=axes('parent',fh(state),'position',[0.0 0 0.20 1]);
                    roidata=squeeze(roi_data{i}(:,:,80))';
                    im2ax=imagesc(ah2ax,(roidata>thr_fig)*this_w);
                    axis xy
                    colormap(ah2ax,cmap_activities);
                    set(ah2ax,'clim',clims);
                    set(im2ax,'alphadata',roidata>thr_fig);
                    set(ah2ax,'visible','off');
                    colormap(ah2ax);
                    % colormap(im2,cmap_activities);
                    
                end

                
                ah_text = axes('parent',fh(state),'position',[0.45 0.2 0.53 0.7]);
                set(ah_text,'visible','off');
                ph=patch(ah_text,[0 1 1 0],[0 0 1 1],[0.95 0.95 0.95]);
                set(ph,'edgecolor',[0.95 0.95 0.95]);
                thxl=get(gca,'xlim');
                thyl=get(gca,'ylim');
                th3=text(thxl(1)-diff(get(gca,'xlim'))*0.0005,thyl(2),sprintf(' State %d',state),'fontweight','bold','verticalalignment','top');
                th3=text(thxl(1)+diff(get(gca,'xlim'))*0.0005,thyl(2),sprintf(' State %d',state),'fontweight','bold','verticalalignment','top');
                th3=text(thxl(1),thyl(2)-diff(get(gca,'ylim'))*0.0005,sprintf(' State %d',state),'fontweight','bold','verticalalignment','top');
                th3=text(thxl(1),thyl(2)+diff(get(gca,'xlim'))*0.0005,sprintf(' State %d',state),'fontweight','bold','verticalalignment','top');
                th3=text(thxl(1),thyl(2)+diff(get(gca,'xlim'))*0.0005,sprintf(' State %d',state),'fontweight','bold','verticalalignment','top','color',cmat(state,:)/255);
                %th=text(thxl(1),thyl(2),sprintf(' Network %d',state),'fontweight','bold','Color',cmat(state,:)/255);
                
                if NSTATES == 10
                    th2=text(thxl(1),thyl(2),sprintf('\n\n%s',network_descriptions{state}),'verticalalignment','top');
                    set(th2,'interpreter','none');
                end
                
                
            end
            
            
            
            % the magic -- make the Big Plot - subdivide into 5-by-2 figs
            % but we need an fh array of fig handles.
            bigf=figure('color',fbgcolor,'visible',VISIBILITY);
            set(bigf,'position',[50 50 800 800])
            fcols=2;
            frows=5;
            offx=0; %% how much space at bottom of figure?
            offy=0.2;
            
            spacing=0.001;
            factorx = (1-offx-(fcols-1)*spacing)/fcols;
            factory = (1-offy-(frows-1)*spacing)/frows;
            
            for i=1:NSTATES %:-1:1
                [fj, fi] = ind2sub([frows fcols],i);
                fj=frows+1-fj;
                
                objs = findobj(fh(i),'type','axes'); 
                
                for i_o = numel(objs):-1:1
                    o=objs(i_o);
                    old_position = get(o, 'position');
                    old_xp = old_position(1);
                    old_yp = old_position(2);
                    old_xs = old_position(3);
                    old_ys = old_position(4);
                    
                    % calculate new position according to my notes:
                    new_xp = old_xp * factorx + (fi-1) * (factorx + spacing) + offx;
                    new_yp = old_yp * factory + (fj-1) * (factory + spacing) + offy;
                    new_xs = old_xs * factorx;
                    new_ys = old_ys * factory;
                    
                    new_position = [new_xp new_yp new_xs new_ys];
                    
                    % then copy that obj into the new figure:
                    newo = copyobj(o, bigf);
                    set(newo,'position',new_position);
                end
            end
            
            
            cb_ax = axes('parent',bigf);
            set(cb_ax,'position',[0.2 0.09 0.6 0.06]);
            cb=colorbar('North');
            set(cb,'Limits',clims);
            set(cb_ax,'clim',clims);
            colormap(cb_ax,cmap_activities);
            set(cb_ax,'visible','off');
            
            annot_ax = axes('parent',bigf);
            set(annot_ax,'position',[0.15 0.03 0.7 0.05],'visible','off');
            set(annot_ax,'xlim',[0 14]);
            set(annot_ax,'ylim',[0 2]);
            
            DESCRIPTORS = {'dDMN','PRE','vDMN','ASN','PSN','LECN','RECN','BGN','AUD','pVIS','hVIS','SMN','VSN','LAN'};
            phs=[];
            for i=1:14
                patch_coords_x = [0 1 1 0] + (i-1);
                patch_coords_y = [0 0 1 1]+1;
                ph=patch(annot_ax,patch_coords_x,patch_coords_y,0);
                phs(end+1) = ph;
                set(ph,'facealpha',0);
                pdiff=0.06;
                switch i
                    case {3, 5, 8, 13}
                        th=text(i-0.5-pdiff/2,1.5,DESCRIPTORS{i});
                        % disp('doingit');
                    case {4, 6, 9, 14}
                        th=text(i-0.5+pdiff/2,1.5,DESCRIPTORS{i});
                        % disp('doingit2');
                    otherwise
                        th=text(i-0.5-pdiff/2,1.5,DESCRIPTORS{i});
                        % disp('doingit3');
                end
                
                
                
                set(th,'fontsize',7,'horizontalalignment','center');
            end
            
            
            % DMN line
            annots_labels={'DMN','SAL','EXEC','SENS'};
            annots_begins = [0, 3, 5, 8];
            annts_ends = [3, 5, 8, 13];
            
            for iannot=1:4
                b=annots_begins(iannot);
                e=annts_ends(iannot);
                t=annots_labels{iannot};
                
                line([b, b, e, e] ,[0.75, 0.5, 0.5, 0.75],'color','k');

                th=text(mean([b, e]), 0.4,t);
                set(th,'horizontalalignment','center','verticalalignment','top','fontsize',7,'color','k');
                % line([b b],[1 2],'linewidth',2,'color','k');
                % line([e e],[1 2],'linewidth',2.5,'color','k');
                % line([e e],[1 2.05],'linewidth',0.5,'color','w');
                % line([b e e b b],[1 1 2 2 1],'linewidth',2,'color','k');
            end
            
            set(phs(3),'XData', get(phs(3),'XData') - [0 pdiff pdiff 0]')
            set(phs(4),'XData', get(phs(4),'XData') + [pdiff 0 0 pdiff]')
            set(phs(5),'XData', get(phs(5),'XData') - [0 pdiff pdiff 0]')
            set(phs(6),'XData', get(phs(6),'XData') + [pdiff 0 0 pdiff]')
            set(phs(8),'XData', get(phs(8),'XData') - [0 pdiff pdiff 0]')
            set(phs(9),'XData', get(phs(9),'XData') + [pdiff 0 0 pdiff]')
            set(phs(13),'XData', get(phs(13),'XData') - [0 pdiff pdiff 0]')
            set(phs(14),'XData', get(phs(14),'XData') + [pdiff 0 0 pdiff]')
            
            
            

            set(cb_ax,'visible','off');
            
            if CLOSE_FIGS
                close(fh);
            end
            cb.Position=cb.Position + [0 0 0 0.009];
            
            cb.Position=cb.Position - [0 0.01 0 0];
            
            
            a_spacingx = 0.1;
            a_spacingy = 0.9;
            
            a1_x = [cb.Position(1) + a_spacingx * cb.Position(3), cb.Position(1)];
            a1_y = cb.Position(2) + (1+a_spacingy) * cb.Position(4) * [1 1];
            
            a2_x = [cb.Position(1) + (1-a_spacingx) * cb.Position(3), cb.Position(1) + cb.Position(3)];
            a2_y = a1_y;
            
            a1_h = annotation('textarrow',a1_x, a1_y,'String','  Low');
            a2_h = annotation('textarrow',a2_x, a2_y,'String','High  ');
            
            
            anno_ax = axes('parent',bigf,'position',[0 0 1 1],'visible','off');
            a3_h = text(cb.Position(1) + 0.5*cb.Position(3),a2_y(1),'Mean','HorizontalAlignment','Center');
            set(a3_h,'Parent',anno_ax);
            
            
            set(cb,'Ticks', [-0.5 0 0.5]);
            set(cb,'TickLabelInterpreter','None');
            set(cb,'TickLabels',{'-50%','0','+50%'});
            
            set(cb,'FontSize',10);
            set(cb,'FontWeight','bold');
            
            
            fh_overview_states=bigf;
            
            
            set(fh_overview_states,'paperunits','centimeters');
            set(fh_overview_states,'papersize',1.25 * [22 16]);
            set(fh_overview_states,'paperposition',1.25 * [0 0 22 16]);
            
            
            saveas(fh_overview_states,[preprocessing_dir '-' analysis_dir '-' num2str(SUMMARY_MEASURES_FILE) '.jpg']);
            
            if CLOSE_FIGS
                close(fh_overview_states);
            end
            
            fprintf('Done making Fig1 for %s, %s, file: %d\n',preprocessing_dir, analysis_dir, SUMMARY_MEASURES_FILE);
            
            % the following can be uncommented toautomatize that one of the
            % figures is copied automatically to a 'figures' directory.
            % fsource=[preprocessing_dir '-' analysis_dir '-' num2str(SUMMARY_MEASURES_FILE) '.jpg'];
            
            % ftarget_file = 'f1_state_10.jpg';

            % ftarget1=[base_cwd '/../figures/' ftarget_file];
            % copyfile(fsource,ftarget1);

            
            
            cd(base_cwd);
            
        end
    end
    
    
    
    
end


