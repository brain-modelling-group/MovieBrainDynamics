% clear all; close all; clc;


ANALYSIS = 'all';



base_cwd = pwd; % run the script from the directory in which it is located.
addpath(base_cwd);
addpath(genpath(base_cwd));
ANALYSIS='all';
load(['valid_inferences_' ANALYSIS '.mat']);
% this_inference = 1241;
scriptpath=pwd;

SCANS_NULLED=5;
CLOSE_FIGS=0;
NSTATES=10;

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
    {'mov1b'},'mov1b'; ... % only mov2 day2
    {'mov1a','mov1b'},'mov1a-1b'; ... % mov1 day and and day 2
    {'resta'},'resta'; ...  % only resta
    {'restb'},'restb'; ...  % only rest day 2
    {'resta','restb'},'resta-b'; ...  % rest a and rest day 2
    {'mov1a','mov1b','resta','restb'},'all'; ...  % all movie and rest day 1 aand day 2 but not repeats
    };


% apply the nice colormap for maximum contrast:
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

% the parula colormap is used for the binary masks (i.e. story annotations)
parulamap=colormap('parula');

cmat=cmat(1:NSTATES,:);
cmat=cmat./255;


% container for all data -- PER SUBJECT
sd=struct(); % subject data

my_cond_names = {'mova','movb','resta','restb'};
new_statemat = struct(); new_statemat(:)=[]; % this is for counting/keeping track of the states.



preprocessings={'normal','aroma','aroma-gsr'};


for i_summarymeasures=valid_inferences %valid_inferences(1:15)
    SUMMARY_MEASURES_FILE=i_summarymeasures;
    
    for i_preprocessing=2 %1:numel(preprocessings)
        preprocessing_dir=preprocessings{i_preprocessing};
        
        
        % figure out which of the analyses we need to grab depending on
        % what the user specified in ANALYSIS:
        i_analyses = find(strcmp(analyses(:, 2),ANALYSIS));
        
        
        
        
        all_anno_axes=[];
        fh=figure('color','w','visible','off');
        set(fh,'position',[50         120        1700         800]);
        if ~CLOSE_FIGS
            set(fh,'visible','on');
        end
        
        scan_names = analyses{i_analyses, 1};
        
        analysis = analyses{i_analyses, 2};
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
        
        % keyboard;
        
        
        analysis_dir=analyses{i_analyses, 2};
        disp(analysis_dir);
        disp(preprocessing_dir);
        
        cd([base_cwd '/../results_' num2str(NSTATES)]);
        cd(preprocessing_dir)
        cd(analysis_dir)
        
        
        load(['HMMrun_rep_' num2str(SUMMARY_MEASURES_FILE) '.mat']); % this will get our vpath.
        load(['Summary_measures_rep_' num2str(SUMMARY_MEASURES_FILE) '.mat']); % this will get our vpath.
        
        switch analysis_dir
            
            case 'mov1a'
                nvols=[535];
            case 'mov1b'
                nvols=[535];
            case 'mov1a-1b'
                nvols=[535 535];
            case 'resta'
                nvols=[220];
            case 'restb'
                nvols=[220];
            case 'resta-b'
                nvols=[220 220];
            case 'all'
                nvols=[535 535 220 220];
        end
        
        
        nvols=nvols-SCANS_NULLED;
        % edit pr
        
        
        
        nsubs=21;
        subs_to_use=1:21;
        for i=1:numel(scan_names)
            subs_to_use = intersect(subs_to_use, slist.(scan_names{i}));
            nsubs=numel(subs_to_use);
        end
        
        
        
        nscans=sum(nvols);
        nsub = numel(vpath)/nscans;
        
        
        % the state path is one big vector of everything - here we
        % subdivide by subject/session again. The b and e are running
        % indices that change number upon each iteration, so that b:e
        % 'slices' the correct partition of vpath.
        mat={};
        e_vpath=0;
        b_vpath=1;
        for invols=1:numel(nvols)
            nvol=nvols(invols);
            e_vpath = e_vpath+ nsub * nvol;
            tot_el=numel(b_vpath:e_vpath);
            mat{end+1} = reshape(vpath(b_vpath:e_vpath),tot_el/nsub,nsub);
            b_vpath=e_vpath+1;
        end
        
        
        %
        %
        % make the figure.
        %
        %
        
        
        fcols=numel(mat);
        frows=nsub;
        
        spacingx=0.012;
        spacingy=0.006;
        bigspacingxl = 0.10;
        bigspacingyl = 0.10;
        
        bigspacingxu = 0.20;
        bigspacingyu = 0.50;
        
        
        factorx = (1-bigspacingxl-bigspacingxu-(fcols-1)*spacingx)/fcols;
        factory = (1-bigspacingyl-bigspacingyu-(frows-1)*spacingy)/frows;
        ah=[];
        im=[];
        
        
        % so now -- some magic to scale the factors, so as to have it
        % even nicer.
        
        extra_factorx=[1 1 1 1 1];
        
        total_scans=[];
        for i_gr=1:numel(mat)
            total_scans(end+1) = size(mat{i_gr}, 1);
        end
        for i_gr=1:numel(mat)
            extra_factorx(i_gr) = 1/(1/numel(mat) / (total_scans(i_gr)/sum(total_scans)));
        end
        
        
        
        
        chosen_subs=subs_to_use;
        
        
        
        % in which sequence to present figures:
        if numel(scan_names) == 4
            new_order=[3 1 4 2];
        else
            new_order=[1 2];
        end
        
        for i_gr2=1:numel(new_order)
            i_gr = new_order(i_gr2);
            
            x_offset = 0;
            if i_gr2 > 1
                for ii_gr=1:(i_gr2-1)
                    x_offset = x_offset + factorx*extra_factorx(new_order(ii_gr)) + spacingx;
                end
            end
            
            for i=1:nsub % *numel(nvols)
                
                
                [fj, fi] = ind2sub([frows fcols],i+numel(subs_to_use)*(i_gr-1));
                
                
                this_sub = i; % chosen_subs(nsub-i+1);
                
                % v=mat(:,nsub-(i-1))';
                v=mat{i_gr}(:,this_sub)';
                
                
                % save v into big data struct:
                sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(subs_to_use(i))]).(regexprep(analyses{i_analyses,1}{i_gr},'-','_')).hmmdata=v;
                
                
                % convert state indices into 0-1 matrix
                % save into big data struct:
                statemat=zeros(numel(v),NSTATES);
                % statemat(:)=NaN;
                for itstate=1:NSTATES
                    statemat(v==itstate,itstate) = 1;
                    if sum(v==itstate)==0
                        statemat(:,itstate)=NaN;
                    end
                end
                sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(subs_to_use(i))]).(regexprep(analyses{i_analyses,1}{i_gr},'-','_')).hmmmat=statemat;
                
                
                
                old_xp = 0; % old_position(1);
                old_yp = 0; % old_position(2);
                old_xs = 1; % old_position(3);
                old_ys = 1; % old_position(4);
                
                new_xp = old_xp * factorx*extra_factorx(i_gr) + x_offset + bigspacingxl;
                new_yp = 1-(old_yp * factory + (fj) * (factory + spacingy) + bigspacingyl);
                new_xs = old_xs * factorx*extra_factorx(i_gr);
                new_ys = old_ys * factory;
                
                new_position = [new_xp new_yp new_xs new_ys];
                
                ah(end+1) = axes('parent',fh,'position',new_position);
                
                
                
                im(i) = imagesc(v);
                set(ah(end),'visible','off');
                set(ah(end),'colormap',cmat);
                set(ah(end),'clim',[0.5 NSTATES+0.5]);
                % set(ah(end),'colormap',cmat);
                
                if i==1
                    % keyboard;
                    yl=get(ah(end),'ylim');
                    xl=get(ah(end),'xlim');
                    tha=text(mean(xl),yl(1),sprintf('%s', regexprep(regexprep(analyses{i_analyses,1}{i_gr},'-','_'),'1','')));
                    set(tha,'fontweight','bold');
                    set(tha,'horizontalalignment','center','verticalalignment','bottom');
                    set(tha,'parent',ah(end));
                end
                
                if i_gr2==1
                    th=text(ah(end),-10,0,sprintf('P%d - %d',i,''));
                    set(th,'horizontalalignment','right','units','normalized');
                    set(th,'position',[-0.001 0.5]);
                end
                
                
                
                
            end
            
            ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-factory+0.005 new_xs 0.010]);
            %set(ah(end),'ytick',[],'yticklabel',[],'ycolor','w');
            set(ah(end),'xlim',[0 numel(v)*2.2]);
            time_v = 2.2*(1:numel(v));
            ph=plot(time_v, ones(length(time_v)));
            set(ph,'visible','on','color','w');
            set(ah(end),'yticklabel',{},'ytick',[]);
            set(ah(end),'xlim',[0 numel(v)*2.2]);
            xlh=xlabel('time (s)');
            % set(xlh,'horizontalalignment','right');
            box off
            
            
            
            
            
            %
            %
            %      EXTRA FIG -- STATES (the most common ones)
            %
            %
            %
            
            
            
            time_lag=9;
            STHR=7;
            OFFSET=7;
            
            
            m=mat{i_gr}';
            
            
            outs=[];
            outc=[];
            for im=1:size(m,2)
                
                b=im-floor(time_lag/2);
                e=im+floor(time_lag/2);
                
                
                if b < 1 || e > size(m,2)
                    out(im)=NaN;
                    
                else
                    
                    
                    all_vals=m(:,b:e);
                    
                    uv=unique(all_vals);
                    ct=[];
                    for iuv=1:numel(uv)
                        ct(end+1)=sum(sum(all_vals==uv(iuv),2)>0);
                    end
                    
                    um=sortrows([uv ct'],2,'descend');
                    
                    outs(end+1,:)=um(1:2,1);
                    outc(end+1,:)=um(1:2,2);
                    
                end
            end
            
            selection=outc(:,1)>STHR;
            
            vout=outs(:,1);
            
            
            vout(selection==0)=NaN;
            
            
            ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-0.05-0.05-0.03+0.02 new_xs 0.05-0.01]);
            
            
            
            inpatch=0;
            vout = [nan(floor(time_lag/2),1);vout];
            for i=2:numel(vout)
                
                
                if ~inpatch
                    line_px=[(i-0.5)*2.2];
                    line_py=[OFFSET];
                end
                
                if ~isnan(vout(i))
                    
                    inpatch=1;
                    
                    if vout(i) ~= vout(i-1) && ~isnan(vout(i-1))
                        line_px(end+1) = (i-0.5)*2.2;
                        line_py(end+1) = OFFSET;
                        line_px(end+1) = line_px(1);
                        line_py(end+1) = OFFSET;
                        draw_it;
                        % keyboard;
                        
                        line_px=[(i+[-0.5])*2.2];
                        
                        line_py=[OFFSET];
                    end
                    
                    
                    line_px(end+1)=(i+[-0.5])*2.2;
                    line_px(end+1)=(i+[ 0.5])*2.2;
                    
                    line_py(end+1)=outc(i-floor(time_lag/2),1)*[1 ];
                    line_py(end+1)=outc(i-floor(time_lag/2),1)*[1 ];
                    
                    val=vout(i);
                    confidence=outc(i-floor(time_lag/2),1);
                    
                    % ph_back=patch(2.2*(i+[-0.5 0.5 0.5 -0.5]), [0 0 1 1]*confidence,[1 1 1]);
                    ph=patch(2.2*(i+[-0.5 0.5 0.5 -0.5]), [0 0 1 1]*confidence,cmat(val,:));
                    set(ph,'linestyle','none');
                    set(gca,'xlim',[0 size(m,2)*2.2]);
                    set(gca,'ylim',[OFFSET-0.2 size(m,1)]);
                    set(gca,'xtick',[],'xticklabel',[]);
                    if i_gr2 > 1
                        set(gca,'ytick',[],'yticklabel',[]);
                    end
                    
                    if i_gr2 == 1
                        set(gca,'ytick',[OFFSET nsub],'yticklabel',{sprintf('%d',round(OFFSET/nsub*100)),'100'});
                    end
                    
                end
                
                
                if isnan(vout(i)) && numel(line_px)>1 && inpatch && i ~=numel(vout)
                    line_px(end+1) = (i-0.5)*2.2;
                    line_py(end+1) = OFFSET;
                    line_px(end+1) = line_px(1);
                    line_py(end+1) = OFFSET;
                    draw_it;
                    inpatch=0;
                    % keyboard;
                end
                
                if i==numel(vout) && inpatch
                    line_px(end+1) = (i+0.5)*2.2;
                    line_py(end+1) = OFFSET;
                    line_px(end+1) = line_px(1);
                    line_py(end+1) = OFFSET;
                    draw_it;
                    inpatch=0;
                end
                
                
                time_which_consistency = [ (1:(numel(vout) + floor(time_lag/2)))'*2.2 [[zeros(floor(time_lag/2), 1); outc(:,1); zeros(floor(time_lag/2),1)] [vout; zeros(floor(time_lag/2), 1)]]];
                new_statemat(1).(my_cond_names{i_gr}) = time_which_consistency;
                
            end
            
            if i_gr2 == 1
                % keyboard;
            end
            
            
            
            %
            %
            % END EXTRA FIG -- STATES (MOST COMMON)
            %
            %
            %
            
            
            
            %
            %
            %
            %
            %  Insert the ANNOTATIONS!
            %
            %
            %
            %
            
            anno_mat = load([base_cwd '/../data/annotations/annotations_matrix.txt']);
            
            this_anno_scan=analyses{i_analyses,1}{i_gr};
            do_anno=0;
            switch this_anno_scan
                case {'mov1a','mov1b'}
                    do_anno=1;
                    anno_vols=6:535;
                otherwise
                    do_anno=0;
            end
            bigmat_labels={'language','changepoint','faces   +','-','scenes   +','-'};
            
            print_anno_text=0;
            if i_gr==1
                print_anno_text=1;
            end
            
            if do_anno
                
                for i_bigmat=[1:6]
                    annov=anno_mat(:,i_bigmat);
                    annov=annov(anno_vols);
                    
                    % newah(end+1) = axes('position',[0.2 1-1/9*(i_bigmat) 0.6 1/10]);
                    ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-0.05-0.1-0.02-0.02*i_bigmat new_xs 0.015]);
                    
                    imagesc(annov');
                    set(ah(end),'visible','off','clim',[0 1]);
                    if print_anno_text
                        th=text(min(get(gca,'xlim'))-0.05*diff(get(gca,'xlim')), mean(get(gca,'ylim')),bigmat_labels{i_bigmat});
                        set(th,'fontsize',9,'horizontalalignment','right');
                    end
                    
                    all_anno_axes(end+1) = ah(end);
                end
                
                
                
            end
            
            
            
            
            
            
            
            
        end
        
        ah(end+1) = axes('parent',fh,'position',[new_xp + new_xs+0.02, bigspacingyu, bigspacingxu, 1-bigspacingyl-bigspacingyu]);
        
        set(ah(end),'visible','off');
        cb=colorbar('West');
        set(cb,'Limits',[0.5 0.5+NSTATES]);
        set(ah(end),'clim',[0.5 0.5+NSTATES]);
        cb.Position=cb.Position + [0 0 0.0018 0];
        cb.Ticks=[1:NSTATES];
        cb.Label.String='HMM States';
        %set(cb,'colormap',cmat);
        
        
        
        colormap(cmat);
        set(all_anno_axes,'colormap',parulamap);
        
        calc_size=0;
        for imat=1:numel(mat)
            calc_size = calc_size+size(mat{imat},1);
        end
        calc_size= calc_size/50;
        calc_size= calc_size + (numel(mat)-1) *0.8;
        
        
        set(fh,'paperunits','centimeters');
        set(fh,'papersize',[10+ calc_size, 20]*1.2);
        set(fh,'paperposition',[0 0 10+ calc_size 20]*1.2);
        
        
        % fig_fname=['Fig2_vpath_' preprocessing_dir '-' analysis '-rep-' num2str(SUMMARY_MEASURES_FILE) '.jpg'];
        
        output_filename = [scriptpath filesep '..' filesep 'figures' filesep 'Fig2_' preprocessing_dir '-' analysis_dir '-' num2str(SUMMARY_MEASURES_FILE) '.jpg'];
        
        
        print('-djpeg','-r600', output_filename);

        
        
        
        
        % Matlab tables are definitely not pandas dataframes:
        Participant = {};
        for i=1:numel(subs_to_use)
            Participant{end+1} = sprintf('Participant_%d',i); 
        end
        Participant=Participant';
        Participant{end+1} = 'Consistency';
        Participant{end+1} = 'MostConsistentState';
        
        
        % now we have made several variables; so we need to call Table with
        % that; i need to say myTable, otherwise it might interfere with
        % another Table already in the namespace. Furthermore, we need to
        % use eval, since Table names columns according to names in the
        % workspace.
        
        % Making Tables in Matlab. It's so simple and elegant, a true joy
        % to behold, and an example to strive towards to. So much better
        % than Pandas dataframes!

        
        %         Circumvent the following error:
        %         Error using writetable (line 124)
        % The data block starting at cell 'A1' exceeds the sheet boundaries by 0 row(s) and 1235 column(s).
        % 
        % Error in Make_Fig2_state_path_figure (line 648)
        %         writetable(myTable, output_filename_values)

        

        for i=1:numel(my_cond_names)
            s='';
            s=[s sprintf('%s = mat{%d}\''\n',my_cond_names{i}, i)];
            s=[s sprintf('%s=[%s; round(new_statemat.%s(:,2)/numel(subs_to_use)*100)\'']\n',my_cond_names{i},my_cond_names{i},my_cond_names{i})];
            s=[s sprintf('%s=[%s; new_statemat.%s(:,3)\'']\n',my_cond_names{i},my_cond_names{i},my_cond_names{i})];
            
            s=[s sprintf('%s = num2str(%s)', my_cond_names{i}, my_cond_names{i})];
            
            eval(s);
            
        end


        
            

        s='';
        s=[s sprintf('myTable = table(')];
        s=[s sprintf('Participant,')];
        my_conditions = fieldnames(new_statemat);
        for i=1:numel(my_conditions)
            s=[s sprintf('%s,',my_conditions{i})];
            
        end
        s(end)=[]; % remove comma
        s=[s sprintf(');')];
        eval(s); % this will make the table.
        % fsource=fig_fname;
        
        
        output_filename_values = [scriptpath filesep '..' filesep 'figures' filesep 'Fig2_' preprocessing_dir '-' analysis_dir '-' num2str(SUMMARY_MEASURES_FILE) '_values.xls'];
        
        myTable;
        writetable(myTable, output_filename_values)
        
        
        
        if CLOSE_FIGS
            close(fh);
        end
        
        
        cd(base_cwd);
        
        
        
        
    end
end

% save this information in a matlab struct to be used later on.
save sd.mat sd




