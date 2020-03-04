% only the traces, not all of the symmary -0 that's fig 3.



cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/

addpath /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/


addpath(genpath(pwd));
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

parulamap=colormap('parula');

cmat=cmat(1:NSTATES,:);
cmat=cmat./255;


% container for all data -- PER SUBJECT
sd=struct(); % subject data

my_cond_names = {'mova','movb','resta','restb'};
new_statemat = struct(); new_statemat(:)=[]; % this is for counting/keeping track of the states.


preprocessings={'normal','aroma','aroma-gsr'};


for i_summarymeasures=13 % 1:15
    SUMMARY_MEASURES_FILE=i_summarymeasures;
    
    for i_preprocessing=2 %1:numel(preprocessings)
        preprocessing_dir=preprocessings{i_preprocessing};
        
        
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
        
        
        
        for i_analyses=9 % 1:size(analyses,1)
            
            
            all_anno_axes=[];
            fh=figure('color','w','visible','off');
            set(fh,'position',[50         120        1700         800]);
            if ~CLOSE_FIGS
                set(fh,'visible','on');
            end
            
            scan_names = analyses{i_analyses, 1};
            
            analysis = analyses{i_analyses, 2}(1:end-4);
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
            
            cd(['../results_new2_' num2str(NSTATES)]);
            cd(preprocessing_dir)
            cd(analysis_dir)
            
            
            load(['Summary_measures_rep_' num2str(SUMMARY_MEASURES_FILE) '.mat']); % this will get our vpath.
            NSTATES=10;
            
            switch analysis_dir
                
                case 'mov1a---2'
                    nvols=[535];
                case 'mov2a---2'
                    nvols=[250];
                case 'mov1b---2'
                    nvols=[535];
                case 'mov1a-1b---2'
                    nvols=[535 535];
                case 'mov1a-2a---2'
                    nvols=[535 250];
                case 'resta---2'
                    nvols=[220];
                case 'restb---2'
                    nvols=[220];
                case 'resta-b---2'
                    nvols=[220 220];
                case 'all---2'
                    nvols=[535 535 220 220];
                case 'allPlus---2'
                    nvols=[535 250 535 220 220];
                case 'sessiona---2'
                    nvols = [535 220];
            end
            nvols=nvols-SCANS_NULLED;
            % edit pr
            
            
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
            
            
            
            nsubs=21;
            subs_to_use=1:21;
            for i=1:numel(scan_names)
                subs_to_use = intersect(subs_to_use, slist.(scan_names{i}));
                nsubs=numel(subs_to_use);
            end
            
            
            
            % size(vpath);
            nscans=sum(nvols);
            nsub = numel(vpath)/nscans;
            
            % this should partition the vpath into subsections... if all is well.
            % but in our case, we should use some more smart code for it.
            
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
            
            % group-by-scan...
            
            
            
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
            
            
            
            new_order=[3 1 4 2];
            
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
                    
                    
                    % convert state indices into ehm.. a 0-1 matrix
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
                        th=text(ah(end),-10,0,sprintf('S%d - %d',chosen_subs(i),''));
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
%                 
%                 %
%                 %
%                 % END EXTRA FIG -- STATES (MOST COMMON)
%                 %
%                 %
%                 %
%                 
%                 
%                 %
%                 %
%                 % EXTRA FIG -- HEART RATE!
%                 %
%                 %
%                 %
%                 MAT_DO_LOWER=1;
%                 LEAVEAXIS=0;
%                 
%                 % do HR:
%                 ce=load('/home/johan/mnt/hpcworking/projects/hmm-movie/process-physiology/tapas_io_toolbox/s_electrophys.mat');
%                 
%                 curr_analysis = analyses{i_analyses,1}{i_gr};
%                 sel=ce.s.hr.(curr_analysis);
%                 
%                 fns=fieldnames(sel);
%                 
%                 the_matrix = [];
%                 sel_subs=[];
%                 for ifns=1:numel(fns)
%                     if sum(str2double(fns{ifns}(2:end))==slist.(curr_analysis))==1
%                         the_matrix(:,end+1) = sel.(fns{ifns});
%                         sel_subs(end+1) = str2double(fns{ifns}(2:end));
%                     end
%                     
%                 end
%                 
%                 % temporal derivative == HRV?
%                 % the_matrix=diff(the_matrix);
%                 % the_matrix=[zeros(1,size(the_matrix,2)); the_matrix];
%                 % ... DO subtract the av HR from each guy:
%                 the_matrix = the_matrix - ones(size(the_matrix,1),1)*mean(the_matrix,1);
%                 
%                 the_matrix(1:5,:)=[];
%                 old_the_matrix=the_matrix;
%                 
%                 ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-0.05-0.05-0.08-0.02+0.03+0.015+0.005 new_xs 0.07-0.02-0.005-0.005]);
%                 
%                 MAT_DO_LOWER=1;
%                 LEAVEAXIS=0;
%                 SCALING_MATRIX=0.5;
%                 QLTHR=0.2;
%                 QUTHR=0.8;
%                 MAT_AX_LABEL='HR';
%                 plot_the_matrix_down
%                 
%                 
%                 % i will be able to safely use here, since we're not doing
%                 % it for all the subjects.
%                 % to_use_fns=regexprep(fns,'s0*','s');
%                 for i=1:numel(sel_subs)
%                     % keyboard;
%                     % individual HR data:
%                     sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(sel_subs(i))]).(analyses{i_analyses,1}{i_gr}).hrdata=old_the_matrix(:,i);
%                     
%                     outmat=zeros(size(old_the_matrix,1),2);
%                     outmat(old_the_matrix(:,i)<quantile(old_the_matrix(:,i),QLTHR),1) = 1;
%                     outmat(old_the_matrix(:,i)>quantile(old_the_matrix(:,i),QUTHR),2) = 1;
%                     
%                     % convert into 0's and 1's...
%                     sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(sel_subs(i))]).(analyses{i_analyses,1}{i_gr}).hrmat=outmat;
%                 end
%                 
%                 
%                 
%                 
%                 
%                 
%                 
%                 % HRV!
%                 ce=load('/home/johan/mnt/hpcworking/projects/hmm-movie/process-physiology/tapas_io_toolbox/s_electrophys.mat');
%                 
%                 curr_analysis = analyses{i_analyses,1}{i_gr};
%                 sel=ce.s.hr.(curr_analysis);
%                 
%                 fns=fieldnames(sel);
%                 
%                 the_matrix = [];
%                 sel_subs=[];
%                 for ifns=1:numel(fns)
%                     if sum(str2double(fns{ifns}(2:end))==slist.(curr_analysis))==1
%                         the_matrix(:,end+1) = sel.(fns{ifns});
%                         sel_subs(end+1) = str2double(fns{ifns}(2:end));
%                     end
%                     
%                 end
%                 
%                 % temporal derivative == HRV?
%                 the_matrix=diff(the_matrix);
%                 the_matrix=[zeros(1,size(the_matrix,2)); the_matrix];
%                 the_matrix(1:5,:)=[];
%                 ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-0.05-0.05-0.08-0.07-0.02+0.04+0.04 new_xs 0.07-0.02-0.005-0.005]);
%                 
%                 
%                 MAT_DO_LOWER=1;
%                 LEAVEAXIS=0;
%                 SCALING_MATRIX=0.5;
%                 QLTHR=0.2;
%                 QUTHR=0.8;
%                 MAT_AX_LABEL='HRV';
%                 plot_the_matrix_down
%                 
%                 
%                 % i will be able to safely use here, since we're not doing
%                 % it for all the subjects.
%                 % to_use_fns=regexprep(fns,'s0*','s');
%                 for i=1:numel(sel_subs)
%                     
%                     % individual HR data:
%                     sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(sel_subs(i))]).(analyses{i_analyses,1}{i_gr}).hrvdata=the_matrix(:,i);
%                     
%                     outmat=zeros(size(the_matrix,1),2);
%                     outmat(the_matrix(:,i)<quantile(the_matrix(:,i),QLTHR),1) = 1;
%                     outmat(the_matrix(:,i)>quantile(the_matrix(:,i),QUTHR),2) = 1;
%                     
%                     % convert into 0's and 1's...
%                     sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(sel_subs(i))]).(analyses{i_analyses,1}{i_gr}).hrvmat=outmat;
%                 end
%                 
%                 
%                 
%                 
%                 % the RVT, too.
%                 ce=load('/home/johan/mnt/hpcworking/projects/hmm-movie/process-physiology/tapas_io_toolbox/s_electrophys.mat');
%                 
%                 
%                 curr_analysis = analyses{i_analyses,1}{i_gr};
%                 sel=ce.s.rvt.(curr_analysis);
%                 
%                 fns=fieldnames(sel);
%                 
%                 the_matrix = [];
%                 sel_subs=[];
%                 for ifns=1:numel(fns)
%                     if sum(str2double(fns{ifns}(3:end))==slist.(curr_analysis))==1
%                         the_matrix(:,end+1) = sel.(fns{ifns});
%                         sel_subs(end+1) = str2double(fns{ifns}(2:end));
%                     end
%                     
%                 end
%                 % the_matrix = log(abs(the_matrix)); % log-transform of this...
%                 the_matrix = the_matrix - ones(size(the_matrix,1),1)*mean(the_matrix,1);
%                 the_matrix(1:5,:)=[];
%                 
%                 
%                 
%                 ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-0.05-0.05-0.05-0.07-0.1-0.02+0.05+0.060 new_xs 0.07-0.02-0.01]);
%                 
%                 MAT_DO_LOWER=1;
%                 LEAVEAXIS=0;
%                 SCALING_MATRIX=0.5;
%                 QLTHR=0.2;
%                 QUTHR=0.8;
%                 MAT_AX_LABEL='RVT';
%                 
%                 plot_the_matrix_down
%                 
%                 % i will be able to safely use here, since we're not doing
%                 % it for all the subjects.
%                 % to_use_fns=regexprep(fns,'s0*','s');
%                 for i=1:numel(sel_subs)
%                     
%                     % individual HR data:
%                     sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(sel_subs(i))]).(analyses{i_analyses,1}{i_gr}).rvtdata=the_matrix(:,i);
%                     
%                     outmat=zeros(size(the_matrix,1),2);
%                     outmat(the_matrix(:,i)<quantile(the_matrix(:,i),QLTHR),1) = 1;
%                     outmat(the_matrix(:,i)>quantile(the_matrix(:,i),QUTHR),2) = 1;
%                     
%                     % convert into 0's and 1's...
%                     sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(sel_subs(i))]).(analyses{i_analyses,1}{i_gr}).rvtmat=outmat;
%                 end
%                 
%                 
%                 
%                 
%                 %
%                 %
%                 % END EXTRA FIG -- PUPIL DIAMETER
%                 %
%                 %
%                 %
%                 
%                 ce=load('/home/johan/mnt/hpcworking/projects/hmm-movie/process-physiology/eye/s.mat');
%                 
%                 curr_analysis = analyses{i_analyses,1}{i_gr};
%                 sel=ce.s.(curr_analysis);
%                 
%                 fns=fieldnames(sel);
%                 
%                 the_matrix = [];
%                 sel_subs=[];
%                 for ifns=1:numel(fns)
%                     if sum(str2double(fns{ifns}(2:end))==slist.(curr_analysis))==1
%                         the_matrix(:,end+1) = sel.(fns{ifns});
%                         sel_subs(end+1) = str2double(fns{ifns}(2:end));
%                     end
%                     
%                 end
%                 the_matrix(1:5,:)=[];
%                 
%                 
%                 
%                 ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-0.05-0.05-0.05-0.07-0.1-0.07-0.02+0.06+0.08-0.005 new_xs 0.07-0.02-0.01]);
%                 
%                 MAT_DO_LOWER=1;
%                 LEAVEAXIS=0;
%                 SCALING_MATRIX=0.5;
%                 QLTHR=0.2;
%                 QUTHR=0.8;
%                 MAT_AX_LABEL='Pupil';
%                 plot_the_matrix_down
%                 
%                 
%                 % i will be able to safely use here, since we're not doing
%                 % it for all the subjects.
%                 % to_use_fns=regexprep(fns,'s0*','s');
%                 for i=1:numel(sel_subs)
%                     
%                     % individual HR data:
%                     sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(sel_subs(i))]).(analyses{i_analyses,1}{i_gr}).eyedata=the_matrix(:,i);
%                     
%                     outmat=zeros(size(the_matrix,1),2);
%                     outmat(the_matrix(:,i)<quantile(the_matrix(:,i),QLTHR),1) = 1;
%                     outmat(the_matrix(:,i)>quantile(the_matrix(:,i),QUTHR),2) = 1;
%                     
%                     % convert into 0's and 1's...
%                     sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(sel_subs(i))]).(analyses{i_analyses,1}{i_gr}).eyemat=outmat;
%                 end
%                 
%                 
%                 
%                 
%                 ce=load('/home/johan/mnt/hpcworking/projects/hmm-movie/process-physiology/gsr/for_ledalab/exports/s_gsr.mat');
%                 
%                 
%                 
%                 
%                 curr_analysis = analyses{i_analyses,1}{i_gr};
%                 sel=ce.s.(curr_analysis);
%                 
%                 fns=fieldnames(sel);
%                 
%                 the_matrix = [];
%                 sel_subs=[];
%                 for ifns=1:numel(fns)
%                     if sum(str2double(fns{ifns}(3:end))==slist.(curr_analysis))==1
%                         the_matrix(:,end+1) = sel.(fns{ifns});
%                         sel_subs(end+1) = str2double(fns{ifns}(2:end));
%                     end
%                     
%                 end
%                 % the_matrix = log(abs(the_matrix)); % log-transform of this...
%                 the_matrix(1:5,:)=[];
%                 
%                 ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-0.05-0.05-0.05-0.07-0.1-0.07-0.02+0.06+0.08-0.045 new_xs 0.07-0.02-0.01]);
%                 
%                 
%                 MAT_DO_LOWER=0;
%                 LEAVEAXIS=1;
%                 SCALING_MATRIX=0.5;
%                 QLTHR=0.2;
%                 QUTHR=0.8;
%                 MAT_AX_LABEL='GSR';
%                 
%                 plot_the_matrix_down
%                 
%                 
%                 % i will be able to safely use here, since we're not doing
%                 % it for all the subjects.
%                 % to_use_fns=regexprep(fns,'s0*','s');
%                 for i=1:numel(sel_subs)
%                     
%                     % individual HR data:
%                     sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(sel_subs(i))]).(analyses{i_analyses,1}{i_gr}).gsrdata=the_matrix(:,i);
%                     
%                     outmat=zeros(size(the_matrix,1),2);
%                     % outmat(the_matrix(:,i)<quantile(the_matrix(:,i),QLTHR),1) = 1;
%                     outmat(the_matrix(:,i)>quantile(the_matrix(:,i),QUTHR)) = 1;
%                     
%                     % convert into 0's and 1's...
%                     sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(sel_subs(i))]).(analyses{i_analyses,1}{i_gr}).gsrmat=outmat;
%                 end
%                 
%                 
%                 
%                 
%                 %
%                 %
%                 %  Insert the ANNOTATIONS!
%                 %
%                 %
%                 %
%                 %
%                 
%                 anno_mat = load('/home/johan/mnt/hpcworking/projects/hmm-movie/rawdata/1_Data/Annotation/annobigmat.mat');
%                 bigmat=anno_mat.bigmat;
%                 
%                 % now -- check the situation... movie or not?
%                 
%                 % check which task we're in: mov1a, mov1b -- do full. ---
%                 % mov2a, cut to last 250 volumes.
%                 this_anno_scan=analyses{i_analyses,1}{i_gr};
%                 do_anno=0;
%                 switch this_anno_scan
%                     case {'mov1a','mov1b'}
%                         do_anno=1;
%                         anno_vols=6:535;
%                     case 'mov2a'
%                         do_anno=1;
%                         anno_vols=291:535;
%                     otherwise
%                         do_anno=0;
%                 end
%                 bigmat_labels={'language','changepoint','plottwist','faces   +','n','-','scenes   +','n','-'};
%                 
%                 print_anno_text=0;
%                 if i_gr==1
%                     print_anno_text=1;
%                 end
%                 
%                 if do_anno
%                     
%                     for i_bigmat=1:9
%                         annov=bigmat(:,i_bigmat);
%                         annov=annov(anno_vols);
%                         
%                         % newah(end+1) = axes('position',[0.2 1-1/9*(i_bigmat) 0.6 1/10]);
%                         ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-0.05-0.05-0.05-0.07-0.1-0.07-0.02+0.06+0.08-0.045-0.02-0.02*i_bigmat new_xs 0.015]);
%                         
%                         imagesc(annov');
%                         set(ah(end),'visible','off','clim',[0 1]);
%                         if print_anno_text
%                             th=text(min(get(gca,'xlim'))-0.05*diff(get(gca,'xlim')), mean(get(gca,'ylim')),bigmat_labels{i_bigmat});
%                             set(th,'fontsize',9,'horizontalalignment','right');
%                         end
%                         
%                         all_anno_axes(end+1) = ah(end);
%                     end
%                     
%                     
%                     
%                 end
%                 
                
                
                
                
                
                
                
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
            
            
            
            fig_fname=['vpath_manuscript-' preprocessing_dir '-' analysis '-rep-' num2str(SUMMARY_MEASURES_FILE) '.jpg'];
            saveas(fh,fig_fname);
            
            fsource=fig_fname;
            
            ftarget_file = 'f2_vpath_sup4.jpg';
            ftarget1=['/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript/output/' ftarget_file];
            ftarget2=['/home/johan/ldrive/Lab_LucaC/10_Johan/1_HMMMovie/' ftarget_file];
            
            copyfile(fsource,ftarget1);
            copyfile(fsource,ftarget2);
            
            cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/figures_for_manuscript
            
            
            if CLOSE_FIGS
                close(fh);
            end
            
            
            % fh_statetime=fh;
            % cd /home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/
            
            
        end
        
    end
end


% save sd.mat sd




% crazy sorting fun, I focus on OFTEN OCCURRING states
% for i=1:size(mat,2);
%     uv=unique(mat(i,:));
%     if numel(uv)>1
%         cnt=[];
%         for j=1:numel(uv)
%             cnt(end+1) = sum(mat(i,:)==uv(j));
%         end
%         tmpm=sortrows([uv' cnt'],'descend');
%     end
% end
%             uv=unique(mat(:));
%             us=[];
%             for i=1:numel(uv)
%                 us(end+1) = sum(mat(:)==uv(i));
%             end
% wow - states 6, 7 and 10 do not occur. Their vibreny path isn'tlong
% enough.



% crazy sorting fun - are people similar??
%             cal_sub=1;
%             % nsubs=21;
%             start_sub_vec=2:nsub;
%
%             chosen_subs=[1];
%
%             for cntr=1:20
%                 dist_vals = [];
%
%                 for i=1:numel(start_sub_vec)
%                     % calculate the corrs...
%                     cal_v = mat(:,chosen_subs(end));
%                     cal_vec_2 = mat(:,start_sub_vec(i));
%                     dist_vals(end+1) = pdist([cal_v';cal_vec_2'],'cityblock');
%                 end
%                 chosen_sub = start_sub_vec(dist_vals==min(dist_vals));
%                 chosen_subs(end+1) = chosen_sub;
%
%                 start_sub_vec(start_sub_vec==chosen_sub)=[];
%             end
%
%
%             % let's try another sorting fun. clustering~
%             opts = statset('Display','final');
%             [idx,C] = kmeans(mat',8,'Distance','sqeuclidean',...
%                 'Replicates',5,'Options',opts);
%
%             chosen_subs_mat=sortrows([(1:nsub)' idx],2); % change as needed...
%             chosen_subs_mat=sortrows(chosen_subs_mat,1);
%             chosen_subs = chosen_subs_mat(:,1);
%             gr_assignment = chosen_subs_mat(:,2);

