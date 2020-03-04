% process the eye data

% after inspection I foudn a couple of file which SEEM to have OK-kind-of
% eyetracking data:

files_with_good_data = {...
'eye_sub_002_mov1_b.mat';...
'eye_sub_003_mov1_b.mat';...
'eye_sub_006_mov1_a.mat';...
'eye_sub_006_mov2_a.mat';...
'eye_sub_007_mov1_b.mat';...
'eye_sub_007_rest_b.mat';...
'eye_sub_008_mov1_a.mat';...
'eye_sub_008_mov1_b.mat';...
'eye_sub_008_mov2_a.mat';...
'eye_sub_009_mov1_b.mat';...
'eye_sub_009_rest_b.mat';...
'eye_sub_010_mov1_b.mat';...
'eye_sub_011_mov1_b.mat';...
'eye_sub_011_rest_b.mat';...
'eye_sub_012_mov1_b.mat';...
'eye_sub_012_rest_b.mat';...
'eye_sub_013_mov1_b.mat';...
'eye_sub_014_mov1_b.mat';...
'eye_sub_016_mov1_a.mat';...
'eye_sub_016_mov2_a.mat';...
'eye_sub_016_rest_a.mat';...
'eye_sub_016_rest_b.mat';...
'eye_sub_017_mov1_a.mat';...
'eye_sub_017_mov2_a.mat';...
'eye_sub_017_rest_a.mat';...
'eye_sub_018_mov1_a.mat';...
'eye_sub_018_mov2_a.mat';...
'eye_sub_020_mov1_a.mat';...
'eye_sub_020_mov1_b.mat';...
'eye_sub_020_mov2_a.mat';...
'eye_sub_020_rest_a.mat';...
'eye_sub_020_rest_b.mat';...
'eye_sub_021_mov1_a.mat';...
'eye_sub_021_mov2_a.mat';...
'eye_sub_021_rest_a.mat';...
};



TR=2.2;
fs=500;

lines={};
while ~feof(fid)
    lines{end+1} = fgetl(fid);
end

s=struct();
notworking={};
for i_f=1:numel(lines)
    
    try
        
        fname = lines{i_f};
        
        if strcmp(fname,'eye_sub_003_mov1_b.mat')
            fname = 'eye_sub_003_mov2_a.mat';
            warning('switch!');
            % wrongly put into wrong directory
        end
        
        % alternate starting time:
        % edf1.edf1.Events.Start.time
        
        
        submatch=regexprep(fname,'.*_([0-9]{3})_(.{4})_([ab]{1}).mat','$1');
        taskmatch=regexprep(fname,'.*_([0-9]{3})_(.{4})_([ab]{1}).mat','$2');
        sessmatch=regexprep(fname,'.*_([0-9]{3})_(.{4})_([ab]{1}).mat','$3');
        
        switch taskmatch
            case 'rest'
                nvol=220;
            case 'mov1'
                nvol=535;
            case 'mov2'
                nvol=250;
        end
        
        
        fbase=fname(1:end-4);
        
        edf1=load(fname);
        
        % find beginning marker?
        % we assume first value of 5 == the MRI trigger.
        % we assume start with ==0;
        
        if edf1.edf1.Events.Input.value(1)==255
            % edf1.edf1.Events.Input.value(1)=[];
        end
        
        if edf1.edf1.Events.Input.value(1) ~= 0
            warning('logfile doesnt have 0 as first input signal?')
            
        end
        
        if sum(edf1.edf1.Events.Input.value==5) < nvol-50
            error('something wrong with the triggers');
            continue
        end
        
        first_mri_trig_ind = find([edf1.edf1.Events.Input.value]==5,1,'first');
        first_mri_trig_time = edf1.edf1.Events.Input.time(first_mri_trig_ind);
        
        start_acquisition_time = first_mri_trig_time;
        stop_acquisition_time = first_mri_trig_time + 2.2*nvol*1000;
        
        
        sample_eyemov_signal_start = find(abs(edf1.edf1.Samples.time-start_acquisition_time>0),1,'first');
        
        if numel(sample_eyemov_signal_start)==0
            sample_eyemov_signal_start = find(abs(edf1.edf1.Samples.time-edf1.edf1.Events.Start.time>0),1,'first');
        end
        
        
        
        sample_eyemov_signal_stop = find(abs(edf1.edf1.Samples.time-stop_acquisition_time>0),1,'first');
        if numel(sample_eyemov_signal_stop)==0
            sample_eyemov_signal_stop = sample_eyemov_signal_start + round(nvol*TR*fs);
        end
        
        
        spupil = edf1.edf1.Samples.pupilSize;
        tpupil = edf1.edf1.Samples.time;
        
        
        
        newv=spupil;
        
        
        % deal with eyeblinks by interpolation:
        % fund NaN's..
        nanvalues = newv< 100;
        
        newv(nanvalues)=NaN;
        
        % expand upon that (eyeblinks have edge artifact stuff:
        marked = find(nanvalues);
        old_marked=marked;
        new_list=[];for i=1:numel(marked);for j=-55:55;new_list(end+1) = marked(i)+j;end;end
        new_list=unique(new_list);
        new_list(new_list<1)=[];
        new_list(new_list>numel(newv))=[];
        newv(new_list)=NaN;
        
        
        
        % getting rid of spikes
        
        
        
        % interpolate:
        sample_points=1:numel(newv);
        [F,TF] = fillmissing(newv,'linear','SamplePoints',sample_points);
        
        F_despiked2 = medfilt1(F, 25);
        
        figure;plot(spupil);
        hold on;
        plot(new_list,1000*ones(1,length(new_list)),'kx');
        plot(old_marked,1000*ones(1,length(old_marked)),'gx');
        ylim([0 5000]);
        plot(F);
        plot(F_despiked2,'m-');
        title(fname,'interpreter','none');
        
        line(sample_eyemov_signal_start*[1 1],get(gca,'ylim'),'color','r','linewidth',2);
        line(sample_eyemov_signal_stop*[1 1],get(gca,'ylim'),'color','r','linewidth',2);
        
        
        
      %  switch fname
         %   case    {...
                    % 'eye_sub_003_mov1_b.mat',...
                    % 'eye_sub_003_rest_b.mat',...
                    % 'eye_sub_008_mov1_a.mat',...
                    % 'eye_sub_008_rest_b.mat',... --> is bad.
                    % 'eye_sub_007_rest_b.mat',...
                    % 'eye_sub_007_mov1_b.mat',...
                    % 'eye_sub_008_mov2_a.mat',...
                    % 'eye_sub_010_mov1_b.mat',...
                    % 'eye_sub_011_rest_b.mat',...
                    % 'eye_sub_012_rest_b.mat',...
                    % 'eye_sub_013_mov1_b.mat',... 
                    % 'eye_sub_020_rest_b.mat',...
                    % 'eye_sub_016_rest_a.mat',...
                    % 'eye_sub_016_rest_b.mat',...
                  %  }
                % keyboard;
     %   end
        
        saveas(gcf,['tr-' taskmatch '-' sessmatch '-' fbase '.fig']);
        close(gcf);
        
        
        % keyboard;
        
        
        
        sely=F_despiked2(sample_eyemov_signal_start:sample_eyemov_signal_stop);
        selt=tpupil(sample_eyemov_signal_start:sample_eyemov_signal_stop);
        
        
        out_vec = [];
        
        for j=1:nvol
            b=round((j-1)*TR*fs)+1;
            e=b + round(TR*fs)-1;
            out_vec(end+1) = mean(sely(b:e));
            % line(b*[1 1],get(gca,'ylim'));
        end
        
        
        fh=figure;
        plot(out_vec);
        title(sprintf('nvol: %d,  %s',nvol,fname),'interpreter','none');
        xlim([0 nvol+1]);
        saveas(gcf,['tr-' taskmatch '-' sessmatch '-' fbase '.jpg']);
        close(gcf);
        
    catch
        
        notworking{end+1} = fname;
    end
    
    
    s.([taskmatch sessmatch]).(['s' submatch]) = out_vec;
    
end


save s.mat s

% t=edf1.Samples.time;
% figure;plot(t,v);
% real_b_time = edf1.Events.Input.time(find(edf1.Events.Input.value==5,1,'first'));
% nvol=535;
% tr=2.2;
% real_end_time = real_b_time + nvol*tr*1000;
% line(real_b_time*[1 1],get(gca,'ylim'),'color','g')
% line(real_end_time*[1 1],get(gca,'ylim'),'color','g')
% samples = intersect(find(t>real_b_time), find(t<real_end_time));
% eyedata=v(samples);
% nanvalues = eyedata==0;
% figure;plot(eyedata)
% nanvalues
% find(nanvalues)
% plot(find(nanvalues),'rx')
% figure;plot(eyedata)
% hold on;
% plot(find(nanvalues),1000*ones(1,length(find(nanvalues))),'rx');
% nanvalues = eyedata< 300;
% nanvalues = eyedata< 100;
% plot(find(nanvalues),1000*ones(1,length(find(nanvalues))),'gx');
% find(nanvalues)
% marked = find(nanvalues);
% new_list=[];for i=1:numel(marked);for j=-10:10;new_list(end+1) = marked(i)+j;end;end
% new_list=unique(new_list);
% plot(new_list,1000*ones(1,length(new_list)),'kx');
% new_list=[];for i=1:numel(marked);for j=-15:15;new_list(end+1) = marked(i)+j;end;end
% help fillmissing
% newv=eyedata;
% newv(new_list)=NaN;
% newv(new_list)
% numel(newv)
% new_list
% new_list=unique(new_list);
% newv(new_list)=NaN;
% new_list(1)
% new_list(end)
% new_list(new_list<1)=[];
% new_list(new_list>numel(newv))=[];
% newv(new_list)=NaN;
% newnewv=fillmissing(newv);
% doc fillmissing
% sample_points=1:numel(newv);
% [F,TF] = fillmissing(newv,'linear','SamplePoints',sample_points);
% figure;plot(F)
% plot(F,'r')
% ls
% cd
% ls
% %-- 04/02/19 10:26:20 AM --%
% addpath(genpath('PhysioNet-Cardiovascular-Signal-Toolbox-1.0.1/
% addpath(genpath('PhysioNet-Cardiovascular-Signal-Toolbox-1.0.1/'));
% cd PhysioNet-Cardiovascular-Signal-Toolbox-1.0.1/
% ls
% cd Demos/
% ls
% Basic_Demo
% cd ..
% ls
% cd ..
% ls
% Basic_Demo
% cd PhysioNet-Cardiovascular-Signal-Toolbox-1.0.1/
% ls
% InitializeHRVparams
% Main_HRV_Analysis
% ver
% ls
% cd process-physiology/
% ls
% cd eye/
% ls
% ls *.m
% edit plotting.m
% cd edf-converter/
% ls
% cd ..
% ls
% pwd
% ls

